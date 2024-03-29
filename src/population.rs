// author: André Amado

//! Population module

use crate::{
    matrix::Matrix,
    indices::Indices,
    fitness_landscape::FitnessLandscape
};

use rand::{
    prelude::*,
    distributions::WeightedIndex,
    seq::SliceRandom,
    self
};
use rand_distr::{Distribution, Binomial, Uniform};

use regex::Regex;

use std::{
    io::{BufReader, BufRead, Write},
    fs::File
};
type BufferedFile = std::io::BufWriter<std::fs::File>;

use clap::ArgMatches;

/// Stores information about allele fixation in a population
#[derive(Clone, Copy)]
pub enum Allele {
    /// Allele not fixed
    NotFixed,
    /// Allele `allele` fixed at generation `time`
    Fixed { time: usize, allele: bool },
}


/// Describes recombination maps
#[derive(Clone)]
pub enum Recombination {
    /// No recombination map
    None,
    /// Constant recombination map (same recombination probability everywhere), the field stores
    /// the recombination probability
    Constant(f64),
    /// Detailed recombination map, with a custom recombination probability at each locus, the
    /// field is a vector containing the full recombination map
    Map(Vec<f64>)
}

impl Recombination {
    /// Returns a summary of the recombination map in a String without whitespace
    pub fn summary(&self) -> String {
        let mut s = String::new();
        match self {
            Self::None        => s.push('0'),
            Self::Constant(r) => s.push_str(format!("{}", r).as_str()),
            Self::Map(map)    => {
                for r in map.iter() { s.push_str(format!(",{}", *r).as_str()) }
                s.remove(0);
            }
        }
        s
    }

    /// Reads and returns the recombination maps encoded in the user input
    ///
    /// # Arguments:
    /// * `matches` matches to user input provided by get_command_line_matches()
    /// * `l` number of loci in the fitness landscape
    pub fn from_args(matches: &ArgMatches, l: usize) -> Result<Vec<Self>, String> {
        let mut recombination_maps = Vec::<Recombination>::new();

        // Check if the option --recombination_rates is present
        if matches.is_present("recombination_rates") {
            for r in matches.get_many("recombination_rates").unwrap() {
                recombination_maps.push(Recombination::Constant(*r))
            }

        // If not, check if the option --recombination_map is present
        } else if matches.is_present("recombination_map") {
            let value = matches.get_one::<String>("recombination_map").unwrap().trim();

            // Check if the value is a sequence of numbers (that can be interpreted as a recombination map)
            if Regex::new(r"(\d+)(,\d+)*").unwrap().is_match(value) {
                let mut recombination_map = Vec::<f64>::new();
                for r in value.split(",") {
                    match r.parse::<f64>() {
                        Err(_) => return Err(format!("Error: could not parse '{}' as a number", r)),
                        Ok(r)  => recombination_map.push(r)
                    }
                }
                recombination_maps.push(Recombination::Map(recombination_map))

            // Otherwise, assume the value is a file name and tries to read the values in it
            } else {
                let file = match File::open(&value) {
                    Err(_) => return Err(format!("Error: file {} not found", &value)),
                    Ok(f)  => f
                };
                let lines = BufReader::new(file).lines();
                for line in lines {
                    let mut recombination_map = Vec::<f64>::new();
                    let line = line.unwrap();
                    let line = line.trim();
                    if line.starts_with("#") || line.is_empty() { continue }

                    for r in line.split(&[' ', ',', '\t', ';'][..]) {
                        match r.trim().parse::<f64>() {
                            Err(_) => return Err(format!("Error: could not parse '{}' as a number", r)),
                            Ok(r)  => recombination_map.push(r)
                        }
                    }
                    recombination_maps.push(Recombination::Map(recombination_map))
                }
            }

        // If none are present, assume there is no recombination
        } else {
            recombination_maps.push(Recombination::None)
        }

        // Validate recombination rates
        for map in recombination_maps.iter() {
            match map {
                Recombination::None => {},
                Recombination::Constant(r) => {
                    if *r < 0. || *r > 0.5 { return Err(format!("Error: invalid recombination rate (r = {}). Recombination rates should be contained in the interval [0, 0.5]", *r)) }
                },
                Recombination::Map(map) => {
                    // test if the recombination map has the correct number of recombination rates
                    if map.len() != l { return Err(format!("Error: recombination map should have {} positions (number of loci), {} found", l, map.len())) }

                    for r in map.iter() {
                        if *r < 0. || *r > 0.5 { return Err(format!("Error: invalid recombination rate (r = {}). Recombination rates should be contained in the interval [0, 0.5]", *r)) }
                    }
                }
            }
        }

        Ok(recombination_maps)
    }
}


/// Stores a population
pub struct Population {
    indices: Indices,
    pop: Vec<usize>,
    occupied_genotypes: Vec<usize>,
    recombination_map: Recombination,
    l: usize,
    size: usize,
    fixation: Vec<Allele>,
}

/// Describes an initial population
pub enum InitialPopulation {
    /// Initial population generated from a neutral site frequency spectrum
    NeutralSFS,
    /// Initial population generated from a neutral site frequency spectrum, with a drift threshold
    /// with the size given by the first field value and the second field probability that a minor
    /// allele is assigned a positive additive effect
    NeutralSFSDriftThreshold(usize, f64),
    /// Initial population generated with a uniform population, where the probaility of a derived
    /// allele is given by the field value
    Uniform(f64),
    /// Monomorphic with a specified genotype
    Monomorphic(Vec<bool>),
    RandomMonomorphic
}

impl InitialPopulation {
    /// Reads the initial population from the program arguments returning the corresponding
    /// InitialPopulation variant
    ///
    /// # Arguments:
    /// * `matches` matches to user input
    pub fn from_args(matches: &ArgMatches) -> Self {
        if matches.contains_id("neutralsfs") {
            if let Some(mut values) = matches.get_many("neutralsfs") {
                let d = if let Some(&d) = values.next() { d as usize } else { 1 };
                let f = if let Some(&f) = values.next() { f } else { 0.5 };

                Self::NeutralSFSDriftThreshold(d, f)
            } else {
                Self::NeutralSFS
            }
        } else if matches.contains_id("monomorphic") {
            let default = String::from("random");
            let arg = matches.get_one::<String>("monomorphic").unwrap_or(&default);

            if arg == "random" {
                Self::RandomMonomorphic
            } else {
                // Check if all characters are 0's and 1's
                if !arg.chars().all(|c| c == '0' || c == '1') {
                    panic!("Error: monomorphic option not recognized ({})", arg)
                }
                Self::Monomorphic(arg.chars().map(|g| g == '1').collect())
            }
        } else if matches.contains_id("uniform") {
            let p = *matches.get_one("uniform").unwrap();
            Self::Uniform(p)
        } else {
            panic!("No initial population found")
        }
    }

    pub fn short_description(&self) -> String {
        match self {
            Self::NeutralSFS => String::from("neutralsfs"),
            Self::NeutralSFSDriftThreshold(d, p) => format!("neutralsfs_dt{}_p{}", d, p),
            Self::Uniform(p) => format!("uniform({})", p),
            Self::Monomorphic(genotype) => {
                let sequence: String = 
                    genotype
                    .iter()
                    .map(|&g| if g {'1'} else {'0'})
                    .collect();
                format!("monomorphic({})", sequence)
            },
            Self::RandomMonomorphic => String::from("random_monomorphic"),
        }
    }
    pub fn shorter_description(&self) -> String {
        match self {
            Self::NeutralSFS => String::from("neutralsfs"),
            Self::NeutralSFSDriftThreshold(_, _) => String::from("neutralsfs"),
            Self::Uniform(_) => String::from("uniform"),
            Self::Monomorphic(_) => String::from("monomorphic"),
            Self::RandomMonomorphic => String::from("monomorphic"),
        }
    }
}

impl Population {
    /// Generates a new population
    ///
    /// # Arguments:
    /// * `l` number of loci
    /// * `size` population size
    /// * `recombination_map` recombination map
    pub fn new(l: usize, size: usize, recombination_map: Recombination) -> Population {
        Population {
            l,
            size,
            indices: Indices::new(l),
            pop: vec![0usize; 2usize.pow(l as u32)],
            occupied_genotypes: Vec::with_capacity(size),
            recombination_map,
            fixation: vec![Allele::NotFixed; l],
        }
    }

    /// Generates a new nonrecombining population
    ///
    /// # Arguments:
    /// * `l` number of loci
    /// * `size` population size
    pub fn without_recombination(l: usize, size: usize) -> Population {
        Population::new(l, size, Recombination::None)
    }

    /// Gets the number of loci
    pub fn get_l(&self) -> usize {
        self.l
    }

    /// Executes a recombination update
    pub fn recombination(&mut self) {
        // Check if a recombination map exists
        let recombination_map = match &self.recombination_map {
            Recombination::Map(map)    => map.clone(),
            Recombination::Constant(r) => vec![*r; self.l],
            Recombination::None        => return
        };

        let mut rng = thread_rng();
        let mut new_population = vec![0usize; self.pop.len()];
        let mut new_occupied_genotypes = Vec::with_capacity(self.size);

        // Generates a list with all the individuals in a random order
        let mut list_of_individuals = Vec::<usize>::with_capacity(self.size);
        for &i in self.occupied_genotypes.iter() {
            for _ in 0..self.pop[i] {
                list_of_individuals.push(i);
            }
        }
        list_of_individuals.shuffle(&mut rng);

        for indices in list_of_individuals.chunks(2) {
            // for odd populations
            if indices.len() == 1 {
                new_occupied_genotypes.push(indices[0]);
            }

            let mut idx1 = indices[0];
            let mut idx2 = indices[1];

            // if individuals are identical recombination skip trivial recombination
            if idx1 == idx2 {
                new_population[idx1] += 2;
                new_occupied_genotypes.push(idx1);
            } else {
                let recombine = (0..self.l).map(|i| -> bool {
                    let sample: f64 = (&mut rng).sample(rand::distributions::Open01);
                    sample < recombination_map[i]
                });

                let mut recombined = false;
                let mut genotype_1 = self.indices.to_sequence(idx1).to_vec();
                let mut genotype_2 = self.indices.to_sequence(idx2).to_vec();

                for (i, rec) in recombine.enumerate() {
                    if rec {
                        recombined = true;
                        let section1 = (&genotype_1[i..self.l]).to_vec();
                        let section2 = (&genotype_2[i..self.l]).to_vec();
                        (&mut genotype_1[i..self.l]).copy_from_slice(&section2[..]);
                        (&mut genotype_2[i..self.l]).copy_from_slice(&section1[..]);
                    }
                }

                if recombined {
                    idx1 = self.indices.to_index(&genotype_1);
                    idx2 = self.indices.to_index(&genotype_2);
                }
                new_population[idx1] += 1;
                new_population[idx2] += 1;
                new_occupied_genotypes.push(idx1);
                new_occupied_genotypes.push(idx2);
            }
        }

        //remove the repeated genotypes in the occupied_genotypes list
        new_occupied_genotypes.sort_unstable();
        new_occupied_genotypes.dedup();

        self.occupied_genotypes = new_occupied_genotypes;
        self.pop = new_population;
    }

    /// Executes a Wright-Fisher update
    ///
    /// # Arguments:
    /// * `landscape` fitness landscape
    pub fn wright_fisher(&mut self, landscape: &FitnessLandscape) {
        let mut rng = thread_rng();

        let fitness: Vec<f64> = self
            .occupied_genotypes
            .iter()
            .map(|&i| -> f64 { landscape.fitness_index(i) * (self.pop[i] as f64) })
            .collect();

        // Zero the old population
        for &i in self.occupied_genotypes.iter() {
            self.pop[i] = 0;
        }

        // Create the new population
        let new_indices = rand::distributions::WeightedIndex::new(&fitness).unwrap();
        for _ in 0..self.size {
            let genotype = self.occupied_genotypes[new_indices.sample(&mut rng)];
            self.pop[genotype] += 1;
        }

        self.filter_extinct_genotypes();
    }

    /// Executes a mutation update where each genotype can mutate at most once
    ///
    /// # Arguments:
    /// * `mutation_probability_per_locus`
    pub fn single_mutation(&mut self, mutation_probability_per_locus: f64) {
        let mut rng = thread_rng();

        let mutation_probability = mutation_probability_per_locus * self.l as f64;
        if mutation_probability > 1. {
            panic!("mutation probability per individual > 1 (p = {})!", mutation_probability)
        }

        let mut mutated_individuals: Vec<usize> = Vec::new();
        for &genotype_index in self.occupied_genotypes.iter() {
            let n_individuals = self.pop[genotype_index];

            let mutations_distribution = Binomial::new(n_individuals as u64, mutation_probability).unwrap();
            let n_mutations = mutations_distribution.sample(&mut rng) as usize;

            mutated_individuals.reserve(n_mutations);
            for _ in 0..n_mutations { mutated_individuals.push(genotype_index); }
            self.pop[genotype_index] -= n_mutations;
        };

        let locus_distribution: Uniform<usize> = Uniform::new(0, self.l);
        for &individual in mutated_individuals.iter() {
            let locus = locus_distribution.sample(&mut rng);

            let mut sequence: Vec<bool> = self.indices.to_sequence(individual).iter().map(|&i| i).collect();
            sequence[locus] = !sequence[locus];

            let index = self.indices.to_index(sequence.as_slice());
            if !self.occupied_genotypes.contains(&index) {
                self.occupied_genotypes.push(index);
            }
            self.pop[index] += 1;
        }

        self.filter_extinct_genotypes();
    }

    /// Executes a mutation update where each genotype can mutate several times
    ///
    /// # Arguments:
    /// * `mutation_probability_per_locus`
    pub fn multiple_mutation(&mut self, mutation_probability_per_locus: f64) {
        if mutation_probability_per_locus > 0. {
            let mut rng = thread_rng();

            let mut mutated_individuals: Vec<usize> = Vec::new();

            let n_mutations_individual = Binomial::new(self.l as u64, mutation_probability_per_locus).unwrap();
            for &genotype_index in self.occupied_genotypes.iter() {
                let n_individuals = self.pop[genotype_index];
                let original_sequence: Vec<bool> = self.indices.to_sequence(genotype_index).iter().map(|&i| i).collect();

                for _ in 0..n_individuals {
                    let n_mutations = n_mutations_individual.sample(&mut rng) as usize;

                    if n_mutations > 0 {
                        // remove the mutated individual from the population
                        self.pop[genotype_index] -= 1;

                        // determine the positions of the mutations
                        let mutation_positions = rand::seq::index::sample(&mut rng, self.l, n_mutations);
                        let mut sequence = original_sequence.clone();
                        for locus in mutation_positions.iter() {
                            sequence[locus] = !sequence[locus];
                        }

                        // compute the index of the new individual
                        let index = self.indices.to_index(sequence.as_slice());

                        // add the new individual to the list
                        mutated_individuals.push(index);
                    }
                }
            }

            // add the new individuals to the population
            for &individual in &mutated_individuals {
                if !self.occupied_genotypes.contains(&individual) {
                    self.occupied_genotypes.push(individual);
                }
                self.pop[individual] += 1;
            }

            self.filter_extinct_genotypes();
        }
    }

    /// Generates an initial population to replace the current population
    ///
    /// # Arguments:
    /// * `initial_population` type of initial population to generate
    pub fn generate_initial_population(&mut self, initial_population: &InitialPopulation) {
        let mut rng = thread_rng();

        {   // clear previous population
            self.occupied_genotypes.clear();
            let size = self.pop.len();
            self.pop.clear();
            self.pop.resize(size, 0);
        }

        match initial_population {
            InitialPopulation::NeutralSFS => {
                self.generate_initial_population(&InitialPopulation::NeutralSFSDriftThreshold(1, 0.5));
                return
            },
            &InitialPopulation::NeutralSFSDriftThreshold(k, f) => {
                // Generates an initial population where the alleles are distributing according
                // to the site frequency spectrum of a neutrally evolving population
                let mut initial_population = Matrix::<bool>::new(self.size, self.l);

                let harmonic_number: f64 = (k..=(self.size-k)).map(|j| 1. / (j as f64)).sum();
                let probability: Vec<f64> = (k..=(self.size-k))
                    .map(|j| 1. / (j as f64 * harmonic_number))
                    .collect();
                let dist = WeightedIndex::new(&probability).unwrap();

                for j in 0..self.l {
                    let order = rand::seq::index::sample(&mut rng, self.size, self.size);
                    let count = dist.sample(&mut rng) + k;

                    let has_positive_additive_effect = rng.gen_bool(f);
                    let minor_allele = if count < self.size/2 {
                        has_positive_additive_effect
                    } else {
                        !has_positive_additive_effect
                    };

                    for (c, i) in order.into_iter().enumerate() {
                        if c < count {
                            initial_population.set(i, j, minor_allele);
                        } else {
                            initial_population.set(i, j, !minor_allele);
                        }
                    }
                }

                for genotype in initial_population.row_iter() {
                    let idx = self.indices.to_index(&genotype);

                    self.pop[idx] += 1;
                    self.occupied_genotypes.push(idx);
                }
            },
            InitialPopulation::Uniform(derived_allele_probability) => {
                // Generates an initial population where each individual has a probability equal to
                // derived_allele_probability of carrying the derived allele form for each allele
                // Minor alleles are encoded have value true
                for _ in 0..self.size {
                    let genotype: Vec<bool> = (0..self.l)
                        .map(|_| (&mut rng).gen::<f64>() < *derived_allele_probability)
                        .collect();

                    let idx = self.indices.to_index(&genotype);

                    self.pop[idx] += 1;
                    self.occupied_genotypes.push(idx);
                }

                // Orders the occupied genotypes and removes duplicates
                self.occupied_genotypes.sort_unstable();
                self.occupied_genotypes.dedup();

                // Recodes the alleles such that the minor allele always has value true
                // First checks if an allele should be recoded
                let mut recode = vec![false; self.l];
                for allele in 0..self.l {
                    let mut count = 0;
                    for &g_idx in &self.occupied_genotypes {
                        let g = self.indices.to_sequence(g_idx);
                        count += if g[allele] { self.pop[g_idx] } else { 0 };
                    }
                    recode[allele] = count > self.size/2;
                }

                // Then recodes all relevant alleles sequentially
                for (allele, &recode_allele) in recode.iter().enumerate() {
                    if recode_allele {
                        let mut new_pop = vec![0usize; 2usize.pow(self.l as u32)];
                        let mut new_occupied_genotypes = Vec::<usize>::new();

                        for &g_idx in &self.occupied_genotypes {
                            let mut g = Vec::<bool>::from(self.indices.to_sequence(g_idx));
                            g[allele] = !g[allele];

                            let new_idx = self.indices.to_index(&g);
                            new_pop[new_idx] = self.pop[g_idx];
                            new_occupied_genotypes.push(new_idx)
                        }

                        self.pop = new_pop;
                        self.occupied_genotypes = new_occupied_genotypes;
                    }
                }
            },
            InitialPopulation::Monomorphic(g) => {
                // Transform the sequence to the corresponding index
                let idx = self.indices.to_index(&g);

                // Add the individuals to the population
                self.occupied_genotypes.push(idx);
                self.pop[idx] = self.size;
            },
            InitialPopulation::RandomMonomorphic => {
                // Pick a random genotype by choosing a random index
                let idx = rng.gen_range(0..self.pop.len());

                // Add the individuals to the population
                self.occupied_genotypes.push(idx);
                self.pop[idx] = self.size;
            }
        }
        // Remove repeated genotypes
        self.occupied_genotypes.sort_unstable();
        self.occupied_genotypes.dedup();

        // Check the alleles that may be fixed by chance in the initial population
        self.fixation = (0..self.l).map(|i| {
            let allele = self.is_fixed(i);
            if let Some(allele) = allele {
                Allele::Fixed { time: 0, allele }
            } else {
                Allele::NotFixed
            }
        }).collect();
    }

    /// Filters out extinct genotypes from the occupied_genotypes
    fn filter_extinct_genotypes(&mut self) {
        // Filter out the extinct genotypes
        self.occupied_genotypes = self
            .occupied_genotypes
            .iter()
            .filter(|&&i| self.pop[i] > 0)
            .map(|&i| i)
            .collect();
    }

    /// Saves current population to a file
    ///
    /// # Arguments:
    /// * `file` where to save the population
    /// * `i` landscape index
    /// * `recombination_map_id` recombination map identifier
    /// * `mutation_probability` mutation probability per locus
    /// * `c` replicate index
    /// * `t` current generation
    pub fn save(&self, file: &mut BufferedFile, i: usize, recombination_map_id: &str, mutation_probability: f64, c: usize, t: usize) -> std::io::Result<()> {
        if t == 0 {
            write!(file, "{}\t{}\t{}\t{}", i, recombination_map_id, mutation_probability, c)?;
        } else {
            write!(file, "{}\t{}\t{}\t{}\t{}", i, recombination_map_id, mutation_probability, c, t)?;
        }
        for &g in self.occupied_genotypes.iter() {
            write!(file, "\t")?;
            for &a in self.indices.to_sequence(g) {
                write!(file, "{}", if a { 1 } else { 0 })?;
            }
            write!(file, "\t{}", self.pop[g])?;
        }
        write!(file, "\n")?;
        Ok(())
    }

    /// Checks if there was any fixation and returns a list of fixed loci and 
    /// corresponding alleles
    ///
    /// # Arguments:
    /// * `t` current generation
    pub fn check_fixation(&mut self, t: usize) -> Vec<(usize, bool)> {
        let mut fixations = Vec::with_capacity(self.l);
        for i in 0..self.l {
            if let Some(fixation) = self.check_fixation_allele(t, i) {
                fixations.push(fixation)
            }
        }
        fixations
    }

    /// Checks if ith allele is fixed
    ///
    /// # Arguments:
    /// * `t` current generation
    #[inline]
    fn is_fixed(&self, i: usize) -> Option<bool> {
        let initial_allele = self.indices.get(self.occupied_genotypes[0], i);
        for j in 1..self.occupied_genotypes.len() {
            if self.indices.get(self.occupied_genotypes[j], i) != initial_allele {
                return None
            }
        }
        Some(initial_allele)
    }

    /// Checks if there was a fixation at the ith allele
    ///
    /// # Arguments:
    /// * `t` current generation
    /// * `i` index of allele to check
    pub fn check_fixation_allele(&mut self, t: usize, i: usize) -> Option<(usize, bool)> {
        // check if the allele is fixed
        let allele = self.is_fixed(i);

        // update the record accordingly
        match self.fixation[i] {
            Allele::NotFixed =>
                if let Some(allele) = allele {
                    self.fixation[i] = Allele::Fixed {
                        time: t,
                        allele,
                    };
                    return Some((i, allele));
                },
            Allele::Fixed { time: _, allele: _ } =>
                if allele == None {
                    self.fixation[i] = Allele::NotFixed;
                }
        }
        None
    }

    /// Gets the fixation status of the ith allele
    ///
    /// # Arguments:
    /// * `i` index of the allele
    pub fn get_fixation(&self, i: usize) -> Allele {
        self.fixation[i]
    }

    /// Count the number of fixed alleles in the population
    pub fn count_fixations(&self) -> usize {
        self.fixation
            .iter()
            .filter(|&x| match x {
                Allele::Fixed { time: _, allele: _ } => true,
                _ => false,
            })
            .count()
    }

    /// Counts the number of genotypes with frequency larger than zero in the population
    pub fn active_genotypes(&self) -> &[usize] {
        &self.occupied_genotypes
    }

    /// Gets the mean fitness of the population
    ///
    /// # Arguments:
    /// * `landscape` fitness landscape
    pub fn mean_fitness(&self, landscape: &FitnessLandscape) -> f64 {
        self.occupied_genotypes
            .iter()
            .map(|&i| -> f64 {
                let n = self.pop[i];
                (n as f64) * landscape.fitness_index(i)
            })
            .sum::<f64>()
            / (self.size as f64)
    }

    /// Gets the variance in fitness of the population
    ///
    /// # Arguments:
    /// * `landscape` fitness landscape
    /// * `f` current population fitness
    pub fn variance_fitness(&self, landscape: &FitnessLandscape, f: f64) -> f64 {
        let f2 = self.occupied_genotypes
            .iter()
            .map(|&i| -> f64 {
                let n = self.pop[i];
                (n as f64) * landscape.fitness_index(i) * landscape.fitness_index(i)
            })
            .sum::<f64>()
            / (self.size as f64);

        f2 - f * f
    }

    /// Gets the median fitness of the population and the corresponding genotype
    ///
    /// # Arguments:
    /// * `landscape` fitness landscape
    pub fn median_fitness(&self, landscape: &FitnessLandscape) -> (usize, f64) {
        let mut gf = self.occupied_genotypes
            .iter()
            .map(|&g| {
                (g, landscape.fitness_index(g))
            })
            .collect::<Vec<(usize, f64)>>();
        gf.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        let mut total = 0;
        let mut gf = gf.iter();
        loop {
            if let Some((g, f)) = gf.next() {
                total += self.pop[*g];
                if total >= self.size / 2 { break (*g, *f) }
            }
        }
    }

    /// Gets the shannon entropy (shannon diversity) of the population
    pub fn shannon_entropy(&self) -> f64 {
        let size = self.size as f64;
        self.occupied_genotypes
            .iter()
            .map(|&i| -> f64 {
                let n = self.pop[i];
                (n as f64) / size
            })
            .map(|f| -> f64 {
                if f > 0. {
                    - f * f.ln()
                } else {
                    0.
                }
            })
            .sum()
    }

    /// Gets the haplotype diversity of the population
    pub fn haplotype_diversity(&self) -> f64 {
        let size = self.size as f64;
        let f2: f64 = self.occupied_genotypes
            .iter()
            .map(|&i| {
                let n = self.pop[i];
                (n as f64) / size
            })
            .map(|f| { f * f })
            .sum();
        1. - f2
    }

    /// Gets the maximum fitness of the population
    ///
    /// # Arguments:
    /// * `landscape` fitness landscape
    pub fn maximum_fitness(&self, landscape: &FitnessLandscape) -> (usize, f64) {
        let (mut max, mut fmax) = (0usize, std::f64::NEG_INFINITY);
        for &g in self.occupied_genotypes.iter() {
            let fitness = landscape.fitness_index(g);
            if fitness > fmax {
                fmax = fitness;
                max = g;
            }
        }
        (max, fmax)
    }

    /// Gets the maximum fitness of the population
    ///
    /// # Arguments:
    /// * `landscape` fitness landscape
    pub fn minimum_fitness(&self, landscape: &FitnessLandscape) -> (usize, f64) {
        let (mut min, mut fmin) = (0usize, std::f64::INFINITY);
        for &g in self.occupied_genotypes.iter() {
            let fitness = landscape.fitness_index(g);
            if fitness < fmin {
                fmin = fitness;
                min = g;
            }
        }
        (min, fmin)
    }

    /// Converts an index to the corresponding sequence
    ///
    /// # Arguments:
    /// * `idx` index of a sequence
    pub fn sequence(&self, idx: usize) -> &[bool] {
        self.indices.to_sequence(idx)
    }

    /// Computes the current major genotype in the population
    pub fn major_genotype(&self) -> Vec<bool> {
        let mut counts = vec![0usize; self.l];
        for &g_idx in self.occupied_genotypes.iter() {
            let n = self.pop[g_idx];
            for (i, &allele) in self.sequence(g_idx).iter().enumerate() {
                if allele { counts[i] += n }
            }
        }
        let threshold = self.size/2;
        counts.iter().map(|&i| i > threshold).collect()
    }
}

impl std::ops::Index<usize> for Population {
    type Output = usize;

    fn index(&self, idx: usize) -> &Self::Output {
        &self.pop[idx]
    }
}
