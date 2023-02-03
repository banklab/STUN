//! Fitness landscape module

#![allow(dead_code)]

use serde::{Deserialize, Serialize};

use std::{
    io::{Write, BufReader, BufWriter, BufRead},
    path::Path,
};
use std::{
    fs::File,
    process::exit
};

use rand::{thread_rng, prelude::IteratorRandom};
use rand_distr::{Distribution, Normal, Uniform};

use crate::lib::fitness_model::FitnessModel;

/// Stores a fitness landscape
#[derive(Serialize, Deserialize)]
pub struct FitnessLandscape {
    model: FitnessModel,
    landscape: Vec<f64>,
    indices: Indices,
    l: usize
}

use crate::lib::matrix::Matrix;
use crate::lib::indices::Indices;

impl FitnessLandscape {
    /// Gets a fitness landscape from the model described in the argument
    ///
    /// # Arguments:
    /// * `model` fitness model
    pub fn from_model(model: &FitnessModel) -> FitnessLandscape {
        let mut rng = thread_rng();
        match model {
            &FitnessModel::Additive {
                mu, sigma, l
            } => {
                let indices = Indices::new(l);
                let size = 2usize.pow(l as u32);

                let normal_additive  = Normal::new(mu, sigma).unwrap();
                let additive_component: Vec<f64> = normal_additive.sample_iter(&mut rng)
                                                                  .take(l).collect();

                let mut landscape = vec![0.; size];
                for g in 0..size {
                    let seq = indices.to_sequence(g);
                    for (i, &allele) in seq.iter().enumerate() {
                        if allele { landscape[g] += additive_component[i] }
                    }
                    landscape[g] = landscape[g].exp();
                }

                FitnessLandscape {
                    model: model.clone(),
                    landscape,
                    indices,
                    l
                }
            },
            &FitnessModel::HouseOfCards {
                mu, sigma, l
            } => {
                let indices = Indices::new(l);
                let size = 2usize.pow(l as u32);

                let normal_epistatic = Normal::new(mu, sigma).unwrap();
                let landscape = normal_epistatic.sample_iter(&mut rng).take(size)
                                                .map(|f| f.exp()).collect();

                FitnessLandscape {
                    model: model.clone(),
                    landscape,
                    indices,
                    l
                }
            },
            &FitnessModel::RoughMountFuji {
                mu_a, sigma_a, sigma_e, l
            } => {
                let indices = Indices::new(l);
                let size = 2usize.pow(l as u32);

                let normal_additive  = Normal::new(mu_a, sigma_a).unwrap();
                let normal_epistatic = Normal::new(0.,   sigma_e).unwrap();

                let additive_component: Vec<f64> = normal_additive.sample_iter(&mut rng).take(l).collect();
                let mut landscape: Vec<f64> = normal_epistatic.sample_iter(&mut rng).take(size).collect();
                for g in 0..size {
                    let seq = indices.to_sequence(g);
                    for (i, &allele) in seq.iter().enumerate() {
                        if allele { landscape[g] += additive_component[i] }
                    }
                    landscape[g] = landscape[g].exp();
                }

                FitnessLandscape {
                    model: model.clone(),
                    landscape,
                    indices,
                    l
                }
            },
            &FitnessModel::NK {
                n, k
            } => {
                let indices = Indices::new(n);
                let size = 2usize.pow(n as u32);

                let uniform = Uniform::new(0_f64, 1_f64/n as f64);

                let mut fitness_contributions = Matrix::<f64>::new(n, 2usize.pow((k+1) as u32));
                for i in 0..n {
                    for j in 0..2usize.pow((k+1) as u32) {
                        fitness_contributions.set(i, j, uniform.sample(&mut rng))
                    }
                }

                let mut landscape = vec![0_f64; size];
                for i in 0..size {
                    let seq = indices.to_sequence(i);

                    for j in 0..n {
                        let mut pos = 0;
                        for ki in 0..k+1 {
                            pos += 2usize.pow(ki as u32) * usize::from(seq[(j + ki) % n])
                        }
                        landscape[i] += fitness_contributions.get(j, pos)
                    }
                }
                FitnessLandscape {
                    model: model.clone(),
                    landscape,
                    indices,
                    l: n
                }
            },
            &FitnessModel::Block {
                n, b, mu, sigma
            } => {
                let indices = Indices::new(n);
                let size = 2usize.pow(n as u32);

                let mut block_boundary = (1..n).choose_multiple(&mut rng, b-1);
                block_boundary.insert(0, 0);
                block_boundary.push(n);

                let normal = Normal::new(mu, sigma).unwrap();

                let mut block_fitness = Vec::<Vec<f64>>::with_capacity(b);
                let mut block_indices = Vec::<Indices>::with_capacity(b);
                for i in 0..b {
                    let block_size = block_boundary[i+1] - block_boundary[i];

                    let fitness: Vec<f64> = normal.sample_iter(&mut rng).take(2usize.pow(block_size as u32)).collect();
                    block_fitness.push(fitness);
                    block_indices.push(Indices::new(block_size));
                }

                let mut landscape = vec![0_f64; size];
                for i in 0..size {
                    let seq = indices.to_sequence(i);

                    for block in 0..b {
                        let intra_block_index = block_indices[block].to_index(&seq[block_boundary[block]..block_boundary[block+1]]);
                        landscape[i] += block_fitness[block][intra_block_index];
                    }
                    landscape[i] = landscape[i].exp()
                }

                FitnessLandscape {
                    model: model.clone(),
                    landscape,
                    indices,
                    l: n
                }
            },
            FitnessModel::Custom {
                ..
            } => {
                println!("Error: custom landscapes cannot be generated!");
                exit(0)
            }
        }
    }

    /// Gets the fitness of a given genotype (access by genotype sequence)
    ///
    /// # Arguments:
    /// * `seq` sequence of the genotype
    pub fn fitness(&self, seq: &[bool]) -> f64 {
        self.landscape[self.indices.to_index(seq)]
    }

    /// Gets the fitness of a given genotype (access by genotype index)
    ///
    /// # Arguments:
    /// * `idx` index of the genotype
    pub fn fitness_index(&self, idx: usize) -> f64 {
        self.landscape[idx]
    }

    /// Gets the Index struct that stores the index information
    pub fn get_indices(&self) -> &Indices {
        &self.indices
    }

    /// Gets a tuple containing the index and fitness of the genotype with the highest fitness
    pub fn maximum(&self) -> (usize, f64) {
        let mut max = 0;
        let mut maxf = f64::NEG_INFINITY;
        for (i, &f) in self.landscape.iter().enumerate() {
            if f > maxf {
                max = i;
                maxf = f;
            }
        }
        (max, maxf)
    }

    /// Gets a tuple containing the index and fitness of the genotype with the lowest fitness
    pub fn minimum(&self) -> (usize, f64) {
        let mut min = 0;
        let mut minf = f64::INFINITY;
        for (i, &f) in self.landscape.iter().enumerate() {
            if f < minf {
                min = i;
                minf = f;
            }
        }
        (min, minf)
    }

    /// Saves a binary representation of the fitness landscape that can be loaded later (this is an
    /// exact and more compact representation of the fitness landscape than the one provided by
    /// save() method)
    ///
    /// # Arguments:
    /// * `filename` name of the file in which to save the fitness landscape
    pub fn save_raw(&self, filename: &str) {
        let file = File::create(filename).unwrap();
        bincode::serialize_into(&file, &self).unwrap();
    }

    /// Loads the binary representation of a fitness landscape
    ///
    /// # Arguments:
    /// * `filename` name of the file containing the fitness landscape
    pub fn load_raw(filename: &str) -> Self {
        let file = match File::open(filename) {
            Ok(f)  => f,
            Err(_) => panic!("File {} does not exist.", filename)
        };
        bincode::deserialize_from(file).unwrap()
    }

    /// Saves a human readable representation of the fitness landscape. It may not be an exact
    /// representation of the fitness landscape due to round error in conversion from binary to
    /// base 10
    ///
    /// # Arguments:
    /// * `filename` name of the file in which to save the fitness landscape
    pub fn save(&self, filename: &str) {
        let path = Path::new(filename);
        let display = path.display();
        let file = match File::create(&path) {
            Err(why) => panic!("couldn't create {}: {}", display, why),
            Ok(file) => file,
        };
        let mut writer = BufWriter::new(file);

        match writer.write(format!("# {}\n", &self.model.long_description().replace("\n", "\n# ")).as_bytes()) {
            Ok(_) => {}
            Err(why) => {
                panic!("couldn't write data to file: {}", why)
            }
        };

        for seq in self.indices.seq_iter() {
            let mut line = String::new();
            for &allele in seq.iter() {
                line = format!("{} {}", line, if allele { 1 } else { 0 });
            }
            line = format!("{} {}\n", line, self.fitness(seq));
            line.remove(0);
            match writer.write(line.as_bytes()) {
                Ok(_) => {}
                Err(why) => {
                    panic!("couldn't write data to file: {}", why)
                }
            };
        }
    }

    /// Loads a human readable version representation of the fitness landscape
    ///
    /// # Arguments:
    /// * `filename` name of the file in which to save the fitness landscape
    /// * `model`
    pub fn load(filename: &str, model: &FitnessModel) -> Self {
        let l = model.get_nloci();
        let file = match File::open(&filename) {
            Ok(f)  => f,
            Err(_) => panic!("Fitness landscape file {} not found.", &filename)
        };

        let indices = Indices::new(l);
        let size = 2usize.pow(l as u32);

        let mut landscape = vec![0.; size];

        let mut current_line = 0;
        for line in BufReader::new(file).lines() {
            current_line += 1;

            // read and trim line
            let line_t = match line {
                Ok(line) => line.trim().to_string(),
                Err(_)   => panic!("Error reading line {} of file {}.", current_line, filename)
            };

            // ignore comments and empty lines
            if line_t.starts_with('#') || line_t.len() == 0 { continue }

            // split the line in individual numbers, ignoring whitespace
            let line_t: Vec<&str> = line_t.split_whitespace().collect();

            // ignore Magellan introductory lines
            if line_t.len() == l { continue }

            if line_t.len() != l+1 {
                panic!("Error reading line {} of file {}. Expected {} numbers, found {}.", current_line, &filename, l+1, line_t.len())
            }
            let mut line_t = line_t.iter();

            // reads each allele in the genotype
            let mut g = Vec::<bool>::with_capacity(l);
            for _ in 0..l {
                let gi = line_t.next().unwrap();
                g.push(match gi {
                    &"0" => false,
                    &"1" => true,
                    _    => panic!("Error reading line {} of file {}. Expected 0 or 1 found {}", current_line, &filename, gi)
                });
            }

            // reads the fitness
            let f = line_t.next().unwrap();
            let f = match f.parse::<f64>() {
                Ok(f)  => f,
                Err(_) => panic!("Error reading line {} of file {}. Could not convert {} to a floating point number.", current_line, &filename, f)
            };

            // store the fitness to the correct genotype
            landscape[indices.to_index(&g)] = f;
        }
        FitnessLandscape {
            model: model.clone(),
            landscape,
            indices,
            l
        }
    }
}
