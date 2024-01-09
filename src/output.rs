// author: André Amado

//! Output module

use crate::{
    fitness_model::FitnessModel,
    fitness_landscape::FitnessLandscape,
    population::{Population, Allele, InitialPopulation},
};

use clap::ArgMatches;

use std::{
    io::{Write, ErrorKind},
    fs,
    convert::TryInto
};
type BufferedFile = std::io::BufWriter<std::fs::File>;

use toml::Value;

///////////////////////////////////////////////////////////////////////////////////////////////////
// auxiliary functions and data structures to write to files and read command line parameters
// unimportant to the Biology
fn open_file(filename: String) -> BufferedFile {
    let path = std::path::Path::new(&filename);
    let file = match std::fs::File::create(&path) {
        Err(why) => panic!("couldn't create {}: {}", path.display(), why),
        Ok(file) => file,
    };
    std::io::BufWriter::new(file)
}

fn write_to_file(file: &mut BufferedFile, contents: &str) {
    match file.write(contents.as_bytes()) {
        Ok(_) => {}
        Err(why) => { panic!("couldn't write data to file: {}", why) }
    };
}

/// Enum that lists the different output types available
#[derive(PartialEq, Clone, Copy)]
pub enum OutputType {
    FitnessMean,
    FitnessVariance,
    FitnessMedian,
    FitnessMaximum,
    FitnessMinimum,
    FitnessMaximumWithGenotype,
    FitnessMinimumWithGenotype,
    FitnessMedianWithGenotype,
    ShannonEntropy,
    HaplotypeDiversity,
    NumberGenotypes,
    MajorGenotype,
    FixationsCount,
}

impl OutputType {
    /// Returns a string containing the corresponding statistics
    ///
    /// # Arguments:
    /// * `population`
    /// * `fitness_landscape`
    pub fn to_string(&self, population: &Population, fitness_landscape: &FitnessLandscape) -> String {
        match self {
            Self::FitnessMean           => format!("\t{}", population.mean_fitness(&fitness_landscape)),
            Self::FitnessVariance       => format!("\t{}", population.variance_fitness(&fitness_landscape, population.mean_fitness(&fitness_landscape))),
            Self::FitnessMedian         => {
                let (_, median_f) = population.median_fitness(&fitness_landscape);
                format!("\t{}", median_f)
            },
            Self::FitnessMaximum        => {
                let (_, maxf) = population.maximum_fitness(&fitness_landscape);
                format!("\t{}", maxf)
            },
            Self::FitnessMinimum        => {
                let (_, minf) = population.minimum_fitness(&fitness_landscape);
                format!("\t{}", minf)
            },
            Self::FitnessMaximumWithGenotype => {
                let (g, maxf) = population.maximum_fitness(&fitness_landscape);
                format!("\t{}\t{}", self.genotype_idx_to_string(&population, g), maxf)
            },
            Self::FitnessMinimumWithGenotype => {
                let (g, minf) = population.minimum_fitness(&fitness_landscape);
                format!("\t{}\t{}", self.genotype_idx_to_string(&population, g), minf)
            },
            Self::FitnessMedianWithGenotype  => {
                let (g, median_f) = population.median_fitness(&fitness_landscape);
                format!("\t{}\t{}", self.genotype_idx_to_string(&population, g), median_f)
            },
            Self::ShannonEntropy        => format!("\t{}", population.shannon_entropy()),
            Self::HaplotypeDiversity    => format!("\t{}", population.haplotype_diversity()),
            Self::NumberGenotypes       => format!("\t{}", population.active_genotypes().len()),
            Self::MajorGenotype         => {
                format!("\t{}", self.genotype_to_string(&population, &population.major_genotype()))
            },
            Self::FixationsCount        => format!("\t{}", population.count_fixations()),
        }
    }

    /// Converts a genotype index into a string containing the sequence
    ///
    /// # Arguments:
    /// * `population`
    /// * `idx` index of the genotype
    fn genotype_idx_to_string(&self, population: &Population, idx: usize) -> String {
        let mut g = String::with_capacity(population.get_l());
        for &allele in population.sequence(idx) {
            g.push(if allele { '1' } else { '0' })
        }
        g
    }

    /// Converts a genotype sequence into a string
    ///
    /// # Arguments:
    /// * `population`
    /// * `idx` index of the genotype
    fn genotype_to_string(&self, population: &Population, genotype: &[bool]) -> String {
        let mut g = String::with_capacity(population.get_l());
        for &allele in genotype.iter() {
            g.push(if allele { '1' } else { '0' })
        }
        g
    }

    /// Returns a string containing the label of the statistic
    pub fn label(&self) -> String {
        match self {
            Self::FitnessMean                => String::from("mean_fitness"),
            Self::FitnessVariance            => String::from("var_fitness"),
            Self::FitnessMedian              => String::from("median_fitness"),
            Self::FitnessMaximum             => String::from("max_fitness"),
            Self::FitnessMinimum             => String::from("min_fitness"),
            Self::FitnessMaximumWithGenotype => String::from("g_max\tmax_fitness"),
            Self::FitnessMinimumWithGenotype => String::from("g_min\tmin_fitness"),
            Self::FitnessMedianWithGenotype  => String::from("g_median\tmedian_fitness"),
            Self::ShannonEntropy             => String::from("shannon_entropy"),
            Self::HaplotypeDiversity         => String::from("haplotype_diversity"),
            Self::NumberGenotypes            => String::from("genotype_count"),
            Self::MajorGenotype              => String::from("major_genotype"),
            Self::FixationsCount             => String::from("fixations_count"),
        }
    }
}

/// Stores the output data information and and pointers to the output files
pub struct Output {
    data:   Vec<OutputType>,
    period: usize,
    population_period: usize,
    mutation_probability: f64,
    file_details:             Option<BufferedFile>,
    file_fixations:           Option<BufferedFile>,
    file_fixation_details:    Option<BufferedFile>,
    file_initial_population:  Option<BufferedFile>,
    file_population:          Option<BufferedFile>,
    file_final:               Option<BufferedFile>
}

impl Output {
    /// Checks if a directory exists and creates it in case it does not
    ///
    /// # Arguments:
    /// * `dir` directory to create
    pub fn create_dir(dir: &str) {
        match std::fs::create_dir(dir) {
            Ok(()) => {},
            Err(e) => {
                match e.kind() {
                    ErrorKind::AlreadyExists    => (), // No problem if the folder already exists!
                    ErrorKind::PermissionDenied => panic!("Error: permission to create directory \"{}\" denied.", dir),
                    _ => panic!("Error: failed to create directory \"{}\".", dir),
                }
            }
        }
    }

    fn generate_filename(name: &str, matches: &ArgMatches, model: &FitnessModel, identifier: &String) -> String {
        name.replace("%m",  &model.name())
            .replace("%M",  &model.short_description())
            .replace("%L",  &model.get_nloci().to_string())
            .replace("%N",  *matches.get_one("population_size").unwrap())
            .replace("%l",  *matches.get_one("nlandscapes").unwrap())
            .replace("%id", identifier)
            .replace("%I",  &InitialPopulation::from_args(&matches).short_description())
            .replace("%i",  &InitialPopulation::from_args(&matches).shorter_description())
    }

    /// Creates an Output struct from the command-line arguments
    ///
    /// # Arguments:
    /// * `matches` matches to the command-line arguments
    /// * `model` fitness model
    /// * `identifier` identifier string
    pub fn from_args(matches: &ArgMatches, model: &FitnessModel, identifier: &String) -> Output {
        let configuration_filename = matches.get_one::<String>("output_configuration").unwrap();
        let configuration_file = if configuration_filename == "default" {
            None
        } else {
            let file = match fs::read_to_string(configuration_filename) {
                Ok(f)  => f,
                Err(_) => panic!("Error: could not read file `{}`", configuration_filename)
            };
            Some(file.parse::<Value>().unwrap())
        };

        let mutation_probability = *matches.get_one("mutation_probability").unwrap_or(&0.);

        let mut period: usize = 1;
        Self::create_dir("data");
        let data = match &configuration_file {
            // If no configuration file is present return the default options
            None => vec![
                        OutputType::FitnessMean,        OutputType::FitnessVariance,
                        OutputType::HaplotypeDiversity, OutputType::FixationsCount,
                    ],
            Some(config) => {
                let mut data = vec![];

                if let Some(details) = config.get("details") {
                    // Checks if a period is defined and reads it
                    if let Some(p) = details.get("period") {
                        if let Some(p) = p.as_integer() {
                            period = p.try_into().unwrap();
                        } else {
                            // If it is defined and cannot be read as an interger give an error
                            println!("Warning: could not recognize a period for saving population details. Period set to 1.");
                        }
                    };

                    // Checks if the save flag is set
                    let save = if let Some(save) = details.get("save") {
                        save.as_bool().unwrap_or(false)
                    } else {
                        false
                    };

                    if save {
                        // Tries to read the options
                        if let Some(options) = details.get("options") {
                            for line in options.as_array().unwrap_or(&vec![]) {
                                let info = match &line.as_str().unwrap().to_lowercase()[..] {
                                    "mean_fitness"                   => OutputType::FitnessMean,
                                    "variance_fitness"               => OutputType::FitnessVariance,
                                    "median_fitness"                 => OutputType::FitnessMedian,
                                    "maximum_fitness"                => OutputType::FitnessMaximum,
                                    "minimum_fitness"                => OutputType::FitnessMinimum,
                                    "maximum_fitness_with_genotype"  => OutputType::FitnessMaximumWithGenotype,
                                    "minimum_fitness_with_genotype"  => OutputType::FitnessMinimumWithGenotype,
                                    "median_fitness_with_genotype"   => OutputType::FitnessMedianWithGenotype,
                                    "shannon_entropy"                => OutputType::ShannonEntropy,
                                    "haplotype_diversity"            => OutputType::HaplotypeDiversity,
                                    "genotypes_count"                => OutputType::NumberGenotypes,
                                    "major_genotype"                 => OutputType::MajorGenotype,
                                    "fixations_count"                => OutputType::FixationsCount,
                                    _ => panic!("Error: unrecognized output: {}", line)
                                };
                                data.push(info);
                            };
                        };
                    };
                }
                data
            }
        };

        let file_details = if data.len() > 0 {
            Self::create_dir("data/population_details");

            // Define the default file name
            let mut filename = if identifier.is_empty(){
                format!("data/population_details/{}.dat", model.short_description())
            } else {
                format!("data/population_details/{}_{}.dat", model.short_description(), identifier)
            };

            // if an alternative file name is given, replace the file name
            if let Some(config) = &configuration_file {
                if let Some(output) = config["details"].get("output_filename") {
                    if let Some(output_filename) = output.as_str() {
                        filename =
                            String::from("data/population_details/")
                            + &Self::generate_filename(output_filename, matches, model, identifier)
                    }
                }
            };

            // open file
            let mut file = open_file(filename);

            // save the labels
            let mut label = String::from("#landscape_id\trecombination_rate\tmutation_probability\treplicate\tt");
            for output in data.iter() {
                label += "\t";
                label += &output.label();
            }
            label += "\n";
            write_to_file(&mut file, &label);

            Some(file)
        } else {
            None
        };


        let mut save_fixations = false;
        let mut fixations_filename = if identifier.is_empty() {
            format!("data/fixations/{}.dat", model.short_description())
        } else {
            format!("data/fixations/{}_{}.dat", model.short_description(), identifier)
        };

        if let Some(conf) = &configuration_file {
            if let Some(fixations) = conf.get("fixations") {
                save_fixations = if let Some(save) = fixations.get("save") {
                    save.as_bool().unwrap_or(false)
                } else {
                    false
                };

                if let Some(name) = fixations.get("output_filename") {
                    fixations_filename =
                        String::from("data/fixations/")
                        + &Self::generate_filename(name.as_str().unwrap(), matches, model, identifier);
                }
            }
        };

        if save_fixations {
            println!("Warning: [fixations] option is deprecated and does not work correctly when paired with a nonzero mutation rate. Please use the [fixation_details] option instead.")
        }

        let file_fixations: Option<std::io::BufWriter<fs::File>> = if save_fixations {
            Self::create_dir("data/fixations");

            let mut file = open_file(fixations_filename);
            write_to_file(&mut file, "#landscape_id\trecombination_rate\tmutation_probability\treplicate\tfixation_time allele\n");

            Some(file)
        } else {
            None
        };


        let mut save_fixation_details = true;
        let mut fixation_details_filename = if identifier.is_empty() {
            format!("data/fixation_details/{}.dat", model.short_description())
        } else {
            format!("data/fixation_details/{}_{}.dat", model.short_description(), identifier)
        };

        if let Some(conf) = &configuration_file {
            if let Some(fixations) = conf.get("fixation_details") {
                save_fixation_details = if let Some(save) = fixations.get("save") {
                    save.as_bool().unwrap_or(false)
                } else {
                    false
                };

                if let Some(name) = fixations.get("output_filename") {
                    fixation_details_filename =
                        String::from("data/fixation_details/")
                        + &Self::generate_filename(name.as_str().unwrap(), matches, model, identifier);
                }
            }
        };

        let file_fixation_details: Option<std::io::BufWriter<fs::File>> = if save_fixation_details {
            Self::create_dir("data/fixation_details");

            let mut file = open_file(fixation_details_filename);
            write_to_file(&mut file, "#landscape_id\trecombination_rate\tmutation_probability\treplicate\tfixation_time\tlocus\tallele\n");

            Some(file)
        } else {
            None
        };


        let mut save_initial_population = true;
        let mut initial_population_filename = if identifier.is_empty(){
            format!("data/initial_population/{}.dat", model.short_description())
        } else {
            format!("data/initial_population/{}_{}.dat", model.short_description(), identifier)
        };

        if let Some(conf) = &configuration_file {
            if let Some(initial_population) = conf.get("initial_population") {
                save_initial_population = if let Some(save) = initial_population.get("save") {
                    save.as_bool().unwrap_or(false)
                } else {
                    false
                };

                if let Some(name) = initial_population.get("output_filename") {
                    initial_population_filename =
                        String::from("data/initial_population/")
                        + &Self::generate_filename(name.as_str().unwrap(), matches, model, identifier);
                }
            }
        };

        let file_initial_population = if save_initial_population {
            Self::create_dir("data/initial_population");

            let mut file = open_file(initial_population_filename);
            write_to_file(&mut file, "#landscape_id\trecombination_rate\tmutation_probability\treplicate\tgenotype\tcount\t...\n");
            
            Some(file)
        } else {
            None
        };

        let mut save_population = false;
        let mut population_period: usize = 0;
        let mut population_filename = if identifier.is_empty(){
            format!("data/periodic_population/{}.dat", model.short_description())
        } else {
            format!("data/periodic_population/{}_{}.dat", model.short_description(), identifier)
        };

        if let Some(conf) = &configuration_file {
            if let Some(population) = conf.get("population") {
                save_population = if let Some(period_temp) = population.get("period") {
                    population_period = period_temp.as_integer().unwrap_or(0).try_into().unwrap();

                    // If period is not 0 return true else return false
                    population_period != 0
                } else {
                    false
                };

                if let Some(name) = population.get("output_filename_periodic_population") {
                    population_filename = String::from("data/periodic_population/") + name.as_str().unwrap();
                }
            }
        };

        let file_population = if save_population {
            Self::create_dir("data/periodic_population");

            let mut file = open_file(population_filename);
            write_to_file(&mut file, "#landscape_id\trecombination_rate\tmutation_probability\treplicate\tgeneration\tgenotype\tcount\t...\n");
            Some(file)
        } else {
            None
        };

        let mut save_final_statistics: bool = true;
        let mut final_statistics_filename = if identifier.is_empty(){
            format!("data/data/{}.dat", model.short_description())
        } else {
            format!("data/data/{}_{}.dat", model.short_description(), identifier)
        };

        if let Some(conf) = &configuration_file {
            if let Some(final_statistics) = conf.get("final_statistics") {
                save_final_statistics = if let Some(save) = final_statistics.get("save") {
                    save.as_bool().unwrap_or(false)
                } else {
                    false
                };

                if let Some(name) = final_statistics.get("output_filename") {
                    final_statistics_filename =
                        String::from("data/data/")
                        + &Self::generate_filename(name.as_str().unwrap(), matches, model, identifier);
                }
            }
        };

        let file_final = if save_final_statistics {
            Self::create_dir("data/data");

            let mut file = open_file(final_statistics_filename);
            write_to_file(&mut file, "#landscape_id\trecombination_rate\tmutation_probability\treplicate\ttfinal\tfound_maximum\tadaptive_load\tfinal_fitness\n");
            Some(file)
        } else {
            None
        };

        Self {
            data,
            period,
            mutation_probability,
            file_details,
            file_fixations,
            file_fixation_details,
            file_initial_population,
            population_period,
            file_population,
            file_final
        }
    }

    /// Saves the population details specified by the user in the configuration file
    /// (option [details])
    ///
    /// # Arguments:
    /// * `i` landscape index
    /// * `recombination` recombination map identifier
    /// * `replicate` replicate index
    /// * `t` generation
    /// * `population`
    /// * `fitness_landscape`
    pub fn save_details(&mut self, i: usize, recombination_map_id: &str, replicate: usize, t: usize, population: &Population, fitness_landscape: &FitnessLandscape) {
        if let Some(file) = self.file_details.as_mut() {
            if t % self.period == 0 {
                let mut data = format!("{}\t{}\t{}\t{}\t{}", i, recombination_map_id, self.mutation_probability, replicate, t);
                for output in self.data.iter() {
                    data += output.to_string(population, fitness_landscape).as_str()
                }
                data += "\n";
                write_to_file(file, &data)
            }
        }
    }

    /// Saves the fixation data specified by the user in the configuration file
    /// (option [fixations])
    ///
    /// # Arguments:
    /// * `i` landscape index
    /// * `recombination` recombination map identifier
    /// * `replicate` replicate index
    /// * `population`
    pub fn save_fixations(&mut self, i: usize, recombination_map_id: &str, replicate: usize, population: &Population) {
        if let Some(file) = self.file_fixations.as_mut() {
            let mut fixations_line = format!("{}\t{}\t{}\t{}", i, recombination_map_id, self.mutation_probability, replicate);
            for j in 0..population.get_l() {
                let allele_info = match population.get_fixation(j) {
                    Allele::NotFixed => "-1 -1".to_string(),
                    Allele::Fixed { time, allele } => {
                        format!("{} {}", time, if allele { 1 } else { 0 })
                    }
                };
                fixations_line = format!("{}\t{}", fixations_line, allele_info);
            }
            fixations_line = format!("{}\n", fixations_line);
            write_to_file(file, &fixations_line);
        }
    }

    /// Saves the fixation details specified by the user in the configuration file
    /// (option [fixation_details])
    ///
    /// # Arguments:
    /// * `i` landscape index
    /// * `recombination` recombination map identifier
    /// * `replicate` replicate index
    /// * `t` current generation
    /// * `fixations` vector containing the fixations in the current generation
    pub fn save_fixation_details(&mut self, i: usize, recombination_map_id: &str, replicate: usize, t: usize, fixations: Vec<(usize, bool)>) {
        if let Some(file) = self.file_fixation_details.as_mut() {
            for (locus, allele) in fixations {
                let fixations_line = format!("{}\t{}\t{}\t{}\t{}\t{}\t{}\n", i, recombination_map_id, self.mutation_probability, replicate, t, locus, if allele { 1 } else { 0 });
                write_to_file(file, &fixations_line);
            }
        }
    }

    /// Saves the initial population as specified by the user in the configuration file
    /// (option [population] -> initial)
    ///
    /// # Arguments:
    /// * `i` landscape index
    /// * `recombination` recombination map identifier
    /// * `replicate` replicate index
    /// * `population`
    pub fn save_initial_population(&mut self, i: usize, recombination_map_id: &str, replicate: usize, population: &Population) {
        if let Some(file) = self.file_initial_population.as_mut() {
            population.save(file, i, recombination_map_id, self.mutation_probability, replicate, 0).unwrap();
        }
    }

    /// Saves the population as specified by the user in the configuration file
    /// (option [population] -> period)
    ///
    /// # Arguments:
    /// * `i` landscape index
    /// * `recombination` recombination map identifier
    /// * `replicate` replicate index
    /// * `t` generation
    /// * `population`
    pub fn save_population(&mut self, i: usize, recombination_map_id: &str, replicate: usize, t: usize, population: &Population) {
        if let Some(file) = self.file_population.as_mut() {
            if t % self.population_period == 0 {
                population.save(file, i, recombination_map_id, self.mutation_probability, replicate, t).unwrap();
            }
        }
    }

    /// Saves the final data as specified by the user in the configuration file
    /// (option [final_statistics])
    ///
    /// # Arguments:
    /// * `i` landscape index
    /// * `recombination` recombination map identifier
    /// * `replicate` replicate index
    /// * `final_genotype` final genotype index
    /// * `fitness_landscape`
    pub fn save_final_data(&mut self, i: usize, recombination_map_id: &str, replicate: usize, t: usize, final_genotype: usize, fitness_landscape: &FitnessLandscape) {
        if let Some(file) = self.file_final.as_mut() {
            let max_sequence = fitness_landscape.get_indices().to_sequence(final_genotype);
            let (max, maxf) = fitness_landscape.maximum();

            let found_maximum = if final_genotype == max { 1 } else { 0 };
            let final_fitness = fitness_landscape.fitness(max_sequence);

            let adaptive_load = (maxf - final_fitness) / maxf;

            let data_line =
                format!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                    i, recombination_map_id, self.mutation_probability, replicate, t, found_maximum, adaptive_load, final_fitness
                );
            write_to_file(file, &data_line);
        }
    }
}
