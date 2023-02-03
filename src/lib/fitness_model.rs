//! Fitness model module

#![allow(dead_code)]

use serde::{Deserialize, Serialize};
use clap::{Command, Arg, ArgMatches, ArgGroup, value_parser, ValueHint};
use std::{
    io::{BufReader, BufRead},
    fs::File,
    process::exit
};
use glob::glob;

/// Stores the information about a fitness model and its parameters
#[derive(Serialize, Deserialize, Clone)]
pub enum FitnessModel {
    Additive {
        mu: f64,
        sigma: f64,
        l: usize
    },
    HouseOfCards {
        mu: f64,
        sigma: f64,
        l: usize
    },
    RoughMountFuji {
        mu_a: f64,
        sigma_a: f64,
        sigma_e: f64,
        l: usize
    },
    NK {
        n: usize,
        k: usize,
    },
    Block {
        n: usize,
        b: usize,
        mu: f64,
        sigma: f64
    },
    Custom {
        filenames: Vec<String>,
        l: usize
    }
}

impl FitnessModel {
    /// Gets a long description of the fitness model
    pub fn long_description(&self) -> String {
        match &self {
            Self::Additive {
                mu, sigma, l
            } => format!("Additive fitness landscape\n{} loci\nmean: {}\nstd: {}", l, mu, sigma),
            Self::HouseOfCards {
                mu, sigma, l
            } => format!("House of Cards fitness landscape\n{} loci\nmean: {}\nstd: {}", l, mu, sigma),
            Self::RoughMountFuji {
                mu_a, sigma_a, sigma_e, l
            } => format!("Rough Mount Fuji fitness landscape\n{} loci\nmean additive component: {}\nstd of additive component: {}\nstd of epistatic component: {}", l, mu_a, sigma_a, sigma_e),
            Self::NK {
                n, k
            } => format!("NK fitness landscape\nN: {}\nK: {}", n, k),
            Self::Block {
                n, b, mu, sigma
            } => format!("Block fitness landscape\nN: {}\nB: {}\nmean fitness contribution: {}\nstd fitness contribution: {}", n, b, mu, sigma),
            Self::Custom {
                ..
            } => format!("Custom fitness landscape"),
        }.to_string()
    }

    /// Gets a short description of the fitness model, containing no whitespace. Useful for
    /// inclusion, e.g., in a file name.
    pub fn short_description(&self) -> String {
        match &self {
            Self::Additive {
                mu, sigma, l
            } => format!("additive_L{}_mu{}_std{}", l, mu, sigma),
            Self::HouseOfCards {
                mu, sigma, l
            } => format!("HoC_L{}_mu{}_std{}", l, mu, sigma),
            Self::RoughMountFuji {
                mu_a, sigma_a, sigma_e, l
            } => format!("RMF_L{}_mu{}_stda{}_stdb{}", l, mu_a, sigma_a, sigma_e),
            Self::NK {
                n, k
            } => format!("NK_N{}_K{}", n, k),
            Self::Block {
                n, b, mu, sigma
            } => format!("block_N{}_B{}_mu{}_std{}", n, b, mu, sigma),
            Self::Custom {
                ..
            } => {
                format!("custom")
            },
        }.to_string()
    }

    /// Gets a short description of the fitness model, containing no whitespace. Useful for
    /// inclusion, e.g., in a file name.
    pub fn name(&self) -> String {
        match &self {
            Self::Additive { .. }       => "additive",
            Self::HouseOfCards { .. }   => "HoC",
            Self::RoughMountFuji { .. } => "RMF",
            Self::NK { .. }             => "NK",
            Self::Block { .. }          => "block",
            Self::Custom { .. }         => "custom",
        }.to_string()
    }

    /// Gets the clap subcommands that allow the user to input the fitness model information
    pub fn get_subcommands() -> [Command<'static>; 6] {[
        Command::new("additive")
            .about("Selects the additive model")
            .allow_negative_numbers(true)
            .subcommand_precedence_over_arg(true)
            .args([
                Arg::new("L")
                    .short('L')
                    .help("Number of loci")
                    .default_value("5")
                    .value_name("number_of_loci")
                    .value_parser(value_parser!(usize)),
                Arg::new("mu")
                    .short('m')
                    .long("mu")
                    .help("Mean")
                    .default_value("0.1")
                    .value_name("mu")
                    .value_parser(value_parser!(f64)),
                Arg::new("sigma")
                    .short('s')
                    .long("sigma")
                    .help("Standard deviation")
                    .default_value("0.1")
                    .value_name("sigma")
                    .value_parser(value_parser!(f64)),
            ])
            .mut_arg("help", |arg| arg.hide(true)),
            Command::new("HoC")
                .about("Selects the House of Cards model")
                .allow_negative_numbers(true)
                .subcommand_precedence_over_arg(true)
                .args([
                    Arg::new("L")
                        .short('L')
                        .help("Number of loci")
                        .default_value("5")
                        .value_name("number_of_loci")
                        .value_parser(value_parser!(usize)),
                    Arg::new("mu")
                        .short('m')
                        .long("mu")
                        .help("Mean")
                        .default_value("0.1")
                        .value_name("mu")
                        .value_parser(value_parser!(f64)),
                    Arg::new("sigma")
                        .short('s')
                        .long("sigma")
                        .help("Standard deviation")
                        .default_value("0.1")
                        .value_name("sigma")
                        .value_parser(value_parser!(f64)),
                ])
                .mut_arg("help", |arg| arg.hide(true)),
        Command::new("RMF")
            .about("Selects the Rough Mount Fuji model")
            .allow_negative_numbers(true)
            .subcommand_precedence_over_arg(true)
            .args([
                Arg::new("L")
                    .short('L')
                    .help("Number of loci")
                    .default_value("5")
                    .value_name("number_of_loci")
                    .value_parser(value_parser!(usize)),
                Arg::new("mu_a")
                    .short('m')
                    .long("mu_a")
                    .help("Mean additive component")
                    .default_value("0.1")
                    .value_name("mu_a")
                    .value_parser(value_parser!(f64)),
                Arg::new("sigma_a")
                    .short('s')
                    .long("sigma_a")
                    .help("Standard deviation of additive component")
                    .default_value("0.1")
                    .value_name("sigma_a")
                    .value_parser(value_parser!(f64)),
                Arg::new("sigma_e")
                    .short('S')
                    .long("sigma_e")
                    .help("Standard deviation of epistatic component")
                    .default_value("0.1")
                    .value_name("sigma_e")
                    .value_parser(value_parser!(f64)),
            ])
            .mut_arg("help", |arg| arg.hide(true)),
        Command::new("NK")
            .about("Selects the NK model")
            .allow_negative_numbers(true)
            .subcommand_precedence_over_arg(true)
            .args([
                Arg::new("N")
                    .short('N')
                    .help("Number of loci")
                    .default_value("5")
                    .value_name("number_of_loci")
                    .value_parser(value_parser!(usize)),
                Arg::new("K")
                    .short('K')
                    .help("Number of interacting loci [0, K-1]")
                    .default_value("0")
                    .value_name("K")
                    .value_parser(value_parser!(usize)),
            ])
            .mut_arg("help", |arg| arg.hide(true)),
            Command::new("block")
                .about("Selects the block model")
                .allow_negative_numbers(true)
                .subcommand_precedence_over_arg(true)
                .args([
                    Arg::new("N")
                        .short('N')
                        .help("Number of loci")
                        .default_value("5")
                        .value_name("number_of_loci")
                        .value_parser(value_parser!(usize)),
                    Arg::new("B")
                        .short('B')
                        .help("Number of blocks")
                        .default_value("1")
                        .value_name("number_of_blocks")
                        .value_parser(value_parser!(usize)),
                    Arg::new("mu")
                        .short('m')
                        .long("mu")
                        .help("Mean value of block fitness components")
                        .default_value("0.")
                        .value_name("mu")
                        .value_parser(value_parser!(f64)),
                    Arg::new("sigma")
                        .short('s')
                        .long("sigma")
                        .help("Standard deviation of block fitness components")
                        .default_value("0.1")
                        .value_name("sigma")
                        .value_parser(value_parser!(f64)),
                ])
                .mut_arg("help", |arg| arg.hide(true)),
        Command::new("custom")
            .about("Selects a custom fitness landscape provided by the user")
            .subcommand_precedence_over_arg(true)
            .args([
                Arg::new("filename")
                    .short('f')
                    .long("file")
                    .help("Name of a file containing the fitness landscape")
                    .value_name("file")
                    .value_hint(ValueHint::FilePath),
                Arg::new("list")
                    .long("list")
                    .help("Name of a file providing a list of fitness landscapes")
                    .value_name("file_with_list"),
            ])
            .group(
                ArgGroup::new("loading")
                         .args(&["filename", "list"])
                         .required(true)
            )
            .mut_arg("help", |arg| arg.hide(true)),
    ]}

    fn all_eq<T: PartialEq>(iter: impl IntoIterator<Item = T>) -> Option<T> {
        let mut iter = iter.into_iter();
        let first = iter.next()?;
        iter.all(|elem| elem == first).then(|| first)
    }

    /// Reads the model information from a clap argument match
    pub fn from_args(matches: &ArgMatches) -> Self {
        if let Some(custom) = matches.subcommand_matches("custom") {
            let filenames = if custom.is_present("filename") {
                let filename = custom.get_one::<String>("filename").unwrap().clone();
                if filename.contains("*") {
                    let file_matches = match glob(&filename) {
                        Ok(m)  => m,
                        Err(_) => {
                            println!("Error: pattern {} could not match any files.", filename);
                            exit(0)
                        }
                    };
                    file_matches.map(|f| f.unwrap().as_os_str().to_os_string().into_string().unwrap())
                                .filter(|f| !f.contains(".bin"))
                                .collect()
                } else {
                    vec![filename]
                }
            } else {
                let file_list = custom.get_one::<String>("list").unwrap();
                // open the file
                let file = match File::open(&file_list) {
                    Ok(f)  => f,
                    Err(_) => {
                        println!("Error: file {} not found.", &file_list);
                        exit(0)
                    }
                };
                // Reads the files from the list, discarding comments and empty lines
                BufReader::new(file).lines().filter_map(|line|
                    if line.as_ref().unwrap().trim().len() == 0 || line.as_ref().unwrap().trim().starts_with('#') {
                        None
                    } else {
                        Some(String::from(line.unwrap().trim()))
                    }
                ).collect()
            };

            let mut l_list = Vec::with_capacity(filenames.len());
            for filename in filenames.iter() {
                // open the file
                let file = match File::open(filename) {
                    Ok(f)  => f,
                    Err(_) => {
                        println!("Error: fitness landscape file {} not found.", filename);
                        exit(0)
                    }
                };

                let mut lines = BufReader::new(file).lines();

                let l = loop {
                    // get a new line
                    let line = match lines.next() {
                        Some(line) => line.unwrap(),
                        None       => panic!("Error reading file {}. Could not determine the number of loci.", &filename)
                    };

                    // ignore if line is empty
                    if line.trim().len() == 0 { continue }

                    // try to extract the first
                    let mut line = line.trim().split_whitespace();
                    let value = match line.next() {
                        Some(value) => value,
                        None        => panic!("Error reading file {}. Could not determine the number of loci.", &filename)
                    };

                    // ignore if the line is a comment
                    if value.starts_with('#') { continue }

                    let first_number = match value.parse::<usize>() {
                        Ok(number) => number,
                        Err(_)     => panic!("Error reading file {}. Could not parse number as integer.", &filename)
                    };

                    break match first_number {
                        0 | 1 => line.count(),     // line with a fitness value
                        2     => line.count() + 1, // Magellan format introductory line, also a valid way to determine loci number
                        _     => panic!("Error reading file {}. Possibly trying to load a file with more than two alleles per locus.", &filename)
                    };
                };
                l_list.push(l);
            }

            println!("Using the following custom fitness landscapes:");
            for f in &filenames { println!("{}", f) }
            println!("");

            if let Some(l) = Self::all_eq(l_list) {
                FitnessModel::Custom {
                    filenames: filenames.clone(), l
                }
            } else {
                println!("Error. Fitness landscapes have inconsistent loci numbers.");
                exit(0)
            }
        } else if let Some(nk) = matches.subcommand_matches("NK") {
            let n = *nk.get_one("N").unwrap();
            let k = *nk.get_one("K").unwrap();

            FitnessModel::NK {
                n, k
            }
        } else if let Some(block) = matches.subcommand_matches("block") {
            let n = *block.get_one("N").unwrap();
            let b = *block.get_one("B").unwrap();
            let mu = *block.get_one("mu").unwrap();
            let sigma = *block.get_one("sigma").unwrap();

            FitnessModel::Block {
                n, b, mu, sigma
            }
        } else if let Some(additive) = matches.subcommand_matches("additive") {
            let mu = *additive.get_one("mu").unwrap();
            let sigma = *additive.get_one("sigma").unwrap();
            let l = *additive.get_one("L").unwrap();
            FitnessModel::Additive {
                mu, sigma, l
            }
        } else if let Some(hoc) = matches.subcommand_matches("HoC") {
            let mu = *hoc.get_one("mu").unwrap();
            let sigma = *hoc.get_one("sigma").unwrap();
            let l = *hoc.get_one("L").unwrap();
            FitnessModel::HouseOfCards {
                mu, sigma, l
            }
        } else if let Some(rmf) = matches.subcommand_matches("RMF") {
            let mu_a = *rmf.get_one("mu_a").unwrap();
            let sigma_a = *rmf.get_one("sigma_a").unwrap();
            let sigma_e = *rmf.get_one("sigma_e").unwrap();
            let l = *rmf.get_one("L").unwrap();
            FitnessModel::RoughMountFuji {
                mu_a, sigma_a, sigma_e, l
            }
        } else {
            panic!("Model not recognized")
        }
    }

    /// Gets the number of loci in a fitness model
    pub fn get_nloci(&self) -> usize {
        match self {
            &Self::Additive { l, .. } => l,
            &Self::HouseOfCards { l, .. } => l,
            &Self::RoughMountFuji { l, .. } => l,
            &Self::NK { n, .. } => n,
            &Self::Block { n, ..} => n,
            &Self::Custom { l, .. } => l,
        }
    }

    /// Gets the list of file names where to find the fitness landscapes
    pub fn get_fitness_landscapes_list(&self, sample_landscapes: usize, identifier: &str) -> Vec<String> {
        match self {
            Self::Custom { filenames, .. } => filenames.clone(),
            _ => {
                (0..sample_landscapes).map(|i| format!(
                    "landscapes/{}_{}{}.bin.fl",
                    self.short_description(), i, identifier
                )).collect()
            }
        }
    }
}
