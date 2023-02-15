// author: André Amado

//! # STUN
//! Forward-time **S**imulation on **TU**nable fit**N**ess landscapes in recombining populations
//!
//! For details on the usage of this program check the associated
//! [manuscript](...) and the
//! [manual](https://bit.ly/stun_manual_binaries)

mod lib;

use glob::glob;
use clap::{Arg, ArgGroup, Command, ArgMatches, value_parser};
use std::{
    process::exit,
    path::Path
};

use crate::lib::{
    fitness_model::FitnessModel,
    population::InitialPopulation,
    fitness_landscape::FitnessLandscape,
    population::{Population, Recombination},
    output::Output,
};


fn main() -> Result<(), String> {
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // configure and read command line parameters
    let matches = get_command_line_matches();
    let model = FitnessModel::from_args(&matches);

    let nreplicates: usize = *matches.get_one("nreplicates").unwrap();
    let mut sample_landscapes: usize = *matches.get_one("nlandscapes").unwrap();

    let empty_string = String::new();
    let identifier = matches.get_one::<String>("identifier").unwrap_or(&empty_string);

    // Generates the landscapes if the model is not custom and they do not already exist
    match &model {
        crate::FitnessModel::Custom { filenames, .. } => {
            sample_landscapes = filenames.len();
            // Issue a warning recommending the usage of a run id for custom models
            if identifier.len() == 0 {
                println!("Warning: it is strongly recommended to define a run identifier when using custom models. It can be defined through the option --id.")
            }
        },
        _ => generate_landscapes(sample_landscapes, &identifier, &model)
    };
    if *matches.get_one::<bool>("generate_only").unwrap() {
        exit(0)
    }

    let fitness_landscape_list = model.get_fitness_landscapes_list(sample_landscapes, &identifier);

    let initial_population = InitialPopulation::from_args(&matches);
    let population_size: usize = *matches.get_one("population_size").unwrap();

    let recombination_map_list = Recombination::from_args(&matches, model.get_nloci())?;
    ///////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // open files to save data
    let mut output = Output::from_args(&matches, &model, &identifier);
    ///////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // Main loop
    // repeat for each landscape
    for (i, fitness_landscape_file) in fitness_landscape_list.iter().enumerate() {
        let fitness_landscape = if fitness_landscape_file.contains(".bin") {
            FitnessLandscape::load_raw(fitness_landscape_file)
        } else {
            FitnessLandscape::load(fitness_landscape_file, &model)
        };

        // run the populations for all the recombination rates listed by the user
        for recombination_map in recombination_map_list.iter() {
            let recombination_map_id = recombination_map.summary();

            // repeat for nreplicates runs
            let mut population = Population::new(model.get_nloci(), population_size, recombination_map.clone());
            for replicate in 0..nreplicates {
                population.generate_initial_population(&initial_population);
                output.save_initial_population(i, &recombination_map_id, replicate, &population);

                output.save_details(i, &recombination_map_id, replicate, 0, &population, &fitness_landscape);

                let mut t = 0;
                let final_genotype = loop {
                    ///////////////////////////////////////////////////////////////////////////////
                    // update the population to the new generation
                    // Biology is here!
                    population.recombination();
                    population.wright_fisher(&fitness_landscape);
                    ///////////////////////////////////////////////////////////////////////////////

                    // check for genotype fixations in the population
                    population.check_fixation(t);
                    t += 1;

                    ///////////////////////////////////////////////////////////////////////////////
                    // save population
                    output.save_population(i, &recombination_map_id, replicate, t, &population);
                    ///////////////////////////////////////////////////////////////////////////////

                    // stop if only one genotype remains
                    let active_genotypes = population.active_genotypes();
                    if active_genotypes.len() == 1 {
                        break active_genotypes[0];
                    }

                    ///////////////////////////////////////////////////////////////////////////////
                    // write results to a file
                    output.save_details(i, &recombination_map_id, replicate, t, &population, &fitness_landscape);
                    ///////////////////////////////////////////////////////////////////////////////
                };
                ///////////////////////////////////////////////////////////////////////////////////
                // write some more results to data files
                output.save_final_data(i, &recombination_map_id, replicate, t, final_genotype, &fitness_landscape);
                output.save_fixations(i, &recombination_map_id, replicate, &population);
                ///////////////////////////////////////////////////////////////////////////////////
            } //End of replicate
        } //End of recombination maps
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////

    Ok(())
}


fn get_command_line_matches() -> ArgMatches {
    let matches = Command::new("Haploid populations with standing variation")
        .version("0.2.0")
        .author("André Amado <andreamado1@gmail.com>")
        .about("Evolution in haploid populations with standing variation.")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand_help_heading("MODELS")
        .subcommand_value_name("MODEL")
        .args([
            Arg::with_name("nreplicates")
                .short('c')
                .long("replicates")
                .value_name("number_of_replicates")
                .help("Number of replicates per fitness landscape and recombination map")
                .default_value("10")
                .value_parser(value_parser!(usize))
                .display_order(1),
            Arg::new("nlandscapes")
                .long("landscapes")
                .short('l')
                .help("Number of fitness landscapes to use (this option is inactivated for custom landscapes)")
                .default_value("1")
                .value_name("value")
                .value_parser(value_parser!(usize))
                .display_order(1),
            Arg::new("identifier")
                .long("id")
                .help("Custom identifier (appended to the end of the filenames)")
                .value_name("id")
                .required(false)
                .display_order(3),
            Arg::with_name("population_size")
                .short('N')
                .long("size")
                .value_name("population_size")
                .help("Population size")
                .default_value("1000")
                .value_parser(value_parser!(usize))
                .display_order(4),
            Arg::with_name("neutralsfs")
                .long("neutralsfs")
                .short('n')
                .min_values(0)
                .max_values(2)
                .value_parser(value_parser!(f64))
                .help("Distributes the initial population according to a neutral site frequency spectrum distribution (can be followed by an optional drift threshold value and a second number for the probability that the minor allele is encoded by 1)")
                .conflicts_with_all(&["uniform", "generate_only"])
                .display_order(5),
            Arg::with_name("uniform")
                .long("uniform")
                .short('u')
                .help("Distributes the initial population uniformly with minor allele probability p")
                .value_name("p")
                .default_value("0.5")
                .value_parser(value_parser!(f64))
                .conflicts_with_all(&["neutralsfs", "generate_only"])
                .display_order(5),
            Arg::with_name("generate_only")
                .long("generate_only")
                .help("Only generates landscapes (does not simulate the adaptive process)")
                .conflicts_with_all(&["uniform", "neutralsfs"])
                .action(clap::ArgAction::SetTrue)
                .display_order(5),
            Arg::with_name("recombination_rates")
                .short('r')
                .long("recombination_rates")
                .value_name("rate")
                .help("List of recombination rates (separated by comma)")
                .takes_value(true)
                .multiple_values(true)
                .use_value_delimiter(true)
                .value_parser(value_parser!(f64))
                .display_order(6),
            Arg::with_name("recombination_map")
                .long("recombination_map")
                .takes_value(true)
                .help("List of recombination probabilities for each allele (separated by comma) or file name for file with lists of recombination probabilities (each recombination map in a line)")
                .display_order(6),
            Arg::with_name("output_configuration")
                .short('o')
                .long("output_conf")
                .help("Configuration file specifying the output details. If not present a default output will be assumed.")
                .default_value("default")
                .value_parser(value_parser!(String))
                .display_order(9),
        ])
        .group(ArgGroup::with_name("initial_population")
            .args(&["neutralsfs", "uniform", "generate_only"])
            .required(true)
            .multiple(false)
        )
        .group(ArgGroup::with_name("recombination")
            .args(&["recombination_map", "recombination_rates"])
            .required(false)
            .multiple(false)
        )
        .subcommands(FitnessModel::get_subcommands())
        .mut_subcommand("help", |subcmd| subcmd.about("Print this message or the help for a given model"))
        .subcommand_precedence_over_arg(true)
        .try_get_matches();

    match matches {
        Err(error)  => {
            match error.kind() {
                clap::ErrorKind::MissingSubcommand => println!("error: a fitness model is required but none was provided"),
                clap::ErrorKind::UnrecognizedSubcommand => println!("error: fitness model not recognized"),
                _ => error.exit()
            }
            exit(1)
        },
        Ok(matches) => matches
    }
}

fn generate_landscapes(sample_landscapes: usize, identifier: &str, model: &FitnessModel) {
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // count preexisting landscapes
    Output::create_dir("landscapes");
    let generic_filename = format!(
        "landscapes/{}_*{}.bin.fl", model.short_description(), identifier
    );
    let existing_landscapes: usize = glob(&generic_filename).unwrap().count();
    ///////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // generate the new landscapes
    let new_landscapes = sample_landscapes - existing_landscapes;

    for i in 0..new_landscapes {
        let filename = format!(
            "landscapes/{}_{}{}.dat",
            model.short_description(), i + existing_landscapes, identifier
        );
        if Path::new(&filename).exists() {
            println!("Landscape '{}' already exists!", filename);
        } else {
            let fitness_landscape = FitnessLandscape::from_model(model);
            fitness_landscape.save(&filename);

            let filename = format!(
                "landscapes/{}_{}{}.bin.fl",
                model.short_description(), i + existing_landscapes, identifier
            );
            fitness_landscape.save_raw(&filename);
        }
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////

    if existing_landscapes > 0 {
        println!("{} landscapes existed matching the specified model", existing_landscapes);
    }
    if new_landscapes > 0 {
        println!("{} new landscapes generated (under the folder \"landscapes\")", new_landscapes);
    }
}
