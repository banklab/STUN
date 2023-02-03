# STUN
Forward-time **S**imulation on **TU**nable fit**N**ess landscapes in recombining
populations

This software runs simulations of adaptive trajectories of populations over
arbitrary fitness landscapes.

If you use this software please cite our associated [manuscript](...)

## Compiling
We provide precompiled versions of the software for [linux, macos and
windows](https://www.dropbox.com/sh/a6hbk3m823h3v43/AAB5wGVw9RIQh4RenOjvKP8Ja?dl=0).

If you wish to compile it yourself
1. install rust if it is not installed (https://www.rust-lang.org/tools/install)
   Note: It is possible that you'll have to restart the shell after installing
   rust
2. clone this repository (or download and unzip it)
3. run `cargo build --release` in the folder of the program code. This creates
   a program called `stun` (that can be found at the folder `target/release/`).
   Internet access is required when program is compiled for the first time,
   because `cargo` should download all required dependencies.
4. (optional) run `cp target/release/stun .` to copy the executable to current
   path

## Example
As an example, to generate 10 rough mount fuji landscapes with 5 loci each and
mu_a = 0.01, sigma_a = 0, sigma_e = 0.05 you would run
```
./stun --generate_only --landscapes 100 RMF -L 5 -m 0.01 -S 0 -s 0.05
```
There are three parts to this command. First, `./stun` the executable. Next, the
general options, in this case `--generate_only` (to specify that the fitness
landscapes should be generated, but no simulation is to be run) and
`--landscapes 100` to define that 100 landscapes should be generated. At last
comes the model specification, `RMF -L 5 -m 0.01 -S 0 -s 0.05`. The model and
its options should always come after all remaining options.

To run a simulation over, e.g., 10 NK fitness landscapes (with N = 8, K = 3),
with 50 replicates, three recombination rates (0.001, 0.01, 0.1), a population
of 200 individuals, initially uniformly distributed over the genotype space, the
user would run the following command
```
./stun --landscapes 10 --replicates 50 --recombination_rates 0.001,0.01,0.1 --size 200 --uniform 0.5 NK -N 8 -K 3
```
or, more compactly,
```
./stun -l 10 -c 50 -r 0.001,0.01,0.1 -N 200 -u 0.5 NK -N 8 -K 3
```
Notice that recombination rates are separated by commas, not spaces. The results
will be stored in the folder `data/`.

## Options
For help on all options available run `./stun help`. More detailed help on all
the options and function of the program are provided in the [manual](https://www.dropbox.com/sh/a6hbk3m823h3v43/AAB5wGVw9RIQh4RenOjvKP8Ja?dl=0) of the
software. Code documentation can be built and open with the command
`cargo doc --open`. This is useful to explore the details of the implementation
or get information to modify or expand the code.

### General options
This options come immediately after the executable name.
* Population size. Option `--size <population_size>` or `-N <population_size>`
* Initial population. Two possible initial populations are available:
  - a population sampled from neutral site frequency spectrum (option
    `--neutralsfs` or `-n`)
  - a population where the alleles in each individual are drawn from a uniform
    distribution with probability p (option `--uniform <p>` or `-u <p>`)
* Number of fitness landscapes. Option `--landscapes <number_of_landscapes>` or
  `-l <number_of_landscapes>`
* Number of replicates. Option `--replicates <number_of_replicates>` or
  `-c <number_of_replicates>`
* Recombination rates. There are two ways of defining the recombination
  - uniform recombination map. Option `--recombination_rates <rate1,rate2,...>`
  or `-r <rate1,rate2,...>`, where rate1, rate2, ... are the recombination rates
  separated by commas
  - arbitrary recombination maps. Option
  `--recombination_map <recombination_map>`, where recombination_map specifies
  the crossover probability in between any two loci (e.g., for a landscape with
  5 loci this could be `--recombination_map 0.05,0.01,0.5,0.01`). Alternatively,
  the user can specify a a file with a list of recombination maps, one per line.
* Output configuration. Option `--output <output_configuration_file>` or
  `-o <output_configuration_file>`. Check the file
  `example_output_configuration.toml` for a detailed example, documenting all
  options available.
* Only generate fitness landscapes, run no simulations. Option `--generate_only`
* Custom identifier, appended to the end of the output file names. Option
  `--id <id>`

For detailed help refer to the [manual](https://www.dropbox.com/sh/a6hbk3m823h3v43/AAB5wGVw9RIQh4RenOjvKP8Ja?dl=0)
of the software.

### Fitness landscape models
This options define the model and come after all other options. For the
options available for each model run `./stun help <model>` (replace model
for the target model, e.g., to obtain information on house of cards model run
`./stun help HoC`). The models available are
* House of cards model. Option `HoC`
* Additive model. Option `additive`
* Rough Mount Fuji model. Option `RMF`
* NK model. Option `NK`
* Block model. Option `block`
* Custom model. Provided by the user, option `custom`

Again, the [manual](https://www.dropbox.com/sh/a6hbk3m823h3v43/AAB5wGVw9RIQh4RenOjvKP8Ja?dl=0)
contains information on all these models, including a short description of each,
bibliography and usage examples.

## Contact
This software was developed at the THEE lab, at University of Bern. You can find
updated contact information at our [website](https://banklab.github.io/people/).
