# Example Configuration for nearly neutral simuation
# Command line to run:  julia n_neutral.jl examples/nn_example
# Creates output file examples/nn_example.csv

const nn_simtype = 1    # nearly neutral infinite alleles
popsize_multiplier_list = [1, 2, 4]
@everywhere const N_list     = [100,200]   # sample size, popsize N=  popsize_multiplier*n
const N_mu_list= [1.0,2.0]   # Mutation rate as a multiple of 1.0/N
const L = 3           # number of loci # not used in infinite alleles
const ngens = 1000      # Generations after burn-in
const burn_in= 2.0    # generations of burn_in as a multiple of N
dfe=dfe_neutral
dfe_str = "dfe_neutral"

