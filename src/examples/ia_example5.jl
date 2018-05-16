# Example Configuration for nearly neutral simuation
# Command line to run:  julia n_neutral.jl examples/nn_example
# Creates output file examples/nn_example.csv

const nn_simtype = 1    # nearly neutral infinite alleles
mu_list_flag =true
type_str="I1"
@everywhere const N_list     = [100,200]   # popsize N
const N_mu_list= [1.0,2.0]   # Mutation rate as a multiple of 1.0/N
const mu_list = [0.1, 0.02]
const L = 1           # number of loci # not used in infinite alleles
const ngens = 1000      # Generations after burn-in
const burn_in= 2.0    # generations of burn_in as a multiple of N
const use_poplist=false   # save a list of populations as a check of correctness
dfe=dfe_neutral
dfe_str = "dfe_neutral"

