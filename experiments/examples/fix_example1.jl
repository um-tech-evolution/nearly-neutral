# Configuration for nearly neutral simuation with fixed selection coefficient
# Example run from command line:  julia n_neutral.jl examples/fix_example1
const nn_simtype = 1    # nearly neutral infinite alleles
const dfe_s_list = [0.0,0.1]
const dfe_modulus = 5 
dfe_str = "dfe_mod with modulus $dfe_modulus"
@everywhere const N_list     = [30]     # sample size n list, popsize N may be larger
const N_mu_list= [2.0]   # Mutation rate as a multiple of 1.0/N
const ngens = 3      # Generations after burn-in
const fix_minimum = 0.45  # Fraction of population needed to "fix"
const burn_in= 2.0    # generations of burn_in as a multiple of 1/mu

