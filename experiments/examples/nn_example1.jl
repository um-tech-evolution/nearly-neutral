# Configuration for nearly neutral simuation
# Example run from command line:  julia n_neutral.jl examples/nn_example1
function dfemod(x) 
  dfe_mod(x,fit_inc=1.2)
end
const nn_simtype = 1    # nearly neutral infinite alleles
const type_str = "MD"
@everywhere const N_list     = [80,160]     # sample size n list, popsize N may be larger
const N_mu_list= [2.0,4.0]   # Mutation rate as a multiple of 1.0/N
const L = 1           # number of loci
const ngens = 10      # Generations after burn-in
const fix_minimum = 0.45  # Fraction of population needed to "fix"
const burn_in= 2.0    # generations of burn_in as a multiple of 1/mu
dfe=dfemod
dfe_str = "dfe_mod fit_inc: 1.2"

