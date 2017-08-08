# Configuration for nearly neutral simuation
# Example run from command line:  julia n_neutral.jl nn_examples/nn_example1
function dfemod(x)
  dfe_mod(x,fit_inc=1.2)
end
const nn_simtype = 1    # nearly neutral infinite alleles
const popsize_multiplier_list = [1,4]  # popsize N is this multiple of sample size n
@everywhere const N_list     = [40,80]     # sample size n list, popsize N may be larger
const N_mu_list= [2.0,4.0]   # Mutation rate as a multiple of 1.0/N
const L = 1           # number of loci
const ngens = 10      # Generations after burn-in
const burn_in= 2.0    # generations of burn_in as a multiple of 1/mu
dfe=dfemod
dfe_str = "dfe_mod fit_inc: 1.2"

