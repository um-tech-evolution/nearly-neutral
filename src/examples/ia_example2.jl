# Configuration for nearly neutral simuation
# Example run from command line:  julia run.jl examples/ia_example2
function dfemod(x)
  dfe_mod(x,fit_inc=1.2)
end
const nn_simtype = 1    # nearly neutral infinite alleles
@everywhere const N_list     = [40,80]     # sample size n list, popsize N may be larger
const N_mu_list= [2.0,4.0]   # Mutation rate as a multiple of 1.0/N
const L = 1           # number of loci
const ngens = 10      # Generations after burn-in
const burn_in= 2.0    # generations of burn_in as a multiple of 1/mu
const use_poplist=false   # save a list of populations as a check of correctness
dfe=dfemod
dfe_str = "dfe_mod fit_inc: 1.2"

