# Example Configuration for nearly neutral infinite sites simuation
# Command line to run:  julia n_neutral.jl examples/in_example1
# Creates output file examples/in_example2.csv
function dfe_disadv( x::Int64 )
  dfe_deleterious( x::Int64; alpha=dfe_disadv_alpha, beta=dfe_disadv_beta )
end

const type_str = "D4"
const mu_list_flag = false
const nn_simtype = 2    # nearly neutral infinite alleles
@everywhere const N_list     = [100,200]   # sample size, popsize N=  popsize_multiplier*n
const N_mu_list= [1.0,2.0]   # Mutation rate as a multiple of 1.0/N
const L = 1           # number of loci # not used in infinite alleles
const ngens = 1000      # Generations after burn-in
const burn_in= 2.0    # generations of burn_in as a multiple of N
dfe = dfe_disadv
dfe_str = "dfe_disadv alpha_disadv=0.2 beta_disadv=1.0"
dfe_disadv_alpha=0.2
dfe_disadv_beta=1.0


