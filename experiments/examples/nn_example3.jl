# Example Configuration for nearly neutral simuation
# Command line to run:  julia n_neutral.jl examples/nn_example3
# Creates output file examples/nn_example3.csv

function dfe_mixd( x::Int64 )
  dfe_mixed( x::Int64; adv_probability=0.01,
    alpha_disadv=dfe_disadv_alpha, beta_disadv=dfe_disadv_beta,
    alpha_adv=dfe_adv_alpha, beta_adv=dfe_adv_beta )
end
const nn_simtype = 1    # nearly neutral infinite alleles
mu_list_flag =true
type_str="I1"
@everywhere const N_list     = [100,200]   # sample size, popsize N=  popsize_multiplier*n
const N_mu_list= [1.0,2.0]   # Mutation rate as a multiple of 1.0/N
const mu_list = [0.1, 0.02]
const L = 1           # number of loci # not used in infinite alleles
const ngens = 1000      # Generations after burn-in
const burn_in= 2.0    # generations of burn_in as a multiple of N
const use_poplist=false  # set to false for "production" runs with large ngens
dfe=dfe_mixd
dfe_adv_prob=0.01
dfe_adv_alpha=1.0
dfe_disadv_alpha=0.2
dfe_adv_beta=0.005
dfe_disadv_beta=0.5
dfe_str = "dfe_mixd advprob=0.01 alpha_adv=1.0 alpha_disadv=0.2 beta_adv=0.005 beta_disadv=0.5"

