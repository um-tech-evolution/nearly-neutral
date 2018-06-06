# Configuration for nearly neutral simuation
function dfe_mixd( x::Int64 )
  dfe_mixed( x::Int64; adv_probability=0.01, 
    alpha_disadv=dfe_disadv_alpha, beta_disadv=dfe_disadv_beta,
    alpha_adv=dfe_adv_alpha, beta_adv=dfe_adv_beta )
end

const nn_simtype = 1    # nearly neutral infinite alleles
type_str = "M3"
const mu_list_flag = false
@everywhere const N_list     = [25,50,100,200,400,800,1600,3200,6400]     # popsize
const N_mu_list= [0.5,1.0,2.0,4.0]   # Mutation rate as a multiple of 1.0/N
const L = 1           # number of loci
const ngens = 400000      # Generations after burn-in
const fix_minimum = 0.45  # Fraction of population needed to "fix"
const burn_in= 2.0    # generations of burn_in as a multiple of 1/mu
dfe = dfe_mixd
dfe_str = "dfe_mixd advprob=0.01 alpha_adv=0.5 alpha_disadv=0.2 beta_adv=0.02 beta_disadv=0.5"
dfe_adv_prob=0.01
dfe_adv_alpha=1.0
dfe_disadv_alpha=0.2
dfe_adv_beta=0.02
dfe_disadv_beta=0.5

const use_poplist=false   # if true, save a list of populations as a check of correctness
