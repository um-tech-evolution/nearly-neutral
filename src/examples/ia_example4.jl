# Configuration for nearly neutral simuation
function dfe_mixd( x::Int64 )
  dfe_mixed( x::Int64; adv_probability=0.01, 
    alpha_disadv=dfe_disadv_alpha, beta_disadv=dfe_disadv_beta,
    alpha_adv=dfe_adv_alpha, beta_adv=dfe_adv_beta )
end

const nn_simtype = 1    # nearly neutral infinite alleles
type_str = "M2"
const mu_list_flag = false
const N_list     = [25,100,400]     # popsize
const N_mu_list= [1.0,2.0]   # Mutation rate as a multiple of 1.0/N
const L = 1           # number of loci
const ngens = 40000      # Generations after burn-in
const burn_in= 2.0    # generations of burn_in as a multiple of 1/mu
const use_poplist=true   # Using poplist does additional correctness checking, but increases memory usage
dfe = dfe_mixd
dfe_str = "dfe_mixd advprob=0.01 alpha_adv=1.0 alpha_disadv=0.2 beta_adv=0.01 beta_disadv=0.5"
dfe_adv_prob=0.01
dfe_adv_alpha=1.0
dfe_disadv_alpha=0.2
dfe_adv_beta=0.01
dfe_disadv_beta=0.5
