# Configuration for nearly neutral simuation
function dfe_mixd( x::Int64 )
  dfe_mixed( x::Int64; adv_probability=0.01, 
    alpha_disadv=dfe_disadv_alpha, beta_disadv=dfe_disadv_beta,
    alpha_adv=dfe_adv_alpha, beta_adv=dfe_adv_beta )
end

const nn_simtype = 2    # nearly neutral infinite sites
type_str = "M1"
const mu_list_flag = true
@everywhere const N_list     = [25,50,100,200,400,800,1600,3200,6400]     # popsize
const mu_list= [0.025,0.01,0.005,0.0025]   # Mutation rate 
const ngens = 500000      # Generations after burn-in
const burn_in= 2.0    # generations of burn_in as a multiple of 1/mu
dfe = dfe_mixd
dfe_str = "dfe_mixd advprob=0.01 alpha_adv=0.2 alpha_disadv=0.5 beta_adv=0.005 beta_disadv=0.5"
dfe_adv_prob=0.01
dfe_adv_alpha=0.5
dfe_disadv_alpha=0.2
dfe_adv_beta=0.01
dfe_disadv_beta=0.5

