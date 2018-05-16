# Configuration for nearly neutral simuation
function dfe_disadv( x::Int64 )
  dfe_deleterious( x::Int64; alpha=dfe_disadv_alpha, theta=dfe_disadv_theta )
end

const nn_simtype = 2    # nearly neutral infinite sites
const type_str = "D2"
const mu_list_flag = true
@everywhere const N_list     = [25,50,100,200,400,800,1600,3200,6400]     # popsize
const mu_list= [0.025,0.01,0.005,0.0025]   # Mutation rate as a multiple of 1.0/N
const ngens = 500000      # Generations after burn-in
const burn_in= 2.0    # generations of burn_in as a multiple of 1/mu
dfe = dfe_disadv
dfe_str = "dfe_disadv alpha_disadv=0.2 theta_disadv=1.0"
dfe_disadv_alpha=0.2
dfe_disadv_theta=1.0

