# Configuration for nearly neutral simuation
function dfe_disadv( x::Int64 )
  dfe_deleterious( x::Int64; alpha=dfe_disadv_alpha, beta=dfe_disadv_beta )
end

const nn_simtype = 1    # nearly neutral infinite alleles
const type_str = "D4"
const mu_list_flag = false
@everywhere const N_list     = [25,50,100,200,400,800,1600,3200,6400]     # popsize
const N_mu_list= [0.5,1.0,2.0,4.0]   # Mutation rate as a multiple of 1.0/N
const L = 1           # number of loci
const ngens = 400000      # Generations after burn-in
const fix_minimum = 0.45  # Fraction of population needed to "fix"
const burn_in= 2.0    # generations of burn_in as a multiple of 1/mu
dfe = dfe_disadv
dfe_str = "dfe_disadv alpha_disadv=0.2 beta_disadv=1.0"
dfe_disadv_alpha=0.2
dfe_disadv_beta=1.0

const use_poplist=false   # if true, save a list of populations as a check of correctness
