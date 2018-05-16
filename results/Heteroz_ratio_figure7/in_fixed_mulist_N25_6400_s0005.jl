# Configuration for nearly neutral simuation
const dfe_s = -0.0005
function dfefixed( x::Int64 )
  dfe_fixed( x::Int64; s = dfe_s )
end

const nn_simtype = 2    # nearly neutral infinite sites
const type_str = "S6"
const mu_list_flag = true
@everywhere const N_list     = [25,50,100,200,400,800,1600,3200,6400]     # popsize
const mu_list= [0.01]   # Mutation rate 
const L = 1           # number of loci
const ngens = 500000      # Generations after burn-in
const burn_in= 2.0    # generations of burn_in as a multiple of 1/mu
dfe = dfefixed
dfe_str = "dfe_fixed s = $(dfe_s)"

