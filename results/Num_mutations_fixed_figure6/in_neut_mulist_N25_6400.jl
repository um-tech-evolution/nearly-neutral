# Configuration for nearly neutral simuation
const nn_simtype = 2    # nearly neutral infinite sites
const type_str = "N2"
const mu_list_flag = true
@everywhere const N_list     = [25, 50, 100,200,400,800,1600, 3200, 6400]     # popsize
const mu_list= [0.025,0.01,0.005,0.0025]   # Mutation rate as a multiple of 1.0/N
const L = 1           # number of loci
const ngens = 500000      # Generations after burn-in
const burn_in= 2.0    # generations of burn_in as a multiple of 1/mu
dfe = dfe_neutral
dfe_str = "dfe_neutral" 

