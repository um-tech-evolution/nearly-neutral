# Configuration for nearly neutral simuation
const nn_simtype = 1    # nearly neutral infinite alleles
const type_str = "N3"
const mu_list_flag = true
@everywhere const N_list     = [25,50,100,200,400,800,1600,3200,6400]     # popsize
const N_mu_list= [0.5,1.0,2.0,4.0]   # Mutation rate as a multiple of 1.0/N
const mu_list= [0.025,0.01,0.005,0.0025]   # Mutation rate
const L = 1           # number of loci
const ngens = 10000      # Generations after burn-in
const fix_minimum = 0.45  # Fraction of population needed to "fix"-
const burn_in= 4.0    # generations of burn_in as a multiple of 1/mu
dfe = dfe_neutral
dfe_str = "dfe_neutral"

const use_poplist=false   # if true, save a list of populations as a check of correctness
