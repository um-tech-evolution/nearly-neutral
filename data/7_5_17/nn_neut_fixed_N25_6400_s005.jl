# Configuration for nearly neutral simuation
const dfe_s_list = [0.0,-0.005]
dfe_str = "dfe_neut_fixed"
const nn_simtype = 1    # nearly neutral infinite alleles
const mu_list_flag = false
@everywhere const N_list     = [25,50,100,200,400,800,1600,3200,6400]     # popsize
const N_mu_list= [0.5,1.0,2.0,4.0]   # Mutation rate as a multiple of 1.0/N
const L = 1           # number of loci
const ngens = 200000      # Generations after burn-in
const fix_minimum = 0.45  # Fraction of population needed to "fix"
const burn_in= 2.0    # generations of burn_in as a multiple of 1/mu

