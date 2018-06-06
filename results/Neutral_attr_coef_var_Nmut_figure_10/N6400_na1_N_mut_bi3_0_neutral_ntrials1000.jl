# Configuration for running spatial simulation
export simtype
@everywhere simtype = 2    
const num_trials = 1000
@everywhere const N_list = [6400,3200,1600,800,400,200,100,50,25]        # Meta-population size
N_mut_list = [0.125, 0.25, 0.5, 1.0]
#mutation_stddev_list = [0.005,0.01,0.02,0.04]
const num_subpops = 1                    # Number of subpopulations
const num_attributes_list = [1]        # number attributes for quantitative representation
const ngens = 1       # Generations after burn-in
const burn_in= 3.0    # generations of burn_in as a multiple of N
const ideal = 1.0
const neutral = true
const fit_max=1.0
const fit_slope=1.0
const additive_error=false
const wrap_attributes=false


