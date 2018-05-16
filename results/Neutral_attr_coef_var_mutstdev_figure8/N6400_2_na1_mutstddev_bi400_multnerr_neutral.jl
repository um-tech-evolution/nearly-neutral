# Configuration for running spatial simulation
export simtype
# simtype=2 means spatial structure by changing the ideal values for attributes
@everywhere simtype = 2    
const num_trials = 1000
@everywhere const N_list = [6400,3200,1600,800,400,200,100,50,25,10,5,2]        # population size
#N_mut_list = [0.125, 0.25, 0.5, 1.0]
mutation_stddev_list = [0.002, 0.005,0.01,0.02,0.04]
const num_subpops = 1                    # Number of subpopulations
const num_attributes_list = [1,5]        # number attributes for quantitative representation
const ngens = 1				# Generations after burn-in
const burn_in= 400    # if an Int, generations of burn_in. If a Float, int_burn_in = burn_in*N+50  
const ideal = 1.0
const neutral = true
const additive_error=false
const wrap_attributes=false


