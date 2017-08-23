# NearlyNeutral.jl

Julia source code for the simulations in support of a paper on the applications of
the nearly neutral theory from population genetics to cultural evolution.  This 
paper is in progress, and does not yet have a defintie title.

## Authors
  + Alden Wright (<alden.wright@umontana.edu>)  primary author of the code.
  + Cheyenne Laue (<cheyenne.laue@umontane.edu>)  primary author of the paper.

## About nearly neutral theory

From [Wikipedia](https://en.wikipedia.org/Nearly_neutral_theory_of_molecular_evolution)

    The nearly neutral theory of molecular evolution is a modification of the neutral theory 
    of molecular evolution that accounts for the fact that not all mutations are either so 
    deleterious such that they can be ignored, or else neutral. Slightly deleterious mutations 
    are reliably purged only when their selection coefficient is greater than one divided by 
    the effective population size. In larger populations, a higher proportion of mutations 
    exceed this threshold for which genetic drift cannot overpower selection, leading to fewer 
    fixation events and so slower molecular evolution.
    
    The nearly neutral theory was proposed by Tomoko Ohta in 1973.[1] 

There two models:  infinite alleles, infinite sites. 

The top-level program for infinite alleles is experiments/n_neutral.jl. 
The top-level program for infinite sites is experiments/n_infsites.jl. 


