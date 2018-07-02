# NearlyNeutral.jl

Julia source code for the simulations in support of a paper on the applications of
the nearly neutral theory from population genetics to cultural evolution.  The simulation
code is in Julia, while code for producing the plots is written in R.
As of June 2018, the code is almost compatible with Julia v. 7 with only a few 
deprecation warnings.

This paper is in progress, and does not yet have a definite title.


## Authors
  + Alden Wright (<alden.wright@umontana.edu>)  primary author of the code.
  + Cheyenne Laue (<cheyenne.laue@umontane.edu>)  first author of the paper.

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

This respository includes code for two models:  infinite alleles, infinite sites. 

To run Julia examples from the command line:
>  cd src
>  julia run.jl examples/ia\_example1    # run infinite alleles model
>  julia run.jl examples/is\_example1    # run infinite sites model
>  julia run.jl examples/ia\_example1 1  # run infinite alleles model using random number seed 1
>  julia run.jl examples/is\_example1 1  # run infinite sites model using random number seed 1

The file run\_tests.sh is a bash script to run all of the examples.  The results can be compared
with the .csv files in  examples/results/.

The top-level results folder includes the R code for producing the plots along with the necessary data.
