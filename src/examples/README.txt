This folder contains example configuration files with .jl extensions, and the resuling .csv output files.
Filenames starting with ia_ are infinite alleles, while filenames starting with is_ are infinite sites.

To run configuration file ia_example1.jl with random number seed 1 from the command line:

[src]: julia run.jl examples/ia_example1 1

If Julia's random number generation has not changed, then the .csv files here should agree with
the output from running the corresponding configuration files.
