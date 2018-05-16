#= Top-llevel file to run nearly neutral simulation from the command line.
Examples: (run from directory evotech/nearly-neutral/src):
> julia run.jl examples/ia_example1
> julia run.jl examples/is_example1
=#
# The following are top-level commands run when the file is executed by julia.
include("dfe_functions.jl")   # These need to be defined to read in the infinite sites configuration files

if length(ARGS) == 0
  simname = "../experiments/examples/ia_example1"
else
  simname = ARGS[1]
  if length(ARGS) >= 2   # second command-line argument is random number seed
    seed = parse(Int,ARGS[2])
    println("random number seed: ",seed)
    srand(seed)
  end
end
println("reading configuration file: ","$(simname).jl")
include("$(simname).jl")
println("opening output CSV file: ","$(simname).csv")
stream = open("$(simname).csv","w")

if !isdefined(:mu_list_flag)
  mu_list_flag=false
end
if nn_simtype==1
  println("infinite alleles model")
  include("run_infalleles.jl")
  ia_run_trials( mu_list_flag )
elseif nn_simtype == 2
  println("infinite sites model")
  include("run_infsites.jl")
  is_run_trials( mu_list_flag )
end

