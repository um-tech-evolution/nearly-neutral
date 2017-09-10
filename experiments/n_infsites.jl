#=
Example run from command line reading from the parameter file  examples/nn_example?.jl (where ? = 1, 2, 3):
$  julia n_infsites.jl examples/in_example?

This run produces the output CSV file:  examples/in_example?.csv
=#
# The following are top-level commands run when the file is executed by julia.
# Front-end for src/infsites.jl
include("../src/InfSites.jl")
if length(ARGS) == 0
  simname = "../experiments/configs/in_example1"
else
  simname = ARGS[1]
  if length(ARGS) >= 2   # second command-line argument is random number seed
    seed = parse(Int,ARGS[2])
    println("seed: ",seed)
    srand(seed)
  end
end
include("$(simname).jl")
println("simname: ",simname)
stream = open("$(simname).csv","w")
println("stream: ",stream)

current_dir = pwd()
date_string = "../data/"*Dates.format(now(),"mm_dd_yy")*"/"
println("date string: ",date_string)
if !isdefined(:mu_list_flag)
  mu_list_flag=false
end

# Note that there are more top-level commands at the end of the file

type infs_result_type
  nn_simtype::Int64
  N::Int64
  N_mu::Float64
  ngens::Int64
  dfe::Function
  dfe_str::AbstractString
  number_mutations::Int64
  number_extinctions::Int64
  number_fixations::Int64
  average_time_to_extinction::Float64
  average_time_to_fixation::Float64
  average_fitness_extinct::Float64
  average_fitness_fixed::Float64
  average_fitness_all::Float64
  stderr_sites_per_gen::Float64
  average_sites_per_gen::Float64
  average_heterozygosity_per_gen::Float64
  stderr_heterozygosity_per_gen::Float64
  count_fixed_del::Int64
  count_fixed_adv::Int64
end

function infs_result( nn_simtype::Int64, N::Int64, N_mu::Float64, ngens::Int64,
     dfe::Function=dfe_neutral, dfe_str::AbstractString="neutral" )
  tr = infs_result_type( nn_simtype, N, N_mu, ngens, dfe, dfe_str, 0, 0, 0, 0.0, 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0 )
  tr
end


function print_infs_result( tr::infs_result_type )
  if tr.nn_simtype == 1
    println("infinite alleleles model")
  elseif tr.nn_simtype == 2
    println("infinite sites model")
  else
    println("nn_simtype: ", tr.nn_simtype)
  end
  println("N: ", tr.N)
  if mu_list_flag
    println("mu: ", tr.N_mu/tr.N)
  else
    println("N_mu: ", tr.N_mu)
  end
  println("ngens: ", tr.ngens)
  println("dfe: ", tr.dfe)
  println("dfe_str: ", tr.dfe_str)
  println("number_extinctions: ", tr.number_extinctions)
  println("number_fixations: ", tr.number_fixations)
  println("count_del_fixed: ", tr.count_fixed_del)
  println("count_adv_fixed: ", tr.count_fixed_adv)
  println("average time to extinction: ", tr.average_time_to_extinction)
  println("average time to fixation: ", tr.average_time_to_fixation)
  println("average fitness_extinct: ", tr.average_fitness_extinct)
  println("average fitness_fixed: ", tr.average_fitness_fixed)
  println("average fitness_all: ", tr.average_fitness_all)
  println("avg per generation count: ",tr.average_sites_per_gen)
  println("stderr per generation count: ",tr.stderr_sites_per_gen)
  println("avg per generation heterozygosity: ",tr.average_heterozygosity_per_gen)
  println("stderr per generation heterozygosity: ",tr.stderr_heterozygosity_per_gen)
  println()
end

function run_trials( mu_list_flag::Bool=false)
  println("stream: ",stream)
  trial = 1
  N = N_list[1]
  #=
  if isdefined(:popsize_multiplier_list)
    n = Int(floor(N*(1//popsize_multiplier_list[1])))
  else
    n=N
    popsize_multiplier_list=[1]
  end
  =#
  if !mu_list_flag
    tr = infs_result( nn_simtype, N, N_mu_list[1], ngens, dfe, dfe_str )
    writeheader(stream, N_list, tr, mu_list_flag=mu_list_flag )
    for N_mu in N_mu_list
      for N in N_list
        tr = infs_result( nn_simtype, N, N_mu, ngens, dfe, dfe_str )
        run_trial( tr )
        writerow(stream, trial, tr )
        print_infs_result( tr )
        trial += 1
      end
    end
  else  # if mu_list_flag
    tr = infs_result( nn_simtype, N, N*mu_list[1], ngens, dfe, dfe_str )
    writeheader(stream, N_list, tr, mu_list_flag=mu_list_flag )
    for mu in mu_list
      for N in N_list
        N_mu = N*mu
        tr = infs_result( nn_simtype, N, N_mu, ngens, dfe, dfe_str )
        run_trial( tr )
        writerow(stream, trial, tr, mu_list_flag=mu_list_flag )
        print_infs_result( tr )
        trial += 1
      end
    end
end  # if !mu_list_flag
end

function run_trial( tr::infs_result_type )
  if tr.nn_simtype == 2
    ic = inf_sites( tr.N, tr.N_mu, tr.ngens, dfe=tr.dfe )
    tr.number_extinctions = length(ic.extinct)
    tr.number_fixations = length(ic.fixed)
    tr.count_fixed_del = ic.count_fixed_del
    tr.count_fixed_adv = ic.count_fixed_adv
    tr.average_time_to_extinction = average_time_to_extinction(ic)
    tr.average_time_to_fixation = average_time_to_fixation(ic)
    tr.average_fitness_extinct = average_fitness_extinct(ic)
    tr.average_fitness_fixed = average_fitness_fixed(ic)
    tr.average_fitness_all = average_fitness_all(ic)
    (tr.average_sites_per_gen, tr.stderr_sites_per_gen) = sites_per_gen(ic)
    (tr.average_heterozygosity_per_gen, tr.stderr_heterozygosity_per_gen) = heterozygosity_per_gen(ic)
    return tr
  else
    println("nn_simtype ",nn_simtype," not implemented")
  end
end

@doc """ function writeheader(), slat_reps::Int64=100000 ) 
Write header line for output CSV file that corresponds to stream.
See the comments for run_simulation for a description of the parameters.
"""
function writeheader(stream::IO, N_list::Vector{Int64}, 
    tr::infs_result_type; mu_list_flag::Bool=false )
  dfe_params = [:dfe_adv_prob, :dfe_adv_alpha, :dfe_adv_theta, :dfe_disadv_prob, :dfe_disadv_alpha, :dfe_disadv_theta]
  param_strings = [
    "# $(string(Dates.today()))",
    "# $((tr.nn_simtype==1)?"infinite alleles model":"infinite_sites_model")",
    "# N_list=\"$(N_list)\"",
    (mu_list_flag ? "# mu_list=\"$(mu_list)\"": "# N_mu_list=\"$(N_mu_list)\""),
    "# ngens=$(tr.ngens)",
    "# dfe=$(tr.dfe)",
    "# dfe_str=$(tr.dfe_str)"
    ]
  for dp in dfe_params
    if isdefined(dp)
      Base.push!(param_strings,"# $(string(dp))=$(eval(dp))")
    end
  end
  write(stream, join(param_strings, "\n"), "\n")
  first_heads = [
    "type",
    #  "n", 
    "N", (mu_list_flag ? "mu":"N_mu")]
  mid_heads = []
  last_heads =
  [ "num_extinct", 
    "num_fixed", 
    "num_adv_fixed",
    "num_del_fixed",
    "fraction_fixed",
    "ave_extinct_time", 
    "ave_fixed_time", 
    "ave_extinct_selcoef", 
    "ave_fixed_selcoef", 
    "ave_all_selcoef",
    "ave_innov_count",
    "stderr_innov_count",
    "ave_heteroz",
    "stderr_heteroz",
  ] 
  line = join(vcat( first_heads, mid_heads, last_heads), ",")
  write(stream, line, "\n")
end

function writerow(stream::IO, trial::Int64, tr::infs_result_type; mu_list_flag::Bool=false )
  first = Any[
    (isdefined(:type_str) ? type_str :trial),
    #tr.n,           # sample size
    tr.N,           # popsize
    mu_list_flag ? tr.N_mu/tr.N : tr.N_mu,          
  ]
  if tr.nn_simtype == 0
    mid = Any[]
  elseif tr.nn_simtype == 2
    mid = Any[
      tr.number_extinctions,
      tr.number_fixations,
      tr.count_fixed_adv,
      tr.count_fixed_del,
      Float64(tr.number_fixations)/(tr.number_extinctions+tr.number_fixations),
      tr.average_time_to_extinction,
      tr.average_time_to_fixation,
      tr.average_fitness_extinct,
      tr.average_fitness_fixed,
      tr.average_fitness_all,
      tr.average_sites_per_gen,
      tr.stderr_sites_per_gen,
      tr.average_heterozygosity_per_gen,
      tr.stderr_heterozygosity_per_gen
    ]
  else 
    println("Error:  nn_simtype must be 2 for inf sites.")
    mid = Any[]
  end
  line = join( vcat( first, mid ), "," )
  write(stream, line, "\n")
end

run_trials( mu_list_flag )
