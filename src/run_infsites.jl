#=
Example run from command line reading from the parameter file  examples/is_example?.jl (where ? = 1, 2, 3):
>  julia run.jl examples/is_example1

This run produces the output CSV file:  examples/in_example?.csv

This file contains functions for the infinite sites model in support of the Nearly Neutral paper
  ??? by Cheyenne Laue, ???, and Alden Wright
Author:  Alden H. Wright, Department of Computer Science, Univeristy of Montana, Missoula, MT 59812 USA
    alden.wright@umontana.edu
=#
using Dates
include("../src/InfSites.jl")

mutable struct infs_result_type
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
    println("infinite alleles model")
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

function is_run_trials( mu_list_flag::Bool=false)
  println("stream: ",stream)
  trial = 1
  N = N_list[1]
  if !mu_list_flag
    tr = infs_result( nn_simtype, N, N_mu_list[1], ngens, dfe, dfe_str )
    is_writeheader(Base.stdout, N_list, tr, mu_list_flag=mu_list_flag )
    is_writeheader(stream, N_list, tr, mu_list_flag=mu_list_flag )
    for N_mu in N_mu_list
      for N in N_list
        tr = infs_result( nn_simtype, N, N_mu, ngens, dfe, dfe_str )
        is_run_trial( tr )
        is_writerow(Base.stdout, trial, tr )
        is_writerow(stream, trial, tr )
        Base.flush(stream)
        #print_infs_result( tr )
        Base.flush(Base.stdout)
        trial += 1
      end
    end
  else  # if mu_list_flag
    tr = infs_result( nn_simtype, N, N*mu_list[1], ngens, dfe, dfe_str )
    is_writeheader(Base.stdout, N_list, tr, mu_list_flag=mu_list_flag )
    is_writeheader(stream, N_list, tr, mu_list_flag=mu_list_flag )
    for mu in mu_list
      for N in N_list
        N_mu = N*mu
        tr = infs_result( nn_simtype, N, N_mu, ngens, dfe, dfe_str )
        is_run_trial( tr )
        is_writerow(Base.stdout, trial, tr, mu_list_flag=mu_list_flag )
        is_writerow(stream, trial, tr, mu_list_flag=mu_list_flag )
        Base.flush(stream)
        #print_infs_result( tr )
        Base.flush(Base.stdout)
        trial += 1
      end
    end
end  # if !mu_list_flag
end

function is_run_trial( tr::infs_result_type )
  if tr.nn_simtype == 2
    sc = InfSites.inf_sites( tr.N, tr.N_mu, tr.ngens, dfe=tr.dfe )
    tr.number_extinctions = length(sc.extinct)
    tr.number_fixations = length(sc.fixed)
    tr.count_fixed_del = sc.count_fixed_del
    tr.count_fixed_adv = sc.count_fixed_adv
    tr.average_time_to_extinction = InfSites.average_time_to_extinction(sc)
    tr.average_time_to_fixation = InfSites.average_time_to_fixation(sc)
    tr.average_fitness_extinct = InfSites.average_fitness_extinct(sc)
    tr.average_fitness_fixed = InfSites.average_fitness_fixed(sc)
    tr.average_fitness_all = InfSites.average_fitness_all(sc)
    (tr.average_sites_per_gen, tr.stderr_sites_per_gen) = InfSites.sites_per_gen(sc)
    (tr.average_heterozygosity_per_gen, tr.stderr_heterozygosity_per_gen) = InfSites.heterozygosity_per_gen(sc)
    return tr
  else
    println("nn_simtype ",nn_simtype," not implemented")
  end
end

@doc """ function writeheader(), slat_reps::Int64=100000 ) 
Write header line for output CSV file that corresponds to stream.
See the comments for run_simulation for a description of the parameters.
"""
function is_writeheader(stream::IO, N_list::Vector{Int64}, 
    tr::infs_result_type; mu_list_flag::Bool=false )
  dfe_params = [:dfe_adv_prob, :dfe_adv_alpha, :dfe_adv_theta, :dfe_disadv_prob, :dfe_disadv_alpha, :dfe_disadv_theta]
  param_strings = [
    "# $(string(Dates.today()))",
    (isdefined(Main,:seed) ? "# random number seed: $(seed)" : "# no random number seed defined"),
    #"# $((tr.nn_simtype==1) ? "infinite alleles model" :"infinite_sites_model")",
    "# infinite_sites_model",
    "# N_list=\"$(N_list)\"",
    (mu_list_flag ? "# mu_list=\"$(mu_list)\"" : "# N_mu_list=\"$(N_mu_list)\""),
    "# ngens=$(tr.ngens)",
    "# dfe=$(tr.dfe)",
    "# dfe_str=$(tr.dfe_str)"
    ]
  for dp in dfe_params
    if isdefined(Main,dp)
      Base.push!(param_strings,"# $(string(dp))=$(eval(dp))")
    end
  end
  write(stream, join(param_strings, "\n"), "\n")
  first_heads = [
    "type",
    #  "n", 
    "N", (mu_list_flag ? "mu" : "N_mu")]
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

function is_writerow(stream::IO, trial::Int64, tr::infs_result_type; mu_list_flag::Bool=false )
  first = Any[
    isdefined(Main,:type_str0) ? type_str : trial,
    #tr.n,           # sample size
    tr.N,           # popsize
    mu_list_flag ? tr.N_mu/tr.N : tr.N_mu,          
  ]
 if tr.nn_simtype == 2
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

