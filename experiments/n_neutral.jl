#=
Example run from command line reading from the parameter file  examples/nn_example1.jl:
$  julia n_neutral.jl examples/nn_example1

This run produces the output CSV file:  examples/nn_example1.csv
=#
# Front end for src/nn_poplist.jl
# The following are top-level commands run when the file is executed by julia.
include("../src/InfAlleles.jl")
if length(ARGS) == 0
  simname = "../experiments/examples/nn_example1"
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
println("nn_simtype: ",nn_simtype)
stream = open("$(simname).csv","w")
println("stream: ",stream)

current_dir = pwd()
date_string = "../data/"*Dates.format(now(),"mm_dd_yy")*"/"
println("date string: ",date_string)
if !isdefined(:mu_list_flag)
  mu_list_flag=false
end

# Note that there are more top-level commands at the end of the file

@doc """ function run_trials()
  Run multiple trials of the infinite alleles model.
  popsize_multiplier_list  is a list of ratios of N to sample size n
          For the nearly neutral paper, sampling is not used, so it is always the case that n==N.
  There are 2 possibilities for iterating over mutation rate.
   (1) mu_list_flag==true:    Then  mu_list  is a list of mutation rates.  
   (2) mu_list_flag==false:   Then  N_mu_list  is a list of N_mu values.  N_mu is the product of N and mu.  This corresponds to constant theta.
"""
function run_trials(popsize_multiplier_list::Vector{Int64}=[1]; mu_list_flag::Bool=false)
  println("stream: ",stream)
  trial = 1
  N = N_list[1]
  n = Int(floor(N*(1//popsize_multiplier_list[1])))
  if !mu_list_flag
    tr = trial_result( nn_simtype, n, N, N_mu_list[1], ngens, burn_in, dfe, dfe_str, use_poplist=use_poplist )
    writeheader(stream, popsize_multiplier_list, N_list, N_mu_list, tr )
    for N_mu in N_mu_list
      for N in N_list
        for psize_m in popsize_multiplier_list
          n = Int(floor(N*(1//psize_m)))
          tr = trial_result( nn_simtype, n, N, N_mu, ngens, burn_in, dfe, dfe_str, use_poplist=use_poplist )
          run_trial( tr )
          writerow(stream, trial, tr )
          trial += 1
          Base.flush(stream)
        end
      end
    end
  else  # if mu_list_flag
    tr = trial_result( nn_simtype, n, N, N*mu_list[1], ngens, burn_in, dfe, dfe_str, use_poplist=use_poplist )
    writeheader(stream, popsize_multiplier_list, N_list, mu_list, tr, mu_list_flag=mu_list_flag )
    for mu in mu_list
      for N in N_list
        for psize_m in popsize_multiplier_list
          n = Int(floor(N*(1//psize_m)))
          N_mu = N*mu
          tr = trial_result( nn_simtype, n, N, N_mu, ngens, burn_in, dfe, dfe_str, use_poplist=use_poplist )
          run_trial( tr )
          writerow(stream, trial, tr, mu_list_flag=mu_list_flag )
          trial += 1
          Base.flush(stream)
        end
      end
    end
  end  # if !mu_list_flag
end

function run_trial( tr::trial_result )
  if tr.nn_simtype == 1
    if tr.use_poplist
      poplist = nn_poplist(tr,combine=false)
      if tr.n < tr.N
        poplist = map(x->sample_population(x,tr.n),poplist)
      end
      add_stats_to_trial_result!( tr, poplist )
    else
      nn_poplist(tr,combine=false)
    end
    print_trial_result( tr )
    Base.flush(STDOUT)
    return tr
  else
    println("nn_simtype ",nn_simtype," not implemented")
  end
end

function stddev(lst)
  N = length(lst)
  sqrt(mapreduce(x->x^2,+,lst)/(N-1) - sum(lst)^2/N/(N-1))
end

function add_stats_to_trial_result!( tr::trial_result, poplist::Vector{Population} )
  pcounts = map(pop_counts64,poplist)
  richness_list = map(x->length(x),pcounts)
  @assert isapprox( tr.average_richness, mean(richness_list) )
  tr.average_richness = mean(richness_list)
  @assert isapprox( tr.stderr_richness, stddev(richness_list)/sqrt(length(richness_list) ) ) 
  tr.stderr_richness = stddev(richness_list)/sqrt(length(richness_list))
  w_homoz_list = map(watterson_homozygosity,pcounts)
  @assert isapprox( tr.w_homoz, mean(w_homoz_list ) )
  tr.w_homoz = mean(w_homoz_list)
  tr.stderr_w_homoz = stddev(w_homoz_list)/sqrt(length(w_homoz_list) )
  IQV_list = map(IQV,pcounts)
  @assert isapprox( tr.IQV, mean(IQV_list) )
  tr.IQV = mean(IQV_list)
  @assert isapprox( tr.stderr_IQV, stddev(IQV_list)/sqrt(length(IQV_list)) )
  tr.stderr_IQV = stddev(IQV_list)/sqrt(length(IQV_list))
end

@doc """ function writeheader()
Write header line for output CSV file that corresponds to stream.
See the comments for run_simulation for a description of the parameters.
"""
function writeheader(stream::IO, popsize_multiplier_list::Vector{Int64}, N_list::Vector{Int64}, 
    mu_list::Vector{Float64}, tr::trial_result; mu_list_flag::Bool=false )
  dfe_params = [:dfe_adv_prob, :dfe_adv_alpha, :dfe_adv_beta, :dfe_disadv_prob, :dfe_disadv_alpha, :dfe_disadv_beta]
  param_strings = [
    "# $(string(Dates.today()))",
    "# $((nn_simtype==1)?"infinite alleles model":"infinite_sites_model")",
    "# popsize_multiplier_list=$(popsize_multiplier_list)",
    "# N_list=$(N_list)",
    mu_list_flag ? "# mu_list=$(mu_list)": # N_mu_list=$(N_mu_list)",
    "# ngens=$(tr.ngens)",
    "# burn_in=$(tr.burn_in)",
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
  last_heads =
  [ "expected_richness",
    "average_richness",
    "stderr_richness",
    "expected_w_heteroz",
    "w_heteroz",
    "stderr_heteroz",
    "IQV",
    "stderr_IQV"
  ] 
  line = join(vcat( first_heads, last_heads), ",")
  write(stream, line, "\n")
end

function writerow(stream::IO, trial::Int64, tr::trial_result; mu_list_flag::Bool=false   )
  first = Any[
    (isdefined(:type_str) ? type_str :trial),
    #tr.n,           # sample size
    tr.N,           # popsize
    mu_list_flag ? tr.N_mu/tr.N:tr.N_mu,          
  ]
  if tr.nn_simtype == 0
    mid = Any[]
  elseif tr.nn_simtype == 1
    mid = Any[
      tr.expected_richness,
      tr.average_richness,
      tr.stderr_richness,
      1.0-tr.expected_w_homoz,
      1.0-tr.w_homoz,
      tr.stderr_w_homoz,
      tr.IQV,
      tr.stderr_IQV
    ]
  end
  line = join( vcat( first, mid ), "," )
  write(stream, line, "\n")
end

# Top-level commands run when the file is loaded.
if isdefined(:popsize_multiplier_list)
  run_trials(popsize_multiplier_list, mu_list_flag=mu_list_flag )
else
  run_trials( mu_list_flag=mu_list_flag)
end
