# Front end for src/nn_poplist.jl
include("../src/InfAlleles.jl")
if length(ARGS) == 0
  simname = "../experiments/examples/nn_example1"
else
  simname = ARGS[1]
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


@doc """ type trial_result
  A universal trial result for all run_trials functions
  Constructors for specific cases are defined below.
"""
type trial_result
  nn_simtype::Int64
  n::Int64    # sample size, must be <= N
  N::Int64    # popsize
  N_mu::Float64  # N*mu, population mutation rate
  ngens::Int64
  burn_in::Float64
  dfe::Function
  dfe_str::AbstractString
  expected_richness::Float64  # sum_{i=0}^{n-1} theta/(theta+i) where theta = 2*N_mu.
  average_richness::Float64   # Average number of traits in populations
  expected_w_homoz::Float64
  w_homoz::Float64
  IQV::Float64
end

function add_expected_richness( tr::trial_result )
  sum = 0.0
  theta = 2*tr.N_mu
  for i = 0:(tr.N-1)
    sum += 1.0/(theta+i)
  end
  tr.expected_richness = theta*sum
end

# Constructor that sets the parameters
function trial_result( nn_symtype::Int64, n::Int64, N::Int64, N_mu::Float64, ngens::Int64,
    burn_in::Float64=2.0, dfe::Function=dfe_neutral, 
    dfe_str::AbstractString="neutral" )
  tr = trial_result( nn_simtype, n, N, N_mu, ngens, burn_in, dfe, dfe_str, 0.0, 0.0, 0.0, 0.0, 0.0 )
  add_expected_richness( tr )
  tr
end


function print_trial_result( tr::trial_result )
  if tr.nn_simtype == 1
    println("\ninfinite alleleles model")
  else
    println("nn_simtype: ", tr.nn_simtype)
  end
  println("n: ", tr.n)
  println("N: ", tr.N)
  println("N_mu: ", tr.N_mu)
  println("mu: ", tr.N_mu/tr.N)
  println("ngens: ", tr.ngens)
  println("burn_in: ", tr.burn_in)
  println("dfe: ", tr.dfe)
  println("dfe_str: ", tr.dfe_str)
  println("expected richness: ",tr.expected_richness)
  println("average richness: ",tr.average_richness)
  println("expected_w_homoz: ",tr.expected_w_homoz)
  println("w_homoz: ",tr.w_homoz)
  println("IQV: ",tr.IQV)
end

function run_trials(popsize_multiplier_list::Vector{Int64}=[1]; mu_list_flag::Bool=false)
  println("stream: ",stream)
  trial = 1
  N = N_list[1]
  n = Int(floor(N*(1//popsize_multiplier_list[1])))
  if !mu_list_flag
    tr = trial_result( nn_simtype, n, N, N_mu_list[1], ngens, burn_in, dfe, dfe_str )
    writeheader(stream, popsize_multiplier_list, N_list, N_mu_list, tr )
    for N_mu in N_mu_list
      for N in N_list
        for psize_m in popsize_multiplier_list
          n = Int(floor(N*(1//psize_m)))
          tr = trial_result( nn_simtype, n, N, N_mu, ngens, burn_in, dfe, dfe_str )
          run_trial( tr )
          writerow(stream, trial, tr )
          trial += 1
        end
      end
    end
  else  # if mu_list_flag
    tr = trial_result( nn_simtype, n, N, N*mu_list[1], ngens, burn_in, dfe, dfe_str )
    writeheader(stream, popsize_multiplier_list, N_list, mu_list, tr, mu_list_flag=mu_list_flag )
    for mu in mu_list
      for N in N_list
        for psize_m in popsize_multiplier_list
          n = Int(floor(N*(1//psize_m)))
          N_mu = N*mu
          tr = trial_result( nn_simtype, n, N, N_mu, ngens, burn_in, dfe, dfe_str )
          run_trial( tr )
          writerow(stream, trial, tr, mu_list_flag=mu_list_flag )
          trial += 1
        end
      end
    end
  end  # if !mu_list_flag
end

function run_trials_mu(popsize_multiplier_list::Vector{Int64}=[1])
  println("stream: ",stream)
  trial = 1
  N = N_list[1]
  n = Int(floor(N*(1//popsize_multiplier_list[1])))
  tr = trial_result( nn_simtype, n, N, N*mu_list[1], ngens, burn_in, dfe, dfe_str )
  writeheader(stream, popsize_multiplier_list, N_list, mu_list, tr )
  for mu in mu_list
    for N in N_list
      for psize_m in popsize_multiplier_list
        n = Int(floor(N*(1//psize_m)))
        N_mu = N*mu
        tr = trial_result( nn_simtype, n, N, N_mu, ngens, burn_in, dfe, dfe_str )
        run_trial( tr )
        writerow(stream, trial, tr )
        trial += 1
      end
    end
  end
end

function run_trial( tr::trial_result )
  if tr.nn_simtype == 1
    poplist = nn_poplist(tr.N,tr.N_mu,tr.ngens,tr.dfe,combine=false)
    if tr.n < tr.N
      poplist = map(x->sample_population(x,tr.n),poplist)
    end
    add_stats_to_trial_result!( tr, poplist )
    print_trial_result( tr )
    return tr
  else
    println("nn_simtype ",nn_simtype," not implemented")
  end
end

function add_stats_to_trial_result!( tr::trial_result, poplist::Vector{Population} )
  pcounts = map(pop_counts64,poplist)
  tr.average_richness = mean(map(x->length(x),pcounts))
  tr.expected_w_homoz = 1.0/(1.0+2.0*tr.N_mu)
  tr.w_homoz = mean(map(watterson_homozygosity,pcounts))
  tr.IQV = mean(map(IQV,pcounts))
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
    "expected_w_heteroz",
    "w_heteroz",
    "IQV"
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
      1.0-tr.expected_w_homoz,
      1.0-tr.w_homoz,
      tr.IQV
    ]
  end
  line = join( vcat( first, mid ), "," )
  write(stream, line, "\n")
end

if isdefined(:popsize_multiplier_list)
  run_trials(popsize_multiplier_list, mu_list_flag=mu_list_flag )
else
  run_trials( mu_list_flag=mu_list_flag)
end
