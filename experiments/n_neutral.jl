# Front end for src/nearly_neutral_poplist.jl
# For the case where the selection coefficient or coefficients are fixed rather than drawn from a gamma distribtuion
include("../src/NeutralCulturalEvolution.jl")
if length(ARGS) == 0
  #simname = "../experiments/nn_configs/nn_example1"
  simname = "../experiments/configs/nn_example"
else
  simname = ARGS[1]
end
include("$(simname).jl")
println("simname: ",simname)
stream = open("$(simname).csv","w")
println("stream: ",stream)

current_dir = pwd()
date_string = "../data/"*Dates.format(now(),"mm_dd_yy")*"/"
try mkdir(date_string) catch end   # create today's directory with no error if it already exists
println("date string: ",date_string)
#Ylist = [2,4,6,8,10]
#println("Ylist: ",Ylist)
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
  fix_minimum::Float64
  dfe::Function
  dfe_str::AbstractString
  expected_richness::Float64  # sum_{i=0}^{n-1} theta/(theta+i) where theta = 2*N_mu.
  average_richness::Float64   # Average number of traits in populations
  expected_w_homoz::Float64
  w_homoz::Float64
  IQV::Float64
  number_active::Int64
  number_extinctions::Int64
  number_fixations::Int64
  extinction_time::Float64
  fixation_time::Float64
  average_fitness_fixed::Float64
  average_fitness_extinct::Float64
  average_fitness_all::Float64
  #t2::Float64  # Average turnover for y=2
  #t4::Float64  # Average turnover for y=4
  #t6::Float64  # Average turnover for y=6
  #t8::Float64  # Average turnover for y=8
  #t10::Float64  # Average turnover for y=10
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
    fix_minimum::Float64=0.5, burn_in::Float64=2.0, dfe::Function=dfe_neutral, 
    dfe_str::AbstractString="neutral" )
  #tr = trial_result( nn_simtype, n, N, N_mu, ngens, burn_in, fix_minimum, dfe, dfe_str, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0,0.0,0.0,0,0,0.0 )
  tr = trial_result( nn_simtype, n, N, N_mu, ngens, burn_in, fix_minimum, dfe, dfe_str, 0.0, 0.0, 0.0, 0.0, 0.0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0 )
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
  println("fix_minimum: ", tr.fix_minimum)
  println("dfe: ", tr.dfe)
  println("dfe_str: ", tr.dfe_str)
  println("expected richness: ",tr.expected_richness)
  println("average richness: ",tr.average_richness)
  println("expected_w_homoz: ",tr.expected_w_homoz)
  println("w_homoz: ",tr.w_homoz)
  println("IQV: ",tr.IQV)
  println("number_active: ", tr.number_active)
  println("number_extinctions: ", tr.number_extinctions)
  println("number_fixations: ", tr.number_fixations)
  println("extinction_time: ", tr.extinction_time)
  println("fixation_time: ", tr.fixation_time)
  println("average_fitness_fixed: ", tr.average_fitness_fixed)
  println("average_fitness_extinct: ", tr.average_fitness_extinct)
  println("average_fitness_all: ", tr.average_fitness_all)
  #println("turnover 2: ",tr.t2)  # These must be compatible with global Ylist
  #println("turnover 4: ",tr.t4)
  #println("turnover 6: ",tr.t6)
  #println("turnover 8: ",tr.t8)
  #println("turnover 10: ",tr.t10)
end

# Convert an innovation_oollection to a trial_result
# The innovation_connection and the trial_result must agree in fields n, N, 
function convert_to_trial_result( ic::innovation_collection, tr::trial_result )
  tr.number_active = length(ic.active)
  tr.number_extinctions = length(ic.extinct)
  tr.number_fixations = length(ic.fixed)
  tr.fix_minimum = ic.fix_minimum
  tr.extinction_time = average_time_to_extinction( ic )
  tr.fixation_time = average_time_to_fixation( ic )
  tr.average_fitness_fixed = average_fitness_fixed( ic::innovation_collection )
  tr.average_fitness_extinct = average_fitness_extinct( ic::innovation_collection )
  tr.average_fitness_all = average_fitness_all( ic::innovation_collection )
end

function run_trials(popsize_multiplier_list::Vector{Int64}=[1]; mu_list_flag::Bool=false)
  println("stream: ",stream)
  trial = 1
  N = N_list[1]
  n = Int(floor(N*(1//popsize_multiplier_list[1])))
  if !mu_list_flag
    tr = trial_result( nn_simtype, n, N, N_mu_list[1], ngens, fix_minimum, burn_in, dfe, dfe_str )
    writeheader(stream, popsize_multiplier_list, N_list, N_mu_list, tr )
    for N_mu in N_mu_list
      for N in N_list
        for psize_m in popsize_multiplier_list
          n = Int(floor(N*(1//psize_m)))
          tr = trial_result( nn_simtype, n, N, N_mu, ngens, fix_minimum, burn_in, dfe, dfe_str )
          run_trial( tr )
          writerow(stream, trial, tr )
          trial += 1
        end
      end
    end
  else  # if mu_list_flag
    tr = trial_result( nn_simtype, n, N, N*mu_list[1], ngens, fix_minimum, burn_in, dfe, dfe_str )
    writeheader(stream, popsize_multiplier_list, N_list, mu_list, tr, mu_list_flag=mu_list_flag )
    for mu in mu_list
      for N in N_list
        for psize_m in popsize_multiplier_list
          n = Int(floor(N*(1//psize_m)))
          N_mu = N*mu
          tr = trial_result( nn_simtype, n, N, N_mu, ngens, fix_minimum, burn_in, dfe, dfe_str )
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
  tr = trial_result( nn_simtype, n, N, N*mu_list[1], ngens, fix_minimum, burn_in, dfe, dfe_str )
  writeheader(stream, popsize_multiplier_list, N_list, mu_list, tr )
  for mu in mu_list
    for N in N_list
      for psize_m in popsize_multiplier_list
        n = Int(floor(N*(1//psize_m)))
        N_mu = N*mu
        tr = trial_result( nn_simtype, n, N, N_mu, ngens, fix_minimum, burn_in, dfe, dfe_str )
        run_trial( tr )
        writerow(stream, trial, tr )
        trial += 1
      end
    end
  end
end

function run_trial( tr::trial_result )
  if tr.nn_simtype == 1
    ic = innovation_collection( tr.N, tr.fix_minimum )
    poplist = nearly_neutral_poplist(tr.N,tr.N_mu,tr.ngens,tr.dfe,combine=false,ic=ic)
    convert_to_trial_result( ic, tr )
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
  #turnover_list = turnover(poplist,Ylist,tr.N_mu)
  #=
  len = minimum(map(length,turnover_list))
  sum_list = zeros(Int64,len)
  for t in turnover_list
    for i = 1:len
      sum_list[i] += t[i]
    end
  end
  =#
  #tr.t2 = len >= 1 ? Float64(sum_list[1])/length(turnover_list) : -1.0
  #tr.t4 = len >= 2 ? Float64(sum_list[2])/length(turnover_list) : -1.0
  #tr.t6 = len >= 3 ? Float64(sum_list[3])/length(turnover_list) : -1.0
  #tr.t8 = len >= 4 ? Float64(sum_list[4])/length(turnover_list) : -1.0
  #tr.t10 = len >= 5 ? Float64(sum_list[5])/length(turnover_list) : -1.0
end

@doc """ function writeheader()
    burn_in::Float64, slat_reps::Int64=100000 ) 
Write header line for output CSV file that corresponds to stream.
See the comments for run_simulation for a description of the parameters.
"""
function writeheader(stream::IO, popsize_multiplier_list::Vector{Int64}, N_list::Vector{Int64}, 
    mu_list::Vector{Float64}, tr::trial_result; mu_list_flag::Bool=false )
  dfe_params = [:dfe_adv_prob, :dfe_adv_alpha, :dfe_adv_theta, :dfe_disadv_prob, :dfe_disadv_alpha, :dfe_disadv_theta]
  param_strings = [
    "# $(string(Dates.today()))",
    "# $((nn_simtype==1)?"infinite alleles model":"infinite_sites_model")",
    "# popsize_multiplier_list=$(popsize_multiplier_list)",
    "# N_list=$(N_list)",
    mu_list_flag ? "# mu_list=$(mu_list)": # N_mu_list=$(N_mu_list)",
    "# ngens=$(tr.ngens)",
    #"# Ylist=$(Ylist)",
    "# burn_in=$(tr.burn_in)",
    "# fix_minimum=$(tr.fix_minimum)",
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
    "IQV",
    "num_extinct", 
    "num_fixed", 
    "fraction_fixed",
    "ave_extinct_time", 
    "ave_fixed_time", 
    "ave_extinct_selcoef", 
    "ave_fixed_selcoef", 
    "ave_all_selcoef",
    #"t2",   # These must be compatible with global Ylist
    #"t4",
    #"t6",
    #"t8",
    #"t10"
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
      tr.IQV,
      tr.number_extinctions,
      tr.number_fixations,
      Float64(tr.number_fixations)/(tr.number_extinctions+tr.number_fixations),
      tr.extinction_time,
      tr.fixation_time,
      tr.average_fitness_extinct,
      tr.average_fitness_fixed,
      tr.average_fitness_all,
      #tr.t2,   # These must be compatible with global Ylist
      #tr.t4,
      #tr.t6,
      #tr.t8,
      #tr.t10
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
