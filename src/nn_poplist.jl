export trial_result, print_trial_result, dfe_fitness, fitness, nn_poplist, pop_counts64, 
    dfe_deleterious,dfe_advantageous, dfe_mixed, dfe_mod, dfe_neutral, dfe_fixed 

#using Distributions

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
  use_poplist::Bool
  expected_richness::Float64  # sum_{i=0}^{n-1} theta/(theta+i) where theta = 2*N_mu.
  average_richness::Float64   # Average number of traits in populations
  stderr_richness::Float64
  expected_w_homoz::Float64
  w_homoz::Float64
  stderr_w_homoz::Float64
  IQV::Float64
  stderr_IQV::Float64
end

# Constructor that sets the parameters
function trial_result( nn_simtype::Int64, n::Int64, N::Int64, N_mu::Float64, ngens::Int64,  
    burn_in::Float64=2.0, dfe::Function=dfe_neutral, 
    dfe_str::AbstractString="neutral"; use_poplist=false )
  tr = trial_result( nn_simtype, n, N, N_mu, ngens, burn_in, dfe, dfe_str, use_poplist, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 )
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
  println("stderr richness: ",tr.stderr_richness)
  println("expected_w_homoz: ",tr.expected_w_homoz)
  println("w_homoz: ",tr.w_homoz)
  println("stderr_w_homoz: ",tr.stderr_w_homoz)
  println("IQV: ",tr.IQV)
  println("stderr_IQV: ",tr.stderr_IQV)
  println("use_poplist: ",tr.use_poplist)
end

function add_expected_richness( tr::trial_result )
  sum = 0.0
  theta = 2*tr.N_mu
  for i = 0:(tr.N-1)
    sum += 1.0/(theta+i)
  end
  tr.expected_richness = theta*sum
end

@doc """ function dfe_fitness( p::Int64, dfe::Function, fitness_table::Dict{Int64,Float64} )
Fitness function as saveed in dictionary fitness_table.
"""
function dfe_fitness( p::Int64, dfe::Function, fitness_table::Dict{Int64,Float64} )
  val = get(fitness_table,p,-1.0)
  #println("val; ",val)
  if val == -1.0   # p is not in fitness table
    val = dfe( p )
    fitness_table[p] = val
  end
  #println("tbl: ",fitness_table)
  val
end 

@doc """ function nn_poplist( N::Int64, mu::Float64, ngens::Int64, dfe::Function; burn_in::Float64=1.0,
    uniform_start::Bool=false )
Note:  dfe is "distribution of fitness effects" function.  
"""
#function nn_poplist( N::Int64, N_mu::Float64, ngens::Int64, dfe::Function; burn_in::Float64=2.0, 
#    uniform_start::Bool=false, nnselect::Int64=1, combine::Bool=true)
function nn_poplist( tr::trial_result; uniform_start::Bool=false, nnselect::Int64=1, combine::Bool=true)
  global fitness_table = Dict{Int64,Float64}()
  #tr.use_poplist=false
  println("tr.use_poplist: ",tr.use_poplist)
  g_limit = 1000000  # upper limit of generations to wait for extinctions and fixations
  int_burn_in = Int(round(tr.burn_in*tr.N/tr.N_mu+50.0))
  #println("int_burn_in: ",int_burn_in)
  mu = tr.N_mu/tr.N
  if uniform_start  # All allele values start with the same value.  Start with a selective sweep.
    new_pop = Int64[1 for i = 1:tr.N] 
    fit = dfe_fitness(1, tr.dfe, fitness_table )  # set fitness of 1 to be tr.dfe(1).
    new_id = 2
  else
    new_pop = collect(1:tr.N)
    new_id = 1
    for i = 1:tr.N
      # Note that dfe_advantageous(i), dfe_deleterious(i), dfe_mixed(i) do not depend on i, but dfe_mod(i) does depend on i
      fit = dfe_fitness(i, tr.dfe, fitness_table )  # set fitness of i to be tr.dfe(i). 
      new_id += 1
    end
  end
  if tr.use_poplist
    poplist= Population[ new_pop ]
  end
  if combine
    pop_result = Population()
  end
  done = false
  # Initialize statistics accumulators
  sum_homoz = 0.0
  sum_sq_homoz = 0.0
  sum_IQV = 0.0
  sum_sq_IQV = 0.0
  sum_richness = 0.0
  sum_sq_richness = 0.0
  gcount = 0

  g = 2
  while !done && g < g_limit+tr.ngens+int_burn_in
    new_pop = deepcopy(new_pop)
    for i in 1:tr.N
      if rand() < mu
        new_pop[i] = new_id
        fit = dfe_fitness( new_id, tr.dfe, fitness_table )  # Set fitness of new_id
        if g > int_burn_in && g <= tr.ngens+int_burn_in
          #println("id: ",new_id,"  fit: ",fit)
        end
        new_id += 1
      end
    end
    new_pop = propsel( new_pop, tr.dfe, fitness_table )

    if g > int_burn_in
      # Accumulate statistics
      pcounts = pop_counts64(new_pop)
      IQVvalue = IQV( pcounts )
      sum_IQV += IQVvalue
      sum_sq_IQV += IQVvalue^2
      homoz = watterson_homozygosity( pcounts )
      #println("g: ",g,"  new_pop: ",new_pop,"  pcounts: ",pcounts,"  homoz: ",homoz)
      #println("length poplist: ",length(poplist))
      sum_homoz += homoz
      sum_sq_homoz += homoz^2
      richness = length(pcounts)
      sum_richness += richness
      sum_sq_richness += richness^2
      #println("g: ",g,"  richness: ",richness)
      gcount += 1
    end

    if tr.use_poplist
      Base.push!( poplist, new_pop )
    end
    if combine && g >= int_burn_in+1 && g <= tr.ngens+int_burn_in
      pop_result = vcat( pop_result, new_pop )
    end
    g += 1
    done = (g > tr.ngens+int_burn_in) 
    if g > int_burn_in
      #println("g: ",g,"  poplist: ", poplist[int_burn_in+1:end])
    end
  end
  #println("gcount: ",gcount)
  tr.expected_w_homoz = 1.0/(1.0+2.0*tr.N_mu)
  tr.w_homoz = sum_homoz/gcount
  tr.stderr_w_homoz = sqrt((sum_sq_homoz/(gcount-1) - sum_homoz^2/gcount/(gcount-1))/gcount)
  tr.average_richness = sum_richness/gcount
  #println("sum_richness: ",sum_richness)
  #println("sum_sq_richness: ",sum_sq_richness)
  tr.stderr_richness = sqrt((sum_sq_richness/(gcount-1) - sum_richness^2/gcount/(gcount-1))/gcount)
  #println("tr.stderr_richness: ", tr.stderr_richness)
  tr.IQV = sum_IQV/gcount
  tr.stderr_IQV = sqrt((sum_sq_IQV/(gcount-1) - sum_IQV^2/gcount/(gcount-1))/gcount)
  if tr.use_poplist
    if combine
      return [pop_result]
    else
      #println("int_burn_in+1: ",int_burn_in+1,"  int_burn_in+tr.ngens: ",int_burn_in+tr.ngens)
      #println("poplist: ", poplist[int_burn_in+1:int_burn_in+tr.ngens])
      return poplist[int_burn_in+1:int_burn_in+tr.ngens]
      #println("poplist: ", poplist)
      ##return poplist
    end
  end
end

@doc """ function pop_counts64( pop::Population )
Returns the sorted frequencies of the alleles of Population pop.
Example:  If pop = [5, 7, 9, 5, 4, 5, 7], then the returned list is [3, 2, 1, 1]
   because there are 3 5's, 2 7's, 1 9, and 1 4.  So the sum of the returned
   list is the length of the population.
"""
function pop_counts64( pop::Population )
  c = Dict{Int64,Int64}()
  for x in pop
    c[x] = get( c, x, 0 ) + 1
  end
  map( x->c[x], sort( unique(pop), by=x->c[x], rev=true ) )
end

function dfe_deleterious( x::Int64; alpha::Float64=0.2, beta::Float64=0.5 )
  global dist_deleterious
  if !isdefined(:dist_deleterious)
    dist_deleterious = Distributions.Gamma(alpha,beta)
  end
  return max(0.1,1.0-rand(dist_deleterious))
end
  
function dfe_advantageous( x::Int64; alpha::Float64=1.0, beta::Float64=0.5 )
  global dist_advantageous
  if !isdefined(:dist_advantageous)
    dist_advantageous = Distributions.Gamma(alpha,beta)
  end
  return 1.0+rand(dist_advantageous)
end

@doc """ function dfe_mixed( x::Int64; adv_probability::Float64=0.2, alpha_disadv::Float64=0.2, alpha_adv::Float64=1.0, beta_disadv::Float64=1.0, beta_adv::Float64=0.01 )
Mixed advantageous and deleterious.
"""

function dfe_mixed( x::Int64; adv_probability::Float64=0.2, alpha_disadv::Float64=0.2, alpha_adv::Float64=1.0, beta_disadv::Float64=1.0, beta_adv::Float64=0.05 )
  global dist_deleterious
  if !isdefined(:dist_deleterious)
    dist_deleterious = Distributions.Gamma(alpha_disadv,beta_disadv)
  end
  global dist_advantageous
  if !isdefined(:dist_advantageous)
    dist_advantageous = Distributions.Gamma(alpha_adv,beta_adv)
  end
  rand() < adv_probability ? 1.0+rand(dist_advantageous) : max(0.01,1.0-rand(dist_deleterious))
end

@doc """ function dfe_mod( )
Increased fitness for every x which is divisible by modulaus.
The increase in fitness is by a multiple of fit_inc.
"""
function dfe_mod( x::Int64; modulus::Int64=5, fit_inc::Float64=1.2 )
  if x % modulus == 0
	  return fit_inc
  else
    return 1.0  # flat fitness for no selection
  end
end

@doc """function dfe_neutral( x::Int64 )
"""
function dfe_neutral( x::Int64 )
  return 1.0   # flat fitness for no selection
end

@doc """function dfe_fixed( )
"""
function dfe_fixed( x::Int64; s::Float64=0.0 )
  return 1.0 + s   # fitness corresponding to selection coefficient s
end

#=  Tests
nn_select = 1
for i = 1:10 println(dfe_fitness(i,dfe_mod)) end
for i = 1:10 println(dfe_fitness(i,dfe_mixed)) end
N=10; N_mu =2.0; ngens = 3;
nn_poplist(N,N_mu,ngens,dfe_mod)
for i = 1:10 println(dfe_fitness(i,x->dfe_mod(x,fit_inc=2.0))) end
nn_poplist(N,N_mu,ngens,x->dfe_mod(x,fit_inc=2.0))
=#
