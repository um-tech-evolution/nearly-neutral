export dfe_fitness, fitness, nn_poplist, pop_counts64, dfe_deleterious,dfe_advantageous, dfe_mixed, 
    dfe_mod, dfe_neutral, dfe_fixed

#using Distributions

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
function nn_poplist( N::Int64, N_mu::Float64, ngens::Int64, dfe::Function; burn_in::Float64=2.0, 
    uniform_start::Bool=false, nnselect::Int64=1, combine::Bool=true)
  global fitness_table = Dict{Int64,Float64}()
  g_limit = 1000000  # upper limit of generations to wait for extinctions and fixations
  int_burn_in = Int(round(burn_in*N/N_mu+50.0))
  mu = N_mu/N
  if uniform_start  # All allele values start with the same value.  Start with a selective sweep.
    poplist= Population[ Int64[1 for i = 1:N] ]
    fit = dfe_fitness(1, dfe, fitness_table )  # set fitness of 1 to be dfe(1).
    new_id = 2
  else
    poplist= Population[ collect(1:N) ]
    new_id = 1
    for i = 1:N
      # Note that dfe_advantageous(i), dfe_deleterious(i), dfe_mixed(i) do not depend on i, but dfe_mod(i) does depend on i
      fit = dfe_fitness(i, dfe, fitness_table )  # set fitness of i to be dfe(i). 
      new_id += 1
    end
  end
  if combine
    pop_result = Population()
  end
  done = false
  g = 2
  while !done && g < g_limit+ngens+burn_in
    for i in 1:N
      if rand() < mu
        poplist[g-1][i] = new_id
        fit = dfe_fitness( new_id, dfe, fitness_table )  # Set fitness of new_id
        if g > int_burn_in && g <= ngens+int_burn_in
          #println("id: ",new_id,"  fit: ",fit)
        end
        new_id += 1
      end
    end
    new_pop = propsel( poplist[g-1], dfe, fitness_table )
    Base.push!( poplist, new_pop )
    # TODO:  delete the next 3 lines and check that nothing is broken.
    #if g >int_burn_in && g <= ngens + int_burn_in
    #  u = unique(new_pop)
    #end
    if combine && g >= int_burn_in+1 && g <= ngens+int_burn_in
      pop_result = vcat( pop_result, new_pop )
    end
    g += 1
    done = (g > ngens+int_burn_in) 
  end
  if combine
    return [pop_result]
  else
    return poplist[int_burn_in+1:int_burn_in+ngens]
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
