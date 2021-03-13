#=
This file contains functions for the infinite alleles model in support of the Nearly Neutral paper.
Author:  Alden H. Wright, Department of Computer Science, Univeristy of Montana, Missoula, MT 59812 USA
    alden.wright@umontana.edu
Example run:
[experiments]$ julia run.jl examples/ia_example1
=#
using QuadGK
export fitness, inf_alleles, pop_counts64, 
    dfe_deleterious,dfe_advantageous, dfe_mixed, dfe_mod, dfe_neutral, dfe_fixed 

integrand(x,theta) = theta/x*(1-x)^(theta-1)

#=
function print_trial_result( tr::InfAlleles.trial_result )
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
  println("mean fitness: ",tr.mean_fitness)
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
=#

#=
function add_expected_richness( tr::trial_result )
  sum = 0.0
  theta = 2*tr.N_mu
  for i = 0:(tr.N-1)
    sum += 1.0/(theta+i)
  end
  tr.expected_richness = theta*sum
end
=#
function add_expected_richness( tr::trial_result )
  theta = 2*tr.N_mu
  println("N: ",tr.N,"  N_mu: ",tr.N_mu,"  integrand(1/tr.N,theta): ",integrand(1/tr.N,theta))
  result, err = quadgk( x->integrand(x,theta), 1.0/tr.N, 1.0 )
  try
    #println("N: ",tr.N,"  N_mu: ",tr.N_mu,"  integrand(1/tr.N,theta): ",integrand(1/tr.N,theta))
    result, err = quadgk( x->integrand(x,theta), 1.0/tr.N, 1.0 )
  catch
    println("computation of integrand failed. expected richness result will be incorrect.")
    result = 1.0
    err = 1.0
  end
  result += theta
  #println("add_expected richness result: ",result,"  err: ",err)
  tr.expected_richness = result
end


@doc """ function dfe_fitness( p::Int64, dfe::Function, fitness_table::Dict{Int64,Float64} )
Fitness function as saveed in dictionary fitness_table.
"""
function dfe_fitness( p::Int64, dfe::Function, fitness_table::Dict{Int64,Float64} )
  val = get(fitness_table,p,-1.0)
  if val == -1.0   # p is not in fitness table
    val = dfe( p )
    fitness_table[p] = val
  end
  val
end 

@doc """ function inf_alleles( )
This function does the iteration over generations for the infinite alleles model.
Note:  dfe is "distribution of fitness effects" function.  
In simulations for the nearly neutral paper, uniform_start and combine are both false.
"""
function inf_alleles( tr::trial_result; uniform_start::Bool=false, combine::Bool=false)
  global fitness_table = Dict{Int64,Float64}()
  g_limit = 1000000  # upper limit of generations to wait for extinctions and fixations
  int_burn_in = Int(round(tr.burn_in*tr.N/tr.N_mu+50.0))
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
  sum_mean_fitness = 0.0
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
        new_id += 1
      end
    end
    new_pop = propsel( new_pop, tr.dfe, fitness_table )
    if g > int_burn_in
      
      # Accumulate statistics
      #println("mean fitness: ", mean( [ fitness_table[p] for p in new_pop ] ) )
      sum_mean_fitness += mean( [ fitness_table[p] for p in new_pop ] ) 
      pcounts = pop_counts64(new_pop)
      IQVvalue = IQV( pcounts )
      sum_IQV += IQVvalue
      sum_sq_IQV += IQVvalue^2
      homoz = watterson_homozygosity( pcounts )
      sum_homoz += homoz
      sum_sq_homoz += homoz^2
      richness = length(pcounts)
      sum_richness += richness
      sum_sq_richness += richness^2
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
  end
  #println("gcount: ",gcount)
  tr.mean_fitness = sum_mean_fitness/gcount
  tr.expected_w_homoz = 1.0/(1.0+2.0*tr.N_mu)
  tr.w_homoz = sum_homoz/gcount
  tr.stderr_w_homoz = sqrt((sum_sq_homoz/(gcount-1) - sum_homoz^2/gcount/(gcount-1))/gcount)
  tr.average_richness = sum_richness/gcount
  tr.stderr_richness = sqrt((sum_sq_richness/(gcount-1) - sum_richness^2/gcount/(gcount-1))/gcount)
  tr.IQV = sum_IQV/gcount
  tr.stderr_IQV = sqrt((sum_sq_IQV/(gcount-1) - sum_IQV^2/gcount/(gcount-1))/gcount)
  if tr.use_poplist
    if combine
      return [pop_result]
    else
      return poplist[int_burn_in+1:int_burn_in+tr.ngens]
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
