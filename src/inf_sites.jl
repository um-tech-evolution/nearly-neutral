export inf_sites
#=
A simplified version of the infinite sites model where sites evolve independently.
Called by nearly-neutral/src/run_infsites.jl
Author:  Alden H. Wright, Department of Computer Science, Univeristy of Montana, Missoula, MT 59812 USA
    alden.wright@umontana.edu
Example run:
[src]$ julia run.jl examples/is_example1
=#

@doc """ function inf_sites( )
Do a infinite sites simulation with popsize N, L loci, mutation rate mu per locus, ngens generations. 
The site update is done in function update_sites() in src/site_collection.jl.
Returns a infs_result_type object, but also prints relevant information.
"""
function inf_sites( N::Int64, N_mu::Float64, ngens::Int64; dfe::Function=dfe_neutral)
  fitness_table = Dict{Int64,Float64}()
  # Create a new empty site collection
  sc = site_collection( N, 1.0 )   # fix_minimum = 1.0 so fixes only if count is N
  g_limit = 100000   # Maximum of the number of extra generations to run while waiting for all sites to fix
  mu = N_mu/N
  new_id = 1
  i = 1
  N_sum = 0
  g_count = 0
  g = 1
  num_gens = 0
  done = false
  while !done && g < g_limit+ngens
    num_mutations = rand(Binomial(N,mu))
    if g <= ngens 
      for j in 1:num_mutations
        fit = dfe_fitness( new_id, dfe, fitness_table )  # Set fitness of new_id
        sc_push!(sc,site(new_id,N,g,fit))  
        #println("new site: id: ",new_id,"  start gen: ",sc.list[new_id].start_gen,"  sum_h: ",sc.list[new_id].sum_heteroz)
        new_id += 1
      end
      g_count += 1
    end
    update_sites!( sc, g, N )
    g += 1
    done = (g > ngens) && length(sc.active) == 0
  end
  #print_summary( sc )   # print a summary of the values stored in site_collection sc
  count_adv_del_fixed( sc )   # sets sc.count_fixed_adv and sc.count_fixed_del
  check_globals_against_site_collection( sc )    # Check that global statistics agree with those in sc.
  sc
end

@doc """ function dfe_fitness( p::Int64, dfe::Function, fitness_table::Dict{Int64,Float64} )
Fitness function as saved in dictionary fitness_table.
"""
function dfe_fitness( p::Int64, dfe::Function, fitness_table::Dict{Int64,Float64} )
  val = get(fitness_table,p,-1.0)
  if val == -1.0   # p is not in fitness table
    val = dfe( p )
    fitness_table[p] = val
  end
  val
end

