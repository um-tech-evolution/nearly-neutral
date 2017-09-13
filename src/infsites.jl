export inf_sites
#=

A simplified version of the infinite sites model where sites evolve independently.
Called by experiments/n_infsites.jl
Author:  Alden H. Wright, Department of Computer Science, Univeristy of Montana, Missoula, MT 59812 USA
    alden.wright@umontana.edu
Called by experiments/n_neutral.jl.
Example run:
[experiments]$ julia n_infsites.jl examples/in_example1
=#

@doc """ function inf_sites( )
Do a infinite sites simulation with popsize N, L loci, mutation rate mu per locus, ngens generations. 
The site update is done in function update_sites() in src/site_collection.jl.
Returns a infs_result_type object, but also prints relevant information.
"""
function inf_sites( N::Int64, N_mu::Float64, ngens::Int64; 
    dfe::Function=dfe_neutral)
  global sum_counts = 0
  global sum_sq_counts = 0
  global sum_heteroz = 0.0
  global sum_sq_heteroz = 0.0
  global sum_gens = 0
  fitness_table = Dict{Int64,Float64}()
  # Create a new empty site collection
  sc = site_collection( N, 1.0 )   # fix_minimum = 1.0 so fixes only if count is N
  #println("inf sum is.sum_generations: ",sc.sum_generations)
  g_limit = 100000   # Maximum of the number of extra generations to run while waiting for all sites to fix
  mu = N_mu/N
  println("N: ",N,"  mu: ",mu,)
  #println("N_mu: ",N_mu)
  #println("ngens: ",ngens)
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
  #println("g_count: ",g_count)
  print_summary(sc)
  #println("inf sum_counts: ",sum_counts)
  #println("inf sum_sq_counts: ",sum_sq_counts)
  #println("inf sum_heteroz: ",sum_heteroz)
  #println("inf sum_sq_heteroz: ",sum_sq_heteroz)
  #println("inf sum_gens: ",sum_gens)
  counts_var = sum_sq_counts/(sum_gens-1) - sum_counts^2/sum_gens/(sum_gens-1)
  heteroz_var = sum_sq_heteroz/(sum_gens-1) - sum_heteroz^2/sum_gens/(sum_gens-1)
  #=
  if counts_var >= 0.0
    println("inf counts_var: ",sqrt(counts_var)/sum_gens) 
  else
    println("inf counts_var: ",-1)
  end
  if heteroz_var >= 0.0
    println("inf heteroz_var: ",sqrt(heteroz_var)/sum_gens) 
  else
    println("inf heteroz_var: ",-1)
  end
  =#
  (sc.count_fixed_adv, sc.count_fixed_del) = count_adv_del_fixed( sc )
  sc
end

