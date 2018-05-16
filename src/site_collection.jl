#=  Type definitions for  site  and  site_collection  and related functions for the infinite sites model.
A site stores statistics for a single site, while a site_collection stores statistics for all of the sites of a trial.
The update_sites1() function is called on every generation to update all sites in the site collection,

GLOBAL VARIABLES:  The global variables declared in this file are used only for correctness checking via the @assert macro.
  There are equivalent fields in site_collection, and results depend only on these fields, not on global variables.
  Thus, all global variables and global variable computations can be removed without affecting functionality.
=#
export site_type, site,  
    site_collection, sc_push!, update_sites!, update_neutral, update_selected,
    print_summary, average_time_to_extinction, average_time_to_fixation, 
    fixed_fraction, average_fitness_fixed, average_fitness_extinct, average_fitness_all,
    sites_per_gen, heterozygosity_per_gen

#=
Stores the properties of an "site" (possibly deleterious or advantageous).
=#

type site_type
  identifier::Int64   # integer for this site, not sure if this is needed
  start_gen::Int64    # the generation (time step) when the site was generated
  final_gen::Int64    # the generation when the site went extinct.  Should be zero while site is evolving
  fitness_coefficient::Float64    # The selection coefficient of the site (allele)
  sum_gens::Int64     # sum of the number of generations where counts and heterozygosity are accumulated
  sum_counts::Int64   # the current accumulated sum of the counts per generation for infinite sites
  sum_sq_counts::Int64   # the current accumulated sum of squares of the counts per generation for infinite sites
  sum_heteroz::Float64  # the current accumulated sum of heterozygosities for infinite sites
  sum_sq_heteroz::Float64  # the current accumulated sum of squares ofheterozygosities for infinite sites
  previous_allele_freq::Int64  # Allele frequence of the previous generation
end

@doc """ function site() 
  Initializes the properties of the site.
  Initial values for sum_counts, sum_sq_counts, sum_heteroz, sum_sq_heteroz are added.  
  This means that the number of generations is one larger than the number of increments in update_innovations!(),
    and 1 is added to sum_gens in sites_per_gen() and heterozygosity_per_gen().
"""
function site( id::Int64, N::Int64, start_gen::Int64, fitness_coef::Float64=1.0 )
  # See the comments on global variables at the beginning of this file.
  global sum_counts += 1
  global sum_sq_counts += 1
  initial_heteroz = 1.0 - watterson_homozygosity([N-1,1])
  global sum_heteroz += initial_heteroz
  global sum_sq_heteroz += initial_heteroz^2
  global sum_gens += 1
  return site_type( id, start_gen, 0, fitness_coef, 1, 1, 1, initial_heteroz, initial_heteroz^2, 1  )
end

#=
Stores a collection of sites.
Sites are partitioned into 3 subsets: active, fixed, and extinct.
=#
type site_collection
  N::Int64                # popsize
  list::Dict{Int64,site_type} #    Stores a collection of sites
  active::IntSet   # indices of actively evolving sites
  fixed::IntSet   # indices of sites that have fixed
  extinct::IntSet   # indices of sites that have gone extinct
  fix_minimum::Float64  # minimum fraction of popsize for fixation for infinite alleles
  sum_gens::Int64     # sum of counts of generations
  sum_counts::Int64   # sum of counts of mutant alleles for infinite sites
  sum_sq_counts::Float64   # sum of counts of mutant alleles for infinite sites
  sum_heteroz::Float64  # sum of heterozygosities for infinite sites
  sum_sq_heteroz::Float64  # sum of heterozygosities squared for infinite sites
  sum_generations::Int64  # Total number of generations corresponding to sum_counts and sum_heteroz
  count_fixed_del::Int64  # Number of fixed mutations that are deleterious (fit_coef < 1.0)
  count_fixed_adv::Int64  # Number of fixed mutations that are advantageous (fit_coef > 1.0)
  in_use::Bool      # If false, not used
end


# Constructor for a new empty site collection with a value for fix_minimum
function site_collection( N::Int64, fix_min::Float64, in_use::Bool=true )  
  global sum_counts = 0
  global sum_sq_counts = 0
  global sum_heteroz = 0.0
  global sum_sq_heteroz = 0.0
  global sum_gens = 0
  site_collection( N, Dict{Int64,site_type}(), IntSet(), IntSet(), IntSet(), fix_min, 0, 0, 0.0, 0.0, 0.0, 0,0,0, in_use )
end

@doc """ sc_push!()
 Adds an site to the collection.
"""
function sc_push!( innov_collection::site_collection, innov::site_type )
  if !innov_collection.in_use 
    return
  end
  innov_collection.list[innov.identifier] = innov 
  if innov.final_gen == 0
    Base.push!( innov_collection.active, innov.identifier )
  elseif innov.history[end] > 0
    Base.push!( innov_collection.fixed, innov.identifier )
  else
    Base.push!( innov_collection.extinct, innov.identifier )
  end
  innov_collection
end

@doc """ fix_test()
The most basic test for fixation, namely that the allele frequency is the population size.
This is used in the infitite sites model.
"""
function fix_test( N::Int64, new_allele_freq::Int64 )
  return N == new_allele_freq
end

@doc """ function update_sites!( )
 Update all active sites, and make some extinct and fixed
 This version is used only by infsites.jl.
"""
function update_sites!( sc::site_collection, g::Int64, N::Int64 )
  if !sc.in_use 
    return
  end
  for index in sc.active  # updates sites to the next generation
    # update_selected does the equivalent of proportional selection
    new_allele_freq = update_selected( index, N, sc.list[index].previous_allele_freq, sc.list[index].fitness_coefficient )
    innov = sc.list[index]
    if new_allele_freq > 0 && new_allele_freq < N
      global sum_counts += new_allele_freq
      global sum_sq_counts += new_allele_freq^2
      innov.sum_gens += 1
      innov.sum_counts += new_allele_freq
      innov.sum_sq_counts += new_allele_freq^2
      heteroz =  1.0 - watterson_homozygosity([N-new_allele_freq, new_allele_freq])
      global sum_heteroz += heteroz
      global sum_sq_heteroz += heteroz^2
      global sum_gens +=1
      innov.sum_heteroz += heteroz
      innov.sum_sq_heteroz += heteroz^2
    end
    innov.previous_allele_freq = new_allele_freq  # save for the next call to update_selected.
    if new_allele_freq == 0 ||  new_allele_freq >= N 
      sc.sum_gens += innov.sum_gens
      sc.sum_counts += innov.sum_counts
      sc.sum_sq_counts += innov.sum_sq_counts
      sc.sum_heteroz += innov.sum_heteroz
      sc.sum_sq_heteroz += innov.sum_sq_heteroz
      sc.sum_generations += innov.sum_gens
      Base.pop!( sc.active,index)
      innov.final_gen = g
      if new_allele_freq == 0  # extinction
        Base.push!( sc.extinct,index)
      elseif new_allele_freq >= N # fixation
        Base.push!( sc.fixed,index)
      end
    end
  end
end

@doc """ function update_neutral() 
  Updates site according to the Wright-Fisher model of drift.  
  Since mutation has already happened by the time a new site is established, there is no mutation.
""" 
function update_neutral( site::Int64, N::Int64, old_allele_freq::Int64 )
  p = Float64(old_allele_freq)/N
  new_allele_freq = rand(Binomial(N,p))
  return new_allele_freq
end 

@doc """ function update_selected() 
  Updates site according to the Wright-Fisher model of drift with selection.  
  Since mutation has already happened by the time a new site is established, there is no mutation.
  This function is called by update_sites!() in site_collection.jl.

""" 
function update_selected( site::Int64, N::Int64, old_allele_freq::Int64, select_coef::Float64 )
  p = min(1.0,(select_coef*Float64(old_allele_freq)/N))
  new_allele_freq = rand(Binomial(N,p))
  return new_allele_freq
end  


function N_inf_sites( sc::site_collection )
  if !sc.in_use 
    return
  end
  sum_N = 0
  for index in sc.active  
    sum_N += sc.list[index].history[end]
  end
  sum_N
end

@doc """ check_globals_against_site_collection()
  Check that the values of global variables approximately agree with the corresponding variables in sc.
"""
function check_globals_against_site_collection( sc::site_collection )
  global sum_counts 
  global sum_sq_counts 
  global sum_heteroz 
  global sum_sq_heteroz 
  global sum_gens 
  @assert isapprox(sc.sum_gens,sum_gens)
  @assert isapprox(sc.sum_counts,sum_counts)
  @assert isapprox(sc.sum_sq_counts,sum_sq_counts)
  @assert isapprox(sc.sum_heteroz,sum_heteroz)
  @assert isapprox(sc.sum_sq_heteroz,sum_sq_heteroz)
end

function print_summary( sc::site_collection; print_lists::Bool=false )
  println("print summary of site collection: ")
  global sum_counts 
  global sum_sq_counts 
  global sum_heteroz 
  global sum_sq_heteroz 
  global sum_gens 

  #global sum_counts
  #println("psum sum_counts: ",sum_counts)
  if !sc.in_use 
    return
  end
  if print_lists
    println("active list: ",sc.active)
    println("fixed list: ",sc.fixed)
    println("extinct list: ",sc.extinct)
  end
  println("number active: ",length(sc.active))
  println("number fixed: ",length(sc.fixed))
  println("number extinct: ",length(sc.extinct))
  println("fixed fraction: ",fixed_fraction(sc))
  println("sum_counts: ",sc.sum_counts,"  glb: ",sum_counts)
  println("sum_sq_counts: ",sc.sum_sq_counts,"  glb: ",sum_sq_counts)
  println("sum_heteroz: ",sc.sum_heteroz,"  glb: ",sum_heteroz)
  println("sum_sq_heteroz: ",sc.sum_sq_heteroz,"  glb: ",sum_sq_heteroz)
  println("sum_gens: ",sc.sum_generations,"  glb: ",sum_gens)
  (avg,stderr) = sites_per_gen(sc) 
  println("avg sites per generation: ",avg)
  println("stderr sites per generation: ",stderr)
  (avg,stderr) = heterozygosity_per_gen(sc) 
  println("avg per generation heterozygosity: ",avg)
  println("stderr per generation heterozygosity: ",stderr)
  println("avg time to fixation: ",average_time_to_fixation(sc))
  println("avg time to extinction: ",average_time_to_extinction(sc))
  println("avg fitness fixed: ",average_fitness_fixed(sc))
  println("avg fitness extinct: ",average_fitness_extinct(sc))
  println("avg fitness all: ",average_fitness_all(sc))
end

function average_time_to_extinction( innov_collection::site_collection )
  if !innov_collection.in_use 
    return 0.0
  end
  if length(innov_collection.extinct) == 0
    println("no extinct sites")
    return -1.0
  end
  #print("average_time_to_extinction:  ")
  sum = 0
  for i in innov_collection.extinct 
    sum += innov_collection.list[i].final_gen - innov_collection.list[i].start_gen
    #println("id: ",innov_collection.list[i].identifier,"  time: ",(innov_collection.list[i].final_gen-innov_collection.list[i].start_gen),"  sum: ",sum)
  end
  return Float64(sum)/length(innov_collection.extinct)
end 

function average_time_to_fixation( innov_collection::site_collection )
  if !innov_collection.in_use 
    return 0.0
  end
  if length(innov_collection.fixed) == 0
    println("no fixed sites")
    return -1.0
  end
  sum = 0
  for i in innov_collection.fixed 
    sum += innov_collection.list[i].final_gen - innov_collection.list[i].start_gen
  end
  return Float64(sum)/length(innov_collection.fixed)
end 

# Fraction of fixed sites out of fixed plus extinct sites
function fixed_fraction( innov_collection::site_collection )
  if !innov_collection.in_use 
    return 0.0
  end
  if  length(innov_collection.extinct) == 0 && length(innov_collection.fixed) == 0
    println("no extinct and no fixed sites")
    return -1.0
  end
  return length(innov_collection.fixed)/Float64(length(innov_collection.fixed)+length(innov_collection.extinct))
end

function average_fitness_fixed( sc::site_collection )
  if !sc.in_use 
    return 0.0
  end
  if  length(sc.fixed) == 0
    return -1.0
  end
  sum = 0.0
  for i in sc.fixed 
    sum += sc.list[i].fitness_coefficient
  end
  return sum/length(sc.fixed)
end

function average_fitness_extinct( sc::site_collection )
  if !sc.in_use 
    return 0.0
  end
  if  length(sc.extinct) == 0
    return -1.0
  end
  sum = 0.0
  for i in sc.extinct 
    sum += sc.list[i].fitness_coefficient
  end
  return sum/length(sc.extinct)
end

function average_fitness_all( sc::site_collection )
  if !sc.in_use 
    return 0.0
  end
  if  length(sc.extinct) == 0
    return -1.0
  end
  sum = 0.0
  for i in sc.extinct 
    sum += sc.list[i].fitness_coefficient
  end
  for i in sc.fixed 
    sum += sc.list[i].fitness_coefficient
  end
  return sum/(length(sc.extinct) + length(sc.fixed))
end

@doc """ function count_adv_del_fixed( )
  Counts the number of fixed sites that are advantageous, and the number that are deleterious at the end of the run.
"""
function count_adv_del_fixed( sc::site_collection )
  #println("count_adv_del_fixed")
  if !sc.in_use 
    return 0.0
  end
  if  length(sc.fixed) == 0
    return 0, 0
  end
  count_adv = 0
  count_del = 0
  for i in sc.fixed 
    fit_coef = sc.list[i].fitness_coefficient 
    #println("i: ",i,"  fit_coef: ",fit_coef)
    if fit_coef > 1.0
      count_adv += 1
    end
    if fit_coef < 1.0
      count_del += 1
    end
  end
  sc.count_fixed_adv = count_adv
  sc.count_fixed_del = count_del
  return count_adv, count_del
end

function sites_per_gen( sc::site_collection )
  global sum_gens 
  global sum_counts
  global sum_sq_counts
  global counts_var = sqrt(sum_sq_counts/(sum_gens-1) - sum_counts^2/sum_gens/(sum_gens-1))/sum_gens
  sumgens = sc.sum_generations
  average = sc.sum_counts/sumgens
  coef_var = sqrt(sc.sum_sq_counts/(sumgens-1) - sc.sum_counts^2/sumgens/(sumgens-1))/sumgens
  @assert isapprox( coef_var, counts_var )
  if coef_var >= 0
    stderr = coef_var
  else
    stderr = -1.0
  end
  return (average, stderr)
end

function heterozygosity_per_gen( sc::site_collection )
  global sum_gens 
  global sum_heteroz
  global sum_sq_heteroz
  heteroz_var = sqrt(sum_sq_heteroz/(sum_gens-1) - sum_heteroz^2/sum_gens/(sum_gens-1))/sum_gens
  sumgens = sc.sum_generations
  average = sc.sum_heteroz/sumgens
  coef_var = sqrt(sc.sum_sq_heteroz/(sumgens-1) - sc.sum_heteroz^2/sumgens/(sumgens-1))/sumgens
  @assert isapprox( coef_var, heteroz_var )
  if coef_var >= 0
    stderr = coef_var
  else
    stderr = -1.0
  end
  return (average, stderr)
end
