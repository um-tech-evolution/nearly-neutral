export innovation_type, innovation,  
    innovation_collection, ic_push!, update_innovations!, update_neutral, update_selected,
    print_summary, average_time_to_extinction, average_time_to_fixation, 
    fixed_fraction, average_fitness_fixed, average_fitness_extinct, average_fitness_all,
    innovations_per_gen, heterozygosity_per_gen


using Base.Test
#=
Stores the properties of an "innovation" (possibly deleterious or advantageous).
=#

type innovation_type
  identifier::Int64   # integer for this innovation, not sure if this is needed
  start_gen::Int64    # the generation (time step) when the innovation was generated
  final_gen::Int64    # the generation when the innovation went extinct.  Should be zero while innovation is evolving
  fitness_coefficient::Float64    # The selection coefficient of the innovation (allele)
  sum_counts::Int64   # the current accumulated sum of the counts per generation for infinite sites
  sum_heteroz::Float64  # the current accumulated sum of heterozygosities for infinite sites
  previous_allele_freq::Int64  # Allele frequence of the previous generation
end

@doc """ function innovation() 
  Initializes the properties of the innovation.
"""
function innovation( id::Int64, N::Int64, start_gen::Int64, fitness_coef::Float64=1.0 )
  #println("innovation: id: ",id,"  N: ",N,"  start_gen: ",start_gen)
  initial_heteroz = 1.0 - watterson_homozygosity([N-1,1])
  return innovation_type( id, start_gen, 0, fitness_coef, 1, initial_heteroz, 1  )
end

#=
Stores a collection of innovations.
Innovations are partitioned into 3 subsets, active, fixed, and extinct.
=#

type innovation_collection
  N::Int64                # popsize
  list::Dict{Int64,innovation_type} #    Stores a collection of innovations
  active::IntSet   # indices of actively evolving innovations
  fixed::IntSet   # indices of innovations that have fixed
  extinct::IntSet   # indices of innovations that have gone extinct
  fix_minimum::Float64  # minimum fraction of popsize for fixation for infinite alleles
  sum_counts::Int64   # sum of counts of mutant alleles for infinite sites
  sum_sq_counts::Float64   # sum of counts of mutant alleles for infinite sites
  sum_heteroz::Float64  # sum of heterozygosities for infinite sites
  sum_sq_heteroz::Float64  # sum of heterozygosities for infinite sites
  sum_generations::Int64  # Total number of generations corresponding to sum_counts and sum_heteroz
  sum_sq_generations::Float64  # Total number of generations corresponding to sum_counts and sum_heteroz
  count_fixed_del::Int64  # Number of fixed mutations that are deleterious (fit_coef < 1.0)
  count_fixed_adv::Int64  # Number of fixed mutations that are advantageous (fit_coef > 1.0)
  in_use::Bool      # If false, not used
end

# Constructor for a new empty innovation collection
function innovation_collection( N::Int64, in_use::Bool=true )  
  #innovation_collection( N, Dict{Int64,innovation_type}(), IntSet(), IntSet(), IntSet(), 1.0, 0, 0.0, 0,0,0, in_use )
  #innovation_collection( N, Dict{Int64,innovation_type}(), IntSet(), IntSet(), IntSet(), 1.0, 0,0, in_use )
  innovation_collection( N, Dict{Int64,innovation_type}(), IntSet(), IntSet(), IntSet(), 1.0, 0, 0, 0.0, 0.0, 0, 0, 0, 0, in_use )
end

# Constructor for a new empty innovation collection with a value for fix_minimum
function innovation_collection( N::Int64, fix_min::Float64, in_use::Bool=true )  
  #innovation_collection( N, Dict{Int64,innovation_type}(), IntSet(), IntSet(), IntSet(), fix_min, 0, 0.0, 0,0,0, in_use )
  #innovation_collection( N, Dict{Int64,innovation_type}(), IntSet(), IntSet(), IntSet(), fix_min, 0,0, in_use )
  innovation_collection( N, Dict{Int64,innovation_type}(), IntSet(), IntSet(), IntSet(), fix_min, 0, 0, 0.0, 0.0, 0,0,0,0, in_use )
end

@doc """ ic_push!()
 Adds an innovation to the collection.
"""
function ic_push!( innov_collection::innovation_collection, innov::innovation_type )
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
  #println("new innovation: id: ",innov.identifier)
  innov_collection
end

@doc """ fix_test()
The most basic test for fixation, namely that the allele frequency is the population size.
This is used in the infitite sites model.
"""
function fix_test( N::Int64, new_allele_freq::Int64 )
  return N == new_allele_freq
end

@doc """ function update_innovations!( )
 Update all active innovations, and make some extinct and fixed
 This version is used only by infsites.jl.
"""
function update_innovations!( ic::innovation_collection, g::Int64, N::Int64 )
  if !ic.in_use 
    return
  end
  #println("active: ",ic.active,"  fixed: ",ic.fixed,"  extinct: ",ic.extinct)
  for index in ic.active  # updates sites to the next generation
    #println("update innovations: id: ",ic.list[index].identifier,"  start gen: ",ic.list[index].start_gen) 
    # update_selected does the equivalent of proportional selection
    new_allele_freq = update_selected( index, N, ic.list[index].previous_allele_freq, ic.list[index].fitness_coefficient )
    #println("id: ",ic.list[index].identifier,"  prev_allele_freq: ",ic.list[index].previous_allele_freq,"  new_allele_freq: ",new_allele_freq)
    innov = ic.list[index]
    if new_allele_freq > 0 && new_allele_freq < N
      innov.sum_counts += new_allele_freq
      innov.sum_heteroz += 1.0 - watterson_homozygosity([N-new_allele_freq, new_allele_freq])
    end
    #println("new_allele_freq: ",new_allele_freq,"  sum_counts: ",innov.sum_counts,"  sum_heteroz: ",innov.sum_heteroz)
    innov.previous_allele_freq = new_allele_freq  # save for the next call to update_selected.
    if new_allele_freq == 0 ||  new_allele_freq >= N 
      ic.sum_counts += innov.sum_counts
      ic.sum_sq_counts += innov.sum_counts^2
      ic.sum_heteroz += innov.sum_heteroz
      ic.sum_sq_heteroz += innov.sum_heteroz^2
      ic.sum_generations += g - innov.start_gen
      ic.sum_sq_generations += (g - innov.start_gen)^2
      Base.pop!( ic.active,index)
      innov.final_gen = g
      if new_allele_freq == 0  # extinction
        Base.push!( ic.extinct,index)
        #println("extinct index:",index,"  gen:",g,"  startg:",innov.start_gen)
      elseif new_allele_freq >= N # fixation
        Base.push!( ic.fixed,index)
        #println("fixed index:",index,"  gen:",g,"  startg:",innov.start_gen)
      end
    end
  end
end

@doc """ function update_innovations!( )
 Update all active innovations, and make some extinct and fixed
 This version is called by the infinite alleles model (either neutral_poplist.jl or nearly_neutral_poplist.jl).
 popcounter is a dictionary that maps alleles to their frequencies in the current generation population
"""
function update_innovations!( ic::innovation_collection, g::Int64, N::Int64, popcounter::Dict{Int64,Int64}, fixation_test::Function=fix_test )
  if !ic.in_use 
    return
  end
  for index in ic.active  # updates sites to the next generation
    new_allele_freq = get( popcounter, index, 0 )
    #ic_update!(ic,index,g,new_allele_freq)
    if new_allele_freq == 0  # extinction
      Base.pop!( ic.active,index)
      Base.push!( ic.extinct,index)
      ic.list[index].final_gen = g
    elseif new_allele_freq >= Int(ceil(ic.fix_minimum*N) )  # fixation
      Base.pop!( ic.active,index)
      Base.push!( ic.fixed,index)
      ic.list[index].final_gen = g
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
  #println("update_neutral: ",site,"  N: ",N,"  new_allele_freq: ",new_allele_freq)
  return new_allele_freq
end 

@doc """ function update_selected() 
  Updates site according to the Wright-Fisher model of drift with selection.  
  Since mutation has already happened by the time a new site is established, there is no mutation.
  This function is called by update_innovations!() in innovation_collection.jl.

""" 
function update_selected( site::Int64, N::Int64, old_allele_freq::Int64, select_coef::Float64 )
  p = min(1.0,(select_coef*Float64(old_allele_freq)/N))
  new_allele_freq = rand(Binomial(N,p))
  #println("update selected: site: ",site,"  old allele_freq: ",old_allele_freq,"  new allele_freq: ",new_allele_freq,"  sel coef: ",select_coef,"  p: ",p)
  return new_allele_freq
end  
function N_inf_sites( ic::innovation_collection )
  if !ic.in_use 
    return
  end
  sum_N = 0
  for index in ic.active  
    sum_N += ic.list[index].history[end]
  end
  sum_N
end

# TODO:  Move to conformist_poplist.jl or delete
function compute_turnovers( pop1::Population, pop2::Population, N_mu::Float64, Ylist::Vector{Int64},
    Zsum_list::Vector{Int64}, count_list::Vector{Int64} )
  i = 1
  for y in Ylist
    if Float64(y) < 5.0*N_mu
      Zsum_list[i] += turnover( pop1, pop2, y )
      count_list[i] += 1
    end
    i+= 1
  end
end

function print_summary( ic::innovation_collection; print_lists::Bool=false )
  if !ic.in_use 
    return
  end
  if print_lists
    println("active list: ",ic.active)
    println("fixed list: ",ic.fixed)
    println("extinct list: ",ic.extinct)
  end
  println("number active: ",length(ic.active))
  println("number fixed: ",length(ic.fixed))
  println("number extinct: ",length(ic.extinct))
  println("fixed fraction: ",fixed_fraction(ic))
  println("sum_counts: ",ic.sum_counts)
  println("sum_heteroz: ",ic.sum_heteroz)
  println("sum_gens: ",ic.sum_generations)
  (avg,stderr) = innovations_per_gen(ic) 
  println("avg innovations per generation: ",avg)
  println("stderr innovations per generation: ",stderr)
  (avg,stderr) = heterozygosity_per_gen(ic) 
  println("avg per generation heterozygosity: ",avg)
  println("stderr per generation heterozygosity: ",stderr)
  println("avg time to fixation: ",average_time_to_fixation(ic))
  println("avg time to extinction: ",average_time_to_extinction(ic))
  println("avg fitness fixed: ",average_fitness_fixed(ic))
  println("avg fitness extinct: ",average_fitness_extinct(ic))
  println("avg fitness all: ",average_fitness_all(ic))
end

function average_time_to_extinction( innov_collection::innovation_collection )
  if !innov_collection.in_use 
    return 0.0
  end
  if length(innov_collection.extinct) == 0
    println("no extinct innovations")
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

function average_time_to_fixation( innov_collection::innovation_collection )
  if !innov_collection.in_use 
    return 0.0
  end
  if length(innov_collection.fixed) == 0
    println("no fixed innovations")
    return -1.0
  end
  sum = 0
  for i in innov_collection.fixed 
    sum += innov_collection.list[i].final_gen - innov_collection.list[i].start_gen
  end
  return Float64(sum)/length(innov_collection.fixed)
end 

# Fraction of fixed innovations out of fixed plus extinct innovations
function fixed_fraction( innov_collection::innovation_collection )
  if !innov_collection.in_use 
    return 0.0
  end
  if  length(innov_collection.extinct) == 0 && length(innov_collection.fixed) == 0
    println("no extinct and no fixed innovations")
    return -1.0
  end
  return length(innov_collection.fixed)/Float64(length(innov_collection.fixed)+length(innov_collection.extinct))
end

function average_fitness_fixed( ic::innovation_collection )
  if !ic.in_use 
    return 0.0
  end
  if  length(ic.fixed) == 0
    return -1.0
  end
  sum = 0.0
  for i in ic.fixed 
    sum += ic.list[i].fitness_coefficient
  end
  return sum/length(ic.fixed)
end

function average_fitness_extinct( ic::innovation_collection )
  if !ic.in_use 
    return 0.0
  end
  if  length(ic.extinct) == 0
    return -1.0
  end
  sum = 0.0
  for i in ic.extinct 
    sum += ic.list[i].fitness_coefficient
  end
  return sum/length(ic.extinct)
end

function average_fitness_all( ic::innovation_collection )
  if !ic.in_use 
    return 0.0
  end
  if  length(ic.extinct) == 0
    return -1.0
  end
  sum = 0.0
  for i in ic.extinct 
    sum += ic.list[i].fitness_coefficient
  end
  for i in ic.fixed 
    sum += ic.list[i].fitness_coefficient
  end
  return sum/(length(ic.extinct) + length(ic.fixed))
end

@doc """ function count_adv_del_fixed( )
  Counts the number of fixed innovations that are advantageous, and the number that are deleterious at the end of the run.
"""
function count_adv_del_fixed( ic::innovation_collection )
  #println("count_adv_del_fixed")
  if !ic.in_use 
    return 0.0
  end
  if  length(ic.fixed) == 0
    return 0, 0
  end
  count_adv = 0
  count_del = 0
  for i in ic.fixed 
    fit_coef = ic.list[i].fitness_coefficient 
    #println("i: ",i,"  fit_coef: ",fit_coef)
    if fit_coef > 1.0
      count_adv += 1
    end
    if fit_coef < 1.0
      count_del += 1
    end
  end
  return count_adv, count_del
end

function innovations_per_gen( ic::innovation_collection )
  sum_gens = ic.sum_generations
  println("avg innov sum_gens: ",sum_gens)
  println("avg innov sum_sq_counts: ",ic.sum_sq_counts)
  average = ic.sum_counts/sum_gens
  var = ((ic.sum_sq_counts/(sum_gens-1) - ic.sum_counts^2/sum_gens/(sum_gens-1))/sum_gens)
  println("avg innov var: ",var)
  stderr = sqrt((ic.sum_sq_counts/(sum_gens-1) - ic.sum_counts^2/sum_gens/(sum_gens-1))/sum_gens)
  return (average, stderr)
end

function heterozygosity_per_gen( ic::innovation_collection )
  sum_gens = ic.sum_generations
  average = ic.sum_heteroz/sum_gens
  var = ((ic.sum_sq_heteroz/(sum_gens-1) - ic.sum_heteroz^2/sum_gens/(sum_gens-1))/sum_gens)
  println("avg hetero var: ",var)
  stderr = sqrt((ic.sum_sq_heteroz/(sum_gens-1) - ic.sum_heteroz^2/sum_gens/(sum_gens-1))/sum_gens)
  return (average, stderr)
end
