# compute the entropy of a probability distribution given as a function p
# defined as a hash table over integers
# IMPORTANT:  One must include "use_populist=true" in the configuration file.

export pop_to_dist, entropy, relative_entropy, mutual_information, sherwin_mi

p = DIST_TYPE( 1 => 0.5, 2=> 0.25, 3=> 0.25 )
p1= DIST_TYPE( 1 => 0.8, 2=> 0.1, 3=> 0.1 )
q0 = DIST_TYPE( 1 => 0.5, 2=> 0.50 )
q1 = DIST_TYPE( 1 => 0.75, 2=> 0.25 )
q2 = DIST_TYPE( 1 => 0.25, 2=> 0.25, 3=>0.5 )
z1 = DIST_TYPE( 1 => 1.0 )
z2 = DIST_TYPE( 2 => 1.0 )
z3 = DIST_TYPE( 3 => 1.0 )
z4 = DIST_TYPE( 3 => 0.9 )

function dist_check( p::DIST_TYPE )
  sum = 0.0
  for x in keys(p)
    if p[x] < 0.0
      error(" value of distribution ", p, " is negative ")
    end
    sum += p[x]
  end
  if !isapprox(sum,1.0)
      error(" sum of distribution ", p, " is not 1.0")
  end
end

# pop can include repeated elements.
# If i is repeated n times in pop, result[i] is Float64(n)/N.
function pop_to_dist( pop::Population )
  N = length(pop)
  result = DIST_TYPE()
  for i in pop
    #pval = get(result,i,0.0)
    result[i] = get(result,i,0.0) + 1.0/N
    #println("i: ",i,"  pval: ",pval,"  result: ",result)
  end
  dist_check(result)
  result
end

# Convert a list of populations to a probability distribution table with one row per population
# There is one column per "allele" of the combined population
# The sum of the entries in the result is 1.0.
# The sum of the entries in row i is length(P[i])/(sum(length(p) for p in P) 
function pops_to_tbl( P::Vector{Population} )
  m = length(P)   # the number of rows in result
  combined_pop = reduce(vcat, P[i] for i = 1:m )
  N = length(combined_pop)
  allele_list = unique(combined_pop)
  n = length(allele_list)   # The number columns in result
  result = zeros(Float64,m,n)
  for i = 1:length(P)
    dist = Dict( allele_list .=> zeros(Float64,n))
    for p in P[i]
      dist[p] += 1.0/N
    end
    result[i,:] = [ dist[a] for a in allele_list ]
  end
  result
end

function entropy( p::DIST_TYPE ) 
  result = 0.0
  for x in keys(p)
    result += - p[x] * log2( p[x] )
    #println("x: ",x, "result_inc: ", p[x] * log2( p[x]),"  result: ",result)
  end
  result
end

function entropy( p::Population )
  entropy( pop_to_dist(p) )
end

# 
function relative_entropy( q::DIST_TYPE, p::DIST_TYPE )
  result = 0.0
  for x in keys(q)  # Assume that any keys of p that aren't in q have value 0 for q
    qval = get(q,x,0.0)
    #result += q[x] * log2( q[x]/p[x] )
    if qval > 0.0
      result += qval * log2( qval/p[x] )
    end
    #println("x: ",x,"  qval: ",qval,"  result: ",result)
  end
  result
end

function relative_entropy( q::Population, p::Population )
  relative_entropy( pop_to_dist(q), pop_to_dist(p) )
end

function relative_entropy( q::Population, pdist::DIST_TYPE)
  relative_entropy( pop_to_dist(q), pdist )
end

# Entropy of a population of allele counts
function counts_entropy( pcounts::Population )
  N = sum(pcounts)
  -sum( map(x->x*log2(x), pc/N ))
end


# Mutual information of a collection of subpopulations by the method given in Sherwin, 
#     Chao, Jost, Smouse  (2017) supplementary info.
# Their formula is  I = (^1 H_\gamma = \overline{^1 H_\alpha ), but my methods gives 1 - I.
function mutual_information( P::Vector{Population} )
  sherwin_mutual_information( pops_to_tbl( P ) )
end

# Mutual information by the "standard" defintion of eq. 2.28 of Cover and Thomas.
#   This is the mutual information of a joint distribution.
function mutual_information( tbl::Array{Float64, 2} )
  (m,n) = size(tbl)
  mr1 = [ sum( tbl[i,:]) for i = 1:m ]
  #println("mr1: ",mr1)
  mr2 = [ sum( tbl[:,j]) for j = 1:n ]
  #println("mr2: ",mr2)
  mi_tbl = zeros(Float64,size(tbl))
  for i = 1:m
    for j = 1:n
      mi_tbl[i,j] = tbl[i,j] != 0.0 ? tbl[i,j]*log2(tbl[i,j]/mr1[i]/mr2[j]) : 0.0
      #println("(i,j): ",(i,j),"  tbl[i,j]: ",tbl[i,j],"  mr1[i]: ",mr1[i]," mr2[j]: ",mr2[j],"  etbl[i,j]: ",mi_tbl[i,j])
    end
  end
  #println("mi_tbl: ",mi_tbl)
  sum(mi_tbl)
end

entf(x) = x != 0 ? -x*log2(x) : 0.0

# tbl is the probability distribution for a total population which is made up of m subpopulations
#    all over the same allele set of size n.
# Row i of tbl corresponds to a probability distribution over subpopulation i.
function sherwin_mutual_information( tbl::Array{Float64,2} )
  (m,n) = size(tbl)
  # weight[i] is the fraction of the total population corresponding to subpopulation i.
  weights = [ sum(tbl[i,:]) for i = 1:m ]  
  s = [ sum( map(entf,tbl[i,:]/weights[i])) for i = 1:m ]
  ent_total = sum(map(entf,[sum(tbl[:,j]) for j = 1:n]))
  ent = ent_total - sum(s[i]*weights[i] for i = 1:m) 
  #println(" -- sh ent: ",ent)
  ent
end
  
