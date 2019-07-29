# compute the entropy of a probability distribution given as a function p
# defined as a hash table over integers
# IMPORTANT:  One must include "use_populist=true" in the configuration file.

#include("types.jl")
export pop_to_dist, entropy, relative_entropy

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

function entropy( p::DIST_TYPE ) 
  result = 0.0
  for x in keys(p)
    result += - p[x] * log2( p[x] )
    #println(x, "  result: ",result)
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

