#=
Estimates of theta = N_e mu based on Watterson "The homozygosity test of neutrality" (1978).

allele_freqs is the vector of allele frequecies of the infinite alleles model population.
=#

export watterson_homozygosity, p_homozygosity, watterson_theta, IQV, w_homoz 

function watterson_homozygosity( allele_freqs::Population )
  if length(allele_freqs) == 0
    return 0.0   # Not sure that this is valid
  end
  n = sum(allele_freqs)
  sum_a_sq = 0
  for a in allele_freqs
    sum_a_sq += a^2
  end
  Float64(sum_a_sq)/Float64(n^2)
end
#=
function watterson_theta( allele_freqs::Config )
  1.0/watterson_homozygosity( allele_freqs ) - 1.0
end
=#

function watterson_theta( allele_freqs::Population )
  1.0/watterson_homozygosity( allele_freqs ) - 1.0
end

@doc """ function IQV()
Allan Wilcox's "index of quantitative variation" defined by formula (7) 
  of Mark E. Madsen's "Neutral Cultural Transmission in Time Averaged
Archaeological Assemblages" arXiv:1024.2043v2 (2012). 
A measure of "evenness".  IQV of a population with a single trait is 0,
while IQV of a population with an even distribution of traits is close to 1.0.
"""
function IQV( allele_freqs::Population )
  k = length(allele_freqs)
  if k == 1
    return 0.0
  end
  Fk = Float64(k)
  wh = watterson_homozygosity( allele_freqs )
  (1.0-wh)*Fk/(Fk-1.0)
end
  
