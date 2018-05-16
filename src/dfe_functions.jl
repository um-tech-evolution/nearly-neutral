#=
Definitions of the dfe (distribution of fitness effects) functions which define
   selection coefficients in terms of Gamma distributions.
Split off from nn_poplist.jl  on 5/14/18.
=#
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
