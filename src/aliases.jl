export Population 

using Distributions
using Random
using QuadGK
import Random.seed!

if !@isdefined(MyInt)
  const MyInt = UInt8 # When called by a fuction in Main.CGP, may have to agree with MyInt defined in CGP.jl
end
const DIST_TYPE = Dict{Any,Float64}
const IDIST_TYPE = Dict{Int64,Float64}
const SDIST_TYPE = Dict{String,Float64}
const DIST2_TYPE = Dict{Tuple{Any,Any},Float64}
const IDIST2_TYPE = Dict{Tuple{Int64,Int64},Float64}
const SDIST2_TYPE = Dict{Tuple{String,String},Float64}
const Population = Vector
const IPopulation = Array{Int,1}
const MIPopulation = Array{MyInt,1}
const SPopulation = Array{String,1}
const FPopulation = Array{Float64,1}
const PopList = Array{Population,1}

