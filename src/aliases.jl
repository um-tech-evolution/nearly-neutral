export Population 

using Distributions
using Random
using QuadGK
import Random.seed!

const DIST_TYPE = Dict{Int64,Float64}
const JOINT_DIST_TYPE = Dict{Tuple{Int64,Int64},Float64}
const Population = Array{Int64,1}
const PopList = Array{Population,1}

