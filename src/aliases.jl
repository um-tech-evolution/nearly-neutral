export Population 

using Distributions
using Random
import Random.seed!

const Population = Array{Int64,1}
const PopList = Array{Population,1}
const DIST_TYPE = Dict{Int64,Float64}

