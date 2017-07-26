#! /usr/local/julia/bin/julia
using DataFrames
if length(ARGS) == 0
  simname = "../experiments/examples/fix_example1"
else
  simname = ARGS[1]
end
include("$(simname).jl")
println("simname: ",simname)

#simname = "examples/fi_example2"
df = readtable("$(simname).csv", makefactors=true, allowcomments=true)
println("df: ",df)
