#! /usr/local/julia/bin/julia
using DataFrames
dfe_neutral(x) = 1.0
if length(ARGS) == 0
  simname = "../experiments/examples/in_example1"
else
  simname = ARGS[1]
end
include("$(simname).jl")
println("simname: ",simname)

#simname = "examples/in_example2"
df = readtable("$(simname).csv", makefactors=true, allowcomments=true)
println("df: ",df)
