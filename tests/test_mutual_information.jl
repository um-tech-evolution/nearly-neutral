include("../src/aliases.jl")
include("../src/entropy.jl")

P1 = Any[1,2,3,4]
P2 = Any[5,6,6,7]
P3 = Any[3,4,5,5]

m11 = mutual_information([P1,P1])
m12 = mutual_information([P1,P2])
m13 = mutual_information([P1,P3])
m23 = mutual_information([P2,P3])
m123 = mutual_information([P1,P2,P3])

println("m11: ",m11,"  m12: ",m12,"  m11: ",m11,"  m23: ",m23,"  m123: ",m123)

Pa = Any[1,1,2]   # For Dewar paper, p1 = 2/3
Pb = Any[1,2,2]   # For Dewar paper, p2 = 1/3
mPaPb =  mutual_information([P1,P2])
println("mPaPb: ",mPaPb)

function build_tbl( pvec::Vector{Float64} )
  m = length(pvec)
  result = zeros(Float64,(m,2))
  for i = 1:m
    result[i,:] = [pvec[i], 1-pvec[i]]
  end
  result/m
end

function build_tbl( parray::Array{Float64,2} )
  (m,n) = size(parray)
  result = zeros(Float64,(m,n+1))
  result[:,1:n] = parray
  diff_vec = [ 1.0-sum(parray[i,:]) for i = 1:m]
  @assert reduce( min, diff_vec ) >= 0.0
  result[:,n+1] = diff_vec
  result/m
end

# Generate a random m by n table whose rows sum to 1/m
function rand_tbl( m::Int64, n::Int64 )
  rtbl = rand(m,n)
  for i = 1:m
    rtbl[i,:] = rtbl[i,:]/sum(rtbl[i,:])/m
  end
  rtbl
end

pvec = [2/3,1/4,2/5,1/2]
ttbl = build_tbl(pvec)
std_mi = mutual_information(ttbl)
sh_mi = sherwin_mutual_information(ttbl)
@assert isapprox( std_mi, sh_mi )

parr = [[0.5 0.2]; [0.1 0.4]; [1/4 1/2]; [1/8 5/8]]
ttbl = build_tbl(parr)
std_mi = mutual_information(ttbl)
sh_mi = sherwin_mutual_information(ttbl)
@assert isapprox( std_mi, sh_mi )

m = 3
n = 5
rtbl = rand_tbl(m,n)
std_mi = mutual_information(rtbl)
sh_mi = sherwin_mutual_information(rtbl)
@assert isapprox( std_mi, sh_mi )


