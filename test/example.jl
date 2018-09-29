include("../src/util.jl")
include("../src/LUMND.jl")

using .LUMND

A = Ad(64,3) # n^3 laplacian
N = size(A,1)
leaf = 64
lvl = Int64(ceil(log(N/leaf)/log(2)))
b = rand(N)
@time tree = LUMND.factorize(A, lvl)
@time x = LUMND.solve(tree, b)
@time xref = A\b
@show norm(A*x-b)
