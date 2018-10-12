include("../src/util.jl")
include("../src/LUMND.jl")

using .LUMND

d = 2
A = Ad(5,d) # n^d laplacian
N = size(A,1)
leaf = 5
lvl = 3 #Int64(ceil(log(N/leaf)/log(2)))
b = rand(N)
(subseps, hrch, dofs) = LUMND.mnd(A, lvl)
# @time tree = LUMND.factorize(A, lvl)
# @time x = LUMND.solve(tree, b)
# @time xref = A\b
# @show norm(A*x-b)
