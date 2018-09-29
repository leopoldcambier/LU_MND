include("../src/util.jl")
include("../src/LUMND.jl")

using .LUMND

for i = 1:1
    d = rand(2:3)
    n = d == 2 ? rand(5:30) : rand(5:10)
    leaf = rand([1,5,10,15])
    lvl = Int64(round(log(n^d / (leaf)) / log(2)))

    d = 2
    n = 4
    lvl = 1
    A = Ad(n,d)
    b = rand(size(A, 1))
    tree = LUMND.factorize(A, lvl)
    x = LUMND.solve(tree, b)
    @printf("%2d | %5d -> %e\n", i, size(A,1), norm(A*x-b))
end
