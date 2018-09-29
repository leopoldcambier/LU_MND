include("../src/util.jl")
include("../src/LUMND.jl")

using .LUMND
using Test
using Random

Random.seed!(0)

@testset "End-to-end accuracy" begin
    @printf("  # |     N lvl         k ->   |Ax-b|   |x-x*|\n")
    for i = 1:10
        d = rand(2:3)
        n = d == 2 ? rand(5:30) : rand(5:10)
        leaf = rand([1,2,5,10,15,20])
        lvl = Int64(round(log(n^d / (leaf)) / log(2)))
        A = Ad(n,d)
        b = rand(size(A, 1))
        tree = LUMND.factorize(A, lvl)
        x = LUMND.solve(tree, b)
        xref = A\b
        k = cond(A, 1)
        @printf("%3d | %5d  %2d  %2.2e -> %4.2e %4.2e\n", i, size(A,1), lvl, k, norm(A*x-b)/norm(b), norm(xref - x)/norm(xref))
        @test(norm(A*x-b)/norm(b)     <     1e-13)
        @test(norm(x-xref)/norm(xref) < k * 1e-13)
    end
end
