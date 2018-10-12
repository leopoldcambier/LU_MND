include("../src/util.jl")
include("../src/LUMND.jl")

using .LUMND
using Test
using Random
using SHA
using MatrixDepot

Random.seed!(0)

@testset "Partition" begin
    for i = 1:50
        d = rand(2:3)
        n = d == 2 ? rand(5:20) : rand(5:10)
        A = Ad(n,d)
        N = size(A,1)
        maxLvl = rand(1:10)
        (subseps, hrch, dofs) = LUMND.mnd(A, maxLvl)
        # Check sizes
        @test length(dofs) == maxLvl
        @test length(hrch) == maxLvl
        @test all([length(dofs[i]) == 2^(maxLvl-i) for i = 1:maxLvl])
        @test all([length(hrch[i]) == 2^(maxLvl-i) for i = 1:maxLvl])
        @test length(subseps) == N
        # Check ND consistency
		for (l1, s1) in LUMND.BinTreeIt(maxLvl)
			for (l2, s2) in LUMND.BinTreeIt(maxLvl)
                d1, d2 = dofs[l1][s1], dofs[l2][s2]
                if ! are_connected(l1, s1, l2, s2, maxLvl)
                    @test norm(A[d1,d2]) == 0.0 
                    @test norm(A[d2,d1]) == 0.0
                end
            end
		end
        # Check hierarchy consistency
		for (l, s) in LUMND.BinTreeIt(maxLvl)
            h = hrch[l][s]
            # One root
            @test length(h.nit[1]) == 1
            # Right depth
            @test length(h.nit) == l
            ss = subseps[dofs[l][s]]
            @test all(sort(unique(ss)) .== 1:length(h.nit[end])) # Will pass if ss = []
            for l_ = 1:length(h.nit)
            	# At least one node per level, all numbered 1 ... n
                l__ = length(h.nit[l_])
                @assert l__ > 0
                @test all(1:l__ .== [n.i for n in h.nit[l_]])
            	# At least one children, except at leaves
                for n in h.nit[l_]
                    if l_ < length(h.nit)
                        @assert length(n.c) > 0
                    else
                        @test length(n.c) == 0
                    end
                end
            end
        end
    end
end

function checkAccuracy(A, lvl)
    b = rand(size(A, 1))
    tree = LUMND.factorize(A, lvl)
    x = LUMND.solve(tree, b)
    xref = A\b
    if issymmetric(A)
        k = cond(A, 1)
    else # Crazy slow otherwise
        k = 100
    end
    resref = norm(A*xref-b)/norm(b)
    @test(norm(A*x-b)/norm(b)     < 10 * resref)
    @test(norm(x-xref)/norm(xref) < k  * 1e-12 )
    return b, k, x, xref
end

@testset "Solve-Smallworld-Gen" begin
    for n in 2:20:500
        A = matrixdepot("smallworld", n) 
        leaf = rand([1,2,5])
        N = size(A,1)
        A = A + sparse(UniformScaling(10.0), N, N) # Make it non-singular
        (I,J,V) = findnz(A)
        sig = rand([1e-1, 1e-2, 1e-3, 1e-4]) # Make it non-symmetric
        V += sig .* randn(size(V)) .* V
        A = sparse(I,J,V,N,N)
        lvl = Int64(round(log(N / (leaf)) / log(2) + 1))
        (b, k, x, xref) = checkAccuracy(A, lvl)
        @printf("Smallworld | %3d -> %4.2e %4.2e\n", N, norm(A*x-b)/norm(b), norm(x-xref)/norm(xref))
    end
end

@testset "Solve-Wathen-SPD" begin
    for n in 1:5:30
        A = matrixdepot("wathen", n) 
        leaf = rand([1,2,5])
        N = size(A,1)
        lvl = Int64(round(log(N / (leaf)) / log(2) + 1))
        (b, k, x, xref) = checkAccuracy(A, lvl)
        @printf("Wathen | %5d -> %4.2e %4.2e\n", N, norm(A*x-b)/norm(b), norm(x-xref)/norm(xref))
    end
end

@testset "Solve-Lapl-SPD" begin
    @printf("  # |     N lvl         k ->   |Ax-b|   |x-x*| hash\n")
    for i = 1:100
        d = rand(2:3)
        n = d == 2 ? rand(5:30) : rand(5:10)
        leaf = rand([1,2,5,10,15,20])
        lvl = Int64(round(log(n^d / (leaf)) / log(2) + 1))
        A = Ad(n,d)
        (b, k, x, xref) = checkAccuracy(A, lvl)
        @printf("%3d | %5d  %2d  %2.2e -> %4.2e %4.2e %d\n", i, size(A,1), lvl, k, norm(A*x-b)/norm(b), norm(xref - x)/norm(xref), hash(x))
    end
end

@testset "Solve-Lapl-Gen" begin
    @printf("  # |     N lvl         k ->   |Ax-b|   |x-x*| hash\n")
    for i = 1:100
        d = rand(2:3)
        n = d == 2 ? rand(5:30) : rand(5:10)
        leaf = rand([1,2,5,10,15,20])
        lvl = Int64(ceil(log(n^d / (leaf)) / log(2)))
        A = Ad(n,d)
        N = size(A,1)
        (I,J,V) = findnz(A)
        sig = rand([1e-1, 1e-2, 1e-3, 1e-4]) # Make it non-symmetric
        V += sig .* randn(size(V)) .* V
        A = sparse(I,J,V,N,N)
        (b, k, x, xref) = checkAccuracy(A, lvl)
        @printf("%3d | %5d  %2d  %2.2e -> %4.2e %4.2e %d\n", i, size(A,1), lvl, k, norm(A*x-b)/norm(b), norm(xref - x)/norm(xref), hash(x))
    end
end

@testset "Solve-Blur-SPD" begin
    for n in 3:5:53
        A = matrixdepot("blur", n) 
        leaf = rand([1,2,5])
        N = size(A,1)
        lvl = Int64(round(log(N / (leaf)) / log(2)))
        (b, k, x, xref) = checkAccuracy(A, lvl)
        @printf("Blur | %2d -> %4.2e %4.2e\n", n, norm(A*x-b)/norm(b), norm(x-xref)/norm(xref))
    end
end

