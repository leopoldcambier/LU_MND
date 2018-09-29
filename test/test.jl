include("../src/util.jl")
include("../src/LUMND.jl")

using .LUMND
using Test
using Random
using SHA 

Random.seed!(0)

@testset "Solve" begin
    @printf("  # |     N lvl         k ->   |Ax-b|   |x-x*| hash\n")
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
        @printf("%3d | %5d  %2d  %2.2e -> %4.2e %4.2e %d\n", i, size(A,1), lvl, k, norm(A*x-b)/norm(b), norm(xref - x)/norm(xref), hash(x))
        @test(norm(A*x-b)/norm(b)     <     1e-13)
        @test(norm(x-xref)/norm(xref) < k * 1e-13)
    end
end

@testset "Partition" begin
    for i = 1:10
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
            # One root, not an anchor
            @test length(h.nit[1]) == 1
            @test h.nit[1][1].anc == false
            # Check anchor invariants
            n = h.nit[1][1]
            stack = [n]
            while length(stack) > 0
                n = pop!(stack)
                if length(n.c) > 0
                    if n.anc
                        @test all([c.anc for c in n.c])
                    end
                end
                for c in n.c
                    push!(stack, c)
                end
            end
        end
    end
end
