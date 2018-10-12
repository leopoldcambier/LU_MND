include("../src/util.jl")
include("../src/LUMND.jl")

using .LUMND

function write_array(io, v)
    if length(v) == 0
        return
    end
    @printf(io, "%d", v[1])
    for v_ in v[2:end]
        @printf(io, ",%d", v_)
    end
end

function order_and_write_file(A, lvl, file_prefix)
    # Compute ordering
    (subseps, hrch, dofs) = LUMND.mnd(A, lvl)

    # Sort dofs
    for l = 1:lvl
        for s = 1:length(dofs[l])
            d_ = dofs[l][s]
            s_ = subseps[d_]
            p_ = sortperm(s_)
            dofs[l][s] = dofs[l][s][p_]
        end
    end

    # Write
    ordering_name = @sprintf("%s_ord_%d.txt", file_prefix, lvl)
    clustering_name = @sprintf("%s_clust_%d.txt", file_prefix, lvl)

    # Write ordering
    file = open(ordering_name, "w")
    @printf(file, "%d %d\n", lvl, 2^lvl-1) # Lvl, # seps
    k = 0 # Base-0
    for l = 1:lvl
        for s = 1:length(dofs[l])
            d_ = dofs[l][s]
            print(file, @sprintf("%d;", k))
            write_array(file, dofs[l][s] .- 1) # Base-0
            println(file)
            k += 1
        end
    end
    close(file)

    # Write partitionning
    file = open(clustering_name, "w")
    @printf(file, "%d %d\n", lvl, 2^lvl-1) # Lvl, # seps
    k = 0 # Base-0
    for l = 1:lvl
        for s = 1:length(dofs[l])
            d_ = dofs[l][s]
            s_ = subseps[d_]
            h_ = hrch[l][s]
            # The leaf-level pointers
            old_ptr = sorted_to_ptr(s_)
            print(file, @sprintf("%d;", k));
            write_array(file, old_ptr .- 1) # Base-0
            print(file, ";")
            # The hierarchy of intervals
            for d = length(h_.nit)-1:-1:1
                ptr = [n.c[1].i for n in h_.nit[d]]
                push!(ptr, length(h_.nit[d+1]) + 1)
                write_array(file, ptr .- 1) # Base-0
                print(file,";")
                old_ptr = ptr
            end
            println(file)
            k += 1
        end
    end
    close(file)

end

n = 20
d = 2
A = Ad(n,d) # n^d laplacian
lvl = 5
prefix = @sprintf("lapl_%d_%d", n, d)
order_and_write_file(A, lvl, prefix)
