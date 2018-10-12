include("../src/util.jl")
include("../src/LUMND.jl")

using .LUMND
using MatrixMarket

function write_array(io, v)
    if length(v) == 0
        return
    end
    @printf(io, "%d", v[1])
    for v_ in v[2:end]
        @printf(io, ",%d", v_)
    end
end

# Inputs:
#   - A, a sparse matrix in CSC format of size N x N
#   - lvl, the number of level in the ND ordering
#   - prefix, something describing the problem. The 2 files (see below) will be saved to
#       prefix_ord_lvl.txt
#       prefix_clust_lvl.txt
#     where lvl are the # of level (lvl)
# Outputs: none. Write the ordering & clustering to file
#
# Separators are numbered from leaves to top, i.e., for lvl=3, we have 2^3-1=7 separators, numbered 0 to 6
#          6
#        4   5
#       0 1 2 3
#
# Description of the format:
#   prefix_ord_lvl.txt
#       Line 1: lvl #sep
#           where lvl is lvl, and #sep the total number of separators (including leaves), i.e., 2^(lvl)-1
#       Line 2...#sep+1
#           line i corresponds to separator i-2
#           Each line is
#               sepid;i0,i1,i2,...,in
#           where sepid is the separator id [0...2^lvl-2]
#           and i0, i1, ... are the degrees of freedom in the original matrix (i.e., rows or columns of A), from 0 to N-1
#
#   prefix_clust_lvl.txt
#       Line 1: lvl #sep
#           where lvl is lvl, and #sep the total number of separators (including leaves), i.e., 2^(lvl)-1
#       Line 2...#sep+1
#           line i corresponds to separator i-2
#           A separator at lvl l as a hierarchy with l+1 levels
#           The hierarchy can be seen as a sequence of nested interval, where the top (root) one is 'everything'
#           Each level as at least one interval, but potentially empty
#           Each line contains
#               sepid;interval 0;interval 1;interval 2;...;interval k
#           Interval 0 indicates how we 'group' the dofs (described in prefix_ord_lvl)
#           For instance, if the dofs are
#                       4,7,2,6,11,15,17,25,41,12
#           And interval 0 is
#                       0,3,4,8,10
#           This means we have 4 clusters at the leaf-level:
#                       4,7,2 - 6 - 11,15,17,25 - 41,12
#           Then if inteval 1 is
#                       0,2,4
#           It means at the next level, the clusters become:
#                       4,7,2,6 - 11,15,17,25,41,12
#           Then finally if interval 2 is
#                       0,2
#           It means at the next level, the clusters become:
#                       4,7,2,6,11,15,17,25,41,12
#
# Everything is stored in base-0 by convention

function order_and_write_file(A, lvl, prefix)

    # Compute ordering
    (subseps, hrch, dofs) = LUMND.mnd(A+A', lvl) # Works on the sparsity pattern of A+A'

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
    ordering_name = @sprintf("%s_ord_%d.txt", prefix, lvl)
    clustering_name = @sprintf("%s_clust_%d.txt", prefix, lvl)
    @printf("Data saved to %s and %s\n", ordering_name, clustering_name)

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

# We can also finally save the matrix
matfile = @sprintf("%s.mtx", prefix)
mmwrite(matfile, A)
@printf("Matrix saved in MatrixMarket format to %s\n", matfile)
