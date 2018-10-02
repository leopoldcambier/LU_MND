##
## Data Structures
##

# Each separator is associated with a tree like
#
#       1
#       |
#       1
#     / | \
#    1  2  3
#   / \ |  |
#   1 2 3  4
#
# The ids are increasing left-right

# A node a the tree
mutable struct NodeIT 
    i::Int64                        # Id
    anc::Bool                       # Wether the point is an anchor
    c::Array{NodeIT,1}              # Childrens
    p::Union{Nothing, NodeIT}       # Parent
    function NodeIT(i::Int64)
        this = new()
        this.i = i
        this.c = NodeIT[]
        this.p = nothing
        this.anc = false
        this
    end
end
Base.show(io::IO, n::NodeIT) = print(io, n.anc ? @sprintf("%d*",n.i) : @sprintf("%d^",n.i))

# The tree
mutable struct IT 
    nit::Array{Array{NodeIT,1},1}       
    function IT(i::Int64)               
        this = new()
        this.nit = [ [ NodeIT(i) ] ]
        this
    end

    function IT()
        this = new()
        this.nit = [ [] ]
        this
    end
end

function show_it(it::IT)
    for l = 1:length(it.nit)
        for n in it.nit[l]
            for c in n.c
                @printf("%s%d-%d (~>%d ; %s)\n", repeat(" ", l), n.i, c.i, c.p.i, n.anc ? "anc" : "--")
            end
        end
    end
end

IP = NamedTuple{(:p, :l, :r),Tuple{Int64,Int64,Int64}} # Parent/Left/Right

function lt_ip(x::IP, y::IP)
    return      (x.p <  y.p) ||
              ( (x.p == y.p) && (x.l <  y.l) ) ||
              ( (x.p == y.p) && (x.l == y.l) && (x.r <  y.r) )
end

ord = Base.ord(isless,identity,false,Base.Forward)

const IP0 = (p=1,l=1,r=1)

# Partition the graph represented by the sparse matrix A
# from level 1 to maxLevel
#
# Each separator is divided in a hierarchy of partitions
# Separator at level 1 <= l <= maxLevel-1 consists of a hierarchy of maxLevel-l+1 nested partitions
# Interiors (leaf level) have 1 nested partition
# Hierarchy is encoded into a integer-tree representing the nested hierarchy
# For instance, a given separator's (lvl maxLevel-2) dofs are numbered
#       subseps[dofs[lvl][sep]] = [1 1 2 3 3 4 5 5 5]
# with a tree stored in 
#       hrch[lvl][sep] =        1
#                             /   \
#                            1     2
#                          / | \   | \
#                         1  2  3  4  5
# indicating a 3-levels hierarchy and the way the merging should be performed
# The nodes are numbered 1 ... |p| and hold a value corresponding to their id in IT.nit, i.e., a node at level l with value i is hold at nit[l][i]
# In particular, the hierarchy is
#                       [1 1 2 3 3 4 5 5 5] 5 partitions
#                       [1 1 1 1 1 2 2 2 2] 2 partitions
#                       [1 1 1 1 1 1 1 1 1] 1 partition - root
#
# Outputs
#       * subseps       subseps[i] = 1 ... #partitions at leaf level for dof i if in a separator ; 1 only for the leaves as they have a 1-hierarchy, always
#       * hrch          hrch[l][s] a (maxLevel-l+1)-levels hierarchy for a separator at level l
#       * dofs          dofs[l][s] the dofs of separator s at level l
function mnd(A::SparseMatrixCSC{Tv,Ti}, maxLevel ; verbose::Bool=false) where {Tv, Ti}

    @assert maxLevel >= 1
    if verbose
        println("Algebraic partitionning of graph with $(size(A,1)) vertices, $maxLevel levels")
    end
    
    n = size(A,1)
    parts = Array{IP,1}(undef, n) # (parent, left, right) for each dof
    fill!(parts,IP0)

    # Prepare stuff to separate
    tosep  = [ collect(Cint, 1:n) ]

    # Temp arrays for Metis
    colptr = Vector{Cint}(undef, n+1)
    rowval = Vector{Cint}(undef, nnz(A))
    pids   = Vector{Cint}(undef, n)

    # Output structures
    dofs    = Vector{Vector{Vector{Int64}}}(undef, maxLevel)
    hrch    = Vector{Vector{IT}}(undef, maxLevel)

    # Do the work
    for lvl = 1:maxLevel

        dofs[lvl] = Vector{Vector{Int64}}(undef, 2^(lvl-1))
        hrch[lvl] = Vector{IT}(undef, 2^(lvl-1))
        for sep = 1:2^(lvl-1)
            hrch[lvl][sep] = IT(1)
        end
        tosep_next = Vector{Vector{Int64}}(undef, 2^lvl)

        if verbose println("Level $lvl, partitioning") end

        for (sep, list) in enumerate(tosep)

            # What is our separator -unique- id
            id = 2^(lvl-1) + (sep - 1)

            # Do the ordering, i.e., cut A[list,list]
            if lvl == maxLevel # Leaf level
                pids[1:length(list)] .= 2
            elseif length(list) <= 1 # Really isn't anything to separate
                pids[1:length(list)] .= 2
            else 
                fill_nodiag_list_sorted!(A, colptr, rowval, list)
                vertex_sep_bissect_fast!(length(list), colptr, rowval, pids)
            end

            # Update parts, so that each node has a proper (left, right) pair
            for (i, side) in zip(list, pids)
                p = parts[i]
                if side == 2 && p.r == id && p.l == id # Separator, brand new
                    p = (p=1, l=2*id, r=2*id+1)
                elseif side == 0 # In the left interior
                    if p.l == id && p.r == id # Interior, brand new
                        p = (p=1, l=2*id, r=2*id)
                    elseif p.l == id # Upper level separator, on its left
                        p = (p=p.p, l=2*id, r=p.r) 
                    elseif p.r == id # Upper level separator, on its right
                        p = (p=p.p, l=p.l, r=2*id) 
                    end
                elseif side == 1 # In the right interior
                    if p.l == id && p.r == id # Interior, brand new
                        p = (p=1, l=2*id+1, r=2*id+1)
                    elseif p.l == id # Upper level separator, on its left
                        p = (p=p.p, l=2*id+1, r=p.r) 
                    elseif p.r == id # Upper level separator, on its right
                        p = (p=p.p, l=p.l, r=2*id+1) 
                    end
                end
                parts[i] = p
            end

            # Find new sep & create initial hierarchy
            mid = findall(i -> parts[i] == (p=1, l=2*id, r=2*id+1), list)
            dofs[lvl][sep] = list[mid]

            # Define new interiors + boundaries
            left  = findall(i -> parts[i].l == 2*id  , list)
            right = findall(i -> parts[i].r == 2*id+1, list)
            tosep_next[2*sep-1] = list[left]
            tosep_next[2*sep  ] = list[right]
        end

        # Go through edges and update hierarchy & parents
        if lvl < maxLevel
            for lvl2 = 1:lvl
                for sep2 = 1:2^(lvl2-1)
                    vals  = parts[dofs[lvl2][sep2]]
                    lr_srt = sortperm(vals, lt=lt_ip) # Sort by parent first - arbitrary after
                    vals_srt = vals[lr_srt]
                    dofs_srt = dofs[lvl2][sep2][lr_srt]
                    u_ptr = sorted_to_ptr(vals_srt)
                    # Create hierarchy
                    it = hrch[lvl2][sep2]
                    push!(it.nit, NodeIT[]) # Create new level
                    if length(dofs_srt) == 0
                        @assert length(it.nit[end-1]) == 1
                        x = NodeIT(1)
                        x.p = it.nit[end-1][1]
                        push!(it.nit[end-1][1].c,x)
                        push!(it.nit[end], x)
                    else
                        for p = 1:length(u_ptr)-1
                            dofs_p = dofs_srt[u_ptr[p]:u_ptr[p+1]-1] # The dofs - this is non-empty by definition of unique 
                            old_p  = parts[dofs_p[1]].p # All same since we sorted by parent 
                            left   = parts[dofs_p[1]].l # ..
                            right  = parts[dofs_p[1]].r # ..
                            # Update partition
                            for i in dofs_p
                                parts[i] = (p=p, l=left, r=right)
                            end
                            # Create node & update hierarchy
                            x     = NodeIT(p)
                            x.p   = it.nit[end-1][old_p]
                            push!(it.nit[end-1][old_p].c, x)
                            push!(it.nit[end], x)
                        end
                    end
                end
            end
        end

        # Next clusters to consider
        tosep = tosep_next
    end

    subseps = map(x -> x.p, parts)
    reverse!(hrch)
    reverse!(dofs)

    return subseps, hrch, dofs

end
