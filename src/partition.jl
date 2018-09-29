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

IP = NamedTuple{(:p, :l, :r, :a),Tuple{Int64,Int64,Int64,Int64}}
# Parent / Left / Right / Anchor lvl

@inline function lt_ip(x::IP, y::IP)
    return      (x.p <  y.p) ||
              ( (x.p == y.p) && (x.l <  y.l) ) ||
              ( (x.p == y.p) && (x.l == y.l) && (x.r <  y.r) ) ||
              ( (x.p == y.p) && (x.l == y.l) && (x.r == y.r) && (x.a < y.a) )
end

ord = Base.ord(isless,identity,false,Base.Forward)

# Fills colptr and rowval (assumed large enough) with A[list, list], excluding the diagonal
# List is assumed to be sorted
# base-0
# Ignore the diagonal
function fill_nodiag_list_sorted!(A::SparseMatrixCSC{Tv,Ti}, colptr::Array{Tj}, rowval::Array{Tj}, list::AbstractArray{Tk}) where {Ti <: Integer, Tj <: Integer, Tk <: Integer, Tv}
    n = length(list)
    colptr[1] = 1
    if n == 0 
        return
    end
    for j = 1:n
        colptr[j+1] = colptr[j] # Initially, nothing
        j_ = list[j]
        i_ = 1
        for k = A.colptr[j_]:(A.colptr[j_+1]-1)
            i = A.rowval[k]
            if i == j_ # Skip diagonal
                continue
            end
            i_ = searchsortedfirst(list, Tk(i), i_, length(list), ord) # Ok since sorted
            if i_ < n+1 && list[i_] == Tk(i) # Yes, at i_
                rowval[colptr[j+1]] = i_
                colptr[j+1] += 1
            end
        end
    end
    nz = colptr[n+1]-1
    colptr[1:n+1] .-= 1
    rowval[1:nz]  .-= 1
    return
end

##
## METIS
##

const metis_options = - ones(Cint, Metis.METIS_NOPTIONS)

function vertex_sep_fast!(n, colptr::Array{Cint,1}, rowval::Array{Cint,1}, part::Array{Cint,1})
    sepSize = zeros(Cint,1)
    n2 = Cint(n)
    err = ccall((:METIS_ComputeVertexSeparator,Metis.libmetis), Cint,
                (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint},     Ptr{Cint}, Ptr{Cint}),
                Ref(n2),    colptr,    rowval,    C_NULL,    metis_options, sepSize,   part)
    err == Metis.METIS_OK || error("METIS_ComputeVertexSeparator returned error code $err")
    return
end

function vertex_sep_bissect_fast!(n, colptr::Array{Cint,1}, rowval::Array{Cint,1}, part::Array{Cint,1})
    nparts = Cint(2)
    one = Cint(1)
    two = Cint(2)
    n2 = Cint(n)
    objval = zeros(Cint,1)
    err = ccall((:METIS_PartGraphRecursive,Metis.libmetis), Int32,
                (Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
                 Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32},
                 Ptr{Int32}, Ptr{Int32}, Ptr{Int32}),
                Ref(n2), Ref(one), colptr, rowval, C_NULL, C_NULL, C_NULL, Ref(two),
                C_NULL, C_NULL, metis_options, objval, part)
    err == Metis.METIS_OK || error("METIS_PartGraphRecursive returned error code $err")
    for i1 = 1:n # update part: those with (0), connected to at least one (1), gets a (2)
        if part[i1] == 1
            continue
        end
        for k0 = colptr[i1]:colptr[i1+1]-1 # 0-based
            k1 = k0 + 1
            j0 = rowval[k1]
            j1 = j0 + 1
            if part[j1] == 1
                part[i1] = 2 # Becomes a sep
                break
            end
        end
    end
end

##
## Partitioning
##

const IP0 = (p=1,l=0,r=0,a=0)

# Partition the graph represented by the sparse matrix A
# from level 1 to maxLevel
# Interiors at level maxLevel are 2^(maxLevel-1)
# Separators at level 1 <= l <= maxLevel-1 are 2^(l-1)
#
# Separators & interiors are numbered
#                       1                 <-- root level 1
#               2               3
#           4       5      6        7     <-- level l, numbered 2^(l-1):(2^l-1)
#         8  9    10 11  12 13    14 15   <-- leaves,  numbered 2^(maxLevel-1):(2^maxLevel-1)
# Total number of clusters = 2^(maxLevel)-1
#
# Each separator is divided in a hierarchy of partitions
# Separator at level 1 <= l <= maxLevel-1 consists of a hierarchy of maxLevel-l+1 nested partitions
# Interiors (leaf level) have 1 nested partition
#                       1                 <-- root lvl 1,               maxLevel     hierarchy
#               2               3         <-- level l,                  maxLevel-l+1 hierarchy
#           4       5      6        7     <-- maxLevel-1,               2-hierarchy
#         8  9    10 11  12 13    14 15   <-- leaves lvl maxLevel,      1-hierarchy
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
# In addition, each node holds a boolean indicating wether its partition is to be considered an anchor or not
#
# Outputs
#       * subseps       subseps[i] = 1 ... #partitions at leaf level for dof i if in a separator ; 1 only for the leaves as they have a 1-hierarchy, always
#       * hrch          hrch[l][s] a (maxLevel-l+1)-levels hierarchy for a separator at level l
#       * dofs          dofs[l][s] the dofs of separator s at level l
function mnd(A::SparseMatrixCSC{Tv,Ti}, maxLevel ; X::Union{Nothing,AbstractArray{Float64,2}}=nothing, verbose::Bool=false) where {Tv, Ti}

    if verbose
        if X == nothing
            println("Algebraic partitionning of graph with $(size(A,1)) vertices, $maxLevel levels")
        else
            println("Geometric partitionning of graph with $(size(A,1)) vertices, $maxLevel levels, in $(size(X,1))D")
        end
    end
    
    # Sep     = (p, l, r, 0   )
    # Leaf    = (p, 0, 0, 0   )
    # Anchors = (p, l, r, lvl )  (anchor is a sep) (anchor starting at level lvl)
    n = size(A,1)
    parts = Array{IP,1}(undef, n)
    fill!(parts,IP0)

    # Prepare ints
    ints  = [ collect(Cint, 1:n)           ] # [   (id, int)     ]
    edges = [ Tuple{Int64,Array{Cint,1}}[] ] # [ [ (side, int) ] ]  1 = left, 2 = right

    # Temp arrays
    colptr = Array{Cint,1}(undef, n+1)
    rowval = Array{Cint,1}(undef, nnz(A))
    pids   = Array{Cint,1}(undef, n)

    # Tag vector
    tag = zeros(Int64, n)

    # Output structures
    L           = (2^maxLevel)-1
    sepk        = 1
    out_dofs    = Array{Array{Int64,1},1}(undef, L)    
    out_hrch    = Array{IT,1}(undef, L)

    # Initial hierarchy
    for i = 1:L
        out_hrch[i] = IT(1)
    end

    for lvl = 1:maxLevel

        if verbose println("Level $lvl") end

        new_ints  = Array{ Array{Cint,1},                       1}(undef, 2^lvl)
        new_edges = Array{ Array{Tuple{Int64,Array{Cint,1}},1}, 1}(undef, 2^lvl)

        ## Break interiors in left/right/mid and update edges (left, right) ids to define hierarchy
        for (id, (dofs_int, edges_sep)) = enumerate(zip(ints, edges))
            ## ND ordering
            list = sort(vcat(dofs_int, map(x -> x[2], edges_sep)...))
            nl = length(list)
            if lvl == maxLevel # Leaf level
                left, right, mid = Cint[], Cint[], list
            elseif length(list) <= 1 # Really isn't anything to separate
                left, right, mid = Cint[], Cint[], list
            else 
                if X == nothing
                    fill_nodiag_list_sorted!(A, colptr, rowval, list)
                    vertex_sep_bissect_fast!(nl, colptr, rowval, pids)
                else
                    vertex_sep_geo_fast!(A, X, list, pids, tag)
                end
                left   = list[pids[1:nl] .== 0]
                right  = list[pids[1:nl] .== 1]
                mid    = list[pids[1:nl] .== 2]
            end
            tag[left]  .= 2*id-1
            tag[right] .= 2*id
            tag[mid]   .= 0
            # Update edges
            for (side, dofs) in edges_sep
                for i in dofs
                    plra = parts[i]
                    if tag[i] == 0 # Anchor in this sep : make it anchor (override potential previous anchors)
                        if side == 1
                            parts[i] = (p=plra.p, l=id, r=0,  a=lvl)
                        else
                            parts[i] = (p=plra.p, l=0,  r=id, a=lvl)
                        end
                    elseif plra.a == 0 # Not anchor before, Not in this sep : regular step
                        if side == 1
                            parts[i] = (p=plra.p, l=tag[i], r=plra.r, a=plra.a)
                        else
                            parts[i] = (p=plra.p, l=plra.l, r=tag[i], a=plra.a)
                        end
                    end # else tag[i] != 0 && plra.a > 1 : Anchor before, Not in this sep : skip
                end
            end
            # Figure out new sep & left/right interiors
            mid_int = filter(i -> parts[i] == IP0, mid)
            fill!(view(parts, mid_int), (p=1,l=2*id-1,r=2*id,a=0))
            # Update int
            left_int  = filter(i -> parts[i] == IP0, left)
            right_int = filter(i -> parts[i] == IP0, right)
            # What is left to work on
            new_ints[2*id-1] = left_int
            new_ints[2*id]   = right_int
            # Output data
            out_dofs[sepk]    = mid_int
            # Initialize empty boundary
            new_edges[2*id-1] = Tuple{Int64,Array{Cint,1}}[]
            new_edges[2*id]   = Tuple{Int64,Array{Cint,1}}[]
            # Next 
            sepk += 1
        end

        # Go through edges and create separators hierarchy
        if lvl < maxLevel
            sepk2 = 1
            for lvl2 = 1:lvl
                for sep2 = 1:2^(lvl2-1)
                    # Get sep
                    dofs = out_dofs[sepk2]
                    # Get unique values, assign each bloc to its interior
                    vals = parts[dofs] # (parent, left, right, anchor lvl)
                    lr_srt = sortperm(vals, lt=lt_ip) # Sort by parent first - arbitrary after
                    vals_srt = vals[lr_srt]
                    dofs_srt = dofs[lr_srt]
                    u_ptr = sorted_to_ptr(vals_srt)
                    # Create hierarchy
                    it     = out_hrch[sepk2]
                    push!(it.nit, NodeIT[]) # Create new level
                    if length(dofs) == 0
                        @assert length(it.nit[end-1]) == 1
                        x = NodeIT(1)
                        x.p = it.nit[end-1][1]
                        push!(it.nit[end-1][1].c,x)
                        push!(it.nit[end], x)
                    else
                        for p = 1:length(u_ptr)-1
                            dofs_p = dofs_srt[u_ptr[p]:u_ptr[p+1]-1] # The dofs - this is non-empty by definition of unique 
                            old_p  = parts[dofs_p[1]].p # All same since we sorted by parent 
                            left   = parts[dofs_p[1]].l
                            right  = parts[dofs_p[1]].r
                            anch   = parts[dofs_p[1]].a
                            # Update partition
                            for i in dofs_p
                                parts[i] = (p=p, l=left, r=right, a=anch)
                            end
                            # Create node & update hierarchy
                            x     = NodeIT(p)
                            x.p   = it.nit[end-1][old_p]
                            x.anc = (anch > 0)
                            push!(it.nit[end-1][old_p].c, x)
                            push!(it.nit[end], x)
                            # If not an anchor, push in new_edges, otherwise won't be touched again
                            if anch == 0 
                                push!(new_edges[left],  (1, dofs_p))
                                push!(new_edges[right], (2, dofs_p))
                            end
                        end
                    end
                    sepk2 += 1
                end
            end
            @assert sepk2-1 == 2^lvl-1
        end

        ints = new_ints
        edges = new_edges

    end

    out_subseps = map(x -> x.p, parts)

    @assert sepk-1 == L

    out_dofs_2 = Vector{Vector{Vector{Int64}}}(undef, maxLevel)
    out_hrch_2 = Vector{Vector{IT}}(undef, maxLevel)

    sepk = 1
    for lvl = 1:maxLevel
        nsep = 2^(lvl-1)
        out_hrch_2[maxLevel-lvl+1] = out_hrch[sepk:(sepk+nsep-1)]
        out_dofs_2[maxLevel-lvl+1] = out_dofs[sepk:(sepk+nsep-1)]
        sepk += nsep
    end

    return out_subseps, out_hrch_2, out_dofs_2

end
