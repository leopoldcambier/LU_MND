Cluster = Tuple{Int64,Int64,Int64}

function lt_cl(x::Cluster, y::Cluster)
    return     (x[1] < y[1]) ||
             ( (x[1] == y[1]) && (x[2] < y[2]) ) ||
             ( (x[1] == y[1]) && (x[2] == y[2]) && (x[3] < y[3]) )
end

mutable struct NDSep # A nested dissection separator

    ## Hierarchy & dofs data
	ptr::Array{Int64,1} 			# A pointer array to the end of each cluster within s at the current stage of the elimination
    hrch                            # A hierarchy, i.e., a tree of integer, monotone left -> right
    start::Int64                    # Where does the separator start
    depth::Int64                    # Where the clustering stand in terms of depth

    ## Data
	Low::Matrix{Float64}            # Ass & Ans, Pivot + Lower-part

    ## Block information, basically pointers to the data above
    ALow::Matrix{StridedMatrix{Float64}}         # ALow[n,i] = ith partition, nth neighbor, below diagonal
    Nbr::Vector{Cluster}                         # Nbr[n]    = nth neighbor (lvl, sep, part)

    ## Range information, for easy backward/forward pass
    NbrRanges::Vector{UnitRange{Int64}}

    function NDSep(start, subs, hrch) 
        if length(subs) >= 1
            ptr = sorted_to_ptr(subs)
        else
            ptr = [1, 1]
        end
        this = new()
        this.start = start
        @assert length(ptr) >= 2
        this.ptr = ptr
        this.hrch = hrch
        this.depth = length(this.hrch.nit)
        this.NbrRanges = Vector{UnitRange{Int64}}()
        this
    end
end

struct NDTree
    t::Vector{Vector{NDSep}}
    p::Vector{Int64}
end


function get_dofs(self::NDSep, i::Int64)
    return (self.start-1) .+ (self.ptr[i]:self.ptr[i+1]-1)
end

function get_dofs(self::NDSep)
    return (self.start-1) .+ (self.ptr[1]:self.ptr[end]-1)
end

function get_size(self::NDSep)
    return (self.ptr[end] - self.ptr[1])
end

function get_size(self::NDSep, i::Int64)
    return (self.ptr[i+1] - self.ptr[i])
end

function nparts(self::NDSep)
    np = length(self.ptr)-1;
    @assert length(self.hrch.nit[self.depth]) == np
    return np
end

function nparts_merged(self::NDSep)
    return length(self.hrch.nit[self.depth-1])
end

function ptr_merged(self::NDSep)
    ptr = self.ptr
    ptr_new = ones(Int64, nparts_merged(self)+1);
    for n in self.hrch.nit[self.depth-1]
        ptr_new[n.i+1] = ptr[n.c[end].i+1]
    end
    return ptr_new
end

# Given i, go up l level in the tree
# Return the range low:high (i \in low:high) 
# corresponding to i's ancestor at i's level
function go_up_down(self::NDSep, i::Int64, l::Int64)
    @assert l >= 0
    if l == 0
        return i:i
    end
    # Go up
    h = self.hrch.nit
    p = h[end][i]
    for _ in 1:l
        p = p.p
    end
    # Go left
    left = p
    for _ in 1:l
        left = left.c[1]
    end
    # Go right
    right = p
    for _ in 1:l
        right = right.c[end]
    end
    return left.i:right.i
end

function parent(self::NDSep, i::Int64)
    self.hrch.nit[self.depth][i].p.i
end

function siblings(self::NDSep, i::Int64) 
    p = self.hrch.nit[self.depth][i].p
    return (p.c[1].i:p.c[end].i, p.i)
end

# BinTree iterator
# Generate
# (1, 1), (1, 2), (1, 3), (1, 4)
#     (2, 1),         (2, 2)
#             (3, 1)
struct BinTreeIt
    maxLevel::Int64
end

function Base.iterate(iter::BinTreeIt, state=(1, 0))
    (l, s) = state
    if l >= iter.maxLevel
        return nothing
    end
    if s >= 2^(iter.maxLevel-l)
        next = (l+1, 1)
    else
        next = (l, s+1)
    end
    return (next, next)
end

function factorize(A, maxLevel)
   
    N = size(A, 1)

    ## Partition
    (seps, subseps, hrch, dofs) = ml_nd_hrch_fast(A, maxLevel, verbose=false);

	## Permute the matrix & initialize an empty tree
	perm = Vector{Int64}(undef, N)
    Tree = Vector{Vector{NDSep}}(undef, maxLevel)
    start = 1
    for lvl = 1:maxLevel
        nsep = length(dofs[lvl])
        Tree[lvl] = Vector{NDSep}(undef, nsep)
        for sep = 1:nsep
            dofs_ = dofs[lvl][sep]
            subs_ = subseps[dofs_]
            perm_ = sortperm(subs_)
            perm[start:(start+length(dofs_)-1)] .= dofs_[perm_]
            Tree[lvl][sep] = NDSep(start, subs_[perm_], hrch[lvl][sep])
            start += length(dofs_)
        end
    end
    A = A[perm,perm]

    ## Create dofs2cl map
    dofs2cl = Vector{Cluster}(undef, N)
    for (lvl, sep) in BinTreeIt(maxLevel)
        for i in dofs[lvl][sep]
            dofs2cl[i] = (lvl, sep, subseps[i])
        end
    end
    dofs2cl = dofs2cl[perm]

    ## Compute sparsity pattern, allocate & assemble
    for (lvl, sep) in BinTreeIt(maxLevel)
        h = Tree[lvl][sep].hrch
        s = length(dofs[lvl][sep])
        np = length(h.nit[end])
        # Get the original neighbors
        (I, t1, t2) = findnz(A[:,get_dofs(Tree[lvl][sep])])
        nbrs = sort(unique(I))
        Nbr = Set(dofs2cl[nbrs])
        # Add the children
        if lvl > 1
            union!(Nbr, Tree[lvl-1][2*sep-1].Nbr)
            union!(Nbr, Tree[lvl-1][2*sep  ].Nbr)
        end
        # Remove too small
        filter!(x -> x[1] >= lvl, Nbr)
        # Add self
        union!(Nbr, [(lvl, sep, i) for i in 1:np])
        # Subsample
        Nbr2 = copy(Nbr)
        for (l_, s_, i_) in Nbr2
            lr = go_up_down(Tree[l_][s_], i_, lvl - 1)
            for i2_ in lr
                push!(Nbr, (l_, s_, i2_))
            end
        end
        # Clean & sort
        Nbr = collect(Nbr);
        sort!(Nbr, lt=lt_cl)
        Tree[lvl][sep].Nbr = Nbr
        # Allocate full data
        nsizes = map(x -> Tree[x[1]][x[2]].ptr[x[3]+1]-Tree[x[1]][x[2]].ptr[x[3]], Nbr)
        n = sum(nsizes)
        Low = zeros(n, s)
        Tree[lvl][sep].Low = Low
        # Get subblocks
        ALow = Matrix{StridedMatrix{Float64}}(undef, (length(Nbr), np))
        sptr = Tree[lvl][sep].ptr
        nptr = sizes_to_ptr(nsizes)
        for i in 1:np
            for n in 1:length(Nbr)
                # Define subblock
                rows = nptr[n]:(nptr[n+1]-1)
                cols = sptr[i]:(sptr[i+1]-1)
                ALow[n,i] = view(Low, rows, cols);  
                # Fill subblock
				(lvln,sepn,in) = Nbr[n]
				rows = get_dofs(Tree[lvln][sepn], in)
				cols = get_dofs(Tree[lvl][sep], i)
				ALow[n,i][:,:] = A[rows, cols]
            end
        end
        Tree[lvl][sep].ALow = ALow
    end

    ## Factorize
    for lvl = 1:maxLevel

        
        ## Eliminate

        for sep = 1:length(dofs[lvl])
            s = Tree[lvl][sep]
            @show lvl, sep, maxLevel
            @show s.ptr
            @show s.Nbr
            @show size(s.ALow)
            @assert size(s.ALow, 2) == 1
            @assert s.Nbr[1] == (lvl, sep, 1)
            # POTF, factor pivot
            Ass = s.ALow[1,1]
            @show size(Ass)
            LAPACK.potrf!('L', Ass)
            push!(s.NbrRanges, get_dofs(s))
            # TRSM, solve panel in place
            for n = 2:length(s.Nbr)
                Ans = s.ALow[n,1]
                (ln,sn,in) = s.Nbr[n]
                push!(s.NbrRanges, get_dofs(Tree[ln][sn],in))
                BLAS.trsm!('R','L','T','N', 1.0, Ass, Ans)
            end
            # GEMM, update schur complement
            for ni = 2:length(s.Nbr)
                for nj = ni:length(s.Nbr)
                    Ais = s.ALow[ni,1]
                    Ajs = s.ALow[nj,1]
                    (li,si,ii) = s.Nbr[ni]
                    n = Tree[li][si]
                    jj = searchsortedfirst(n.Nbr, s.Nbr[nj], lt=lt_cl)
                    Aji = n.ALow[jj,ii]
                    BLAS.gemm!('N', 'T', - 1.0, Ajs, Ais, 1.0, Aji) 
                end
            end
        end

        ## Merge

        # Change the data layout
        for lvl2 = lvl+1:maxLevel
            for sep2 = 1:length(dofs[lvl2])
                # Take care of the columns
                s = Tree[lvl2][sep2]
                ptr2 = ptr_merged(s)
                # Take care of the rows
                n_i1, n_b1, n_merged = 1, 1, 1
                while n_b1 <= length(s.Nbr)
                    # What are n's siblings ?
                    (l_, s_, p_) = s.Nbr[n_b1]
                    (n_siblg, parent) = siblings(Tree[l_][s_], p_)
                    # All the siblings should follow
                    n_b2 = n_b1
                    n_i2 = n_i1
                    for i = 1:length(n_siblg)
                        n_p_  = n_siblg[i]
                        @assert s.Nbr[n_b2] == (l_, s_, n_p_)
                        n_b2 += 1
                        n_i2 += get_size(Tree[l_][s_], n_p_)
                    end
                    # Update the row in Low and in Nbr
                    s.Nbr[n_merged] = (l_, s_, parent)
                    for c = 1:length(ptr2)-1
                        s.ALow[n_merged, c] = view(s.Low, n_i1:(n_i2-1), ptr2[c]:(ptr2[c+1]-1))
                    end
                    # Next
                    n_merged += 1
                    n_i1 = n_i2
                    n_b1 = n_b2
                end
                # Resize
                s.ALow = s.ALow[1:n_merged-1, 1:length(ptr2)-1]
                s.Nbr  = s.Nbr[1:n_merged-1]
            end
        end

        # Update ptr & depth
        for lvl2 = lvl+1:maxLevel
            for sep2 = 1:length(dofs[lvl2])
                s = Tree[lvl2][sep2]
                s.ptr = ptr_merged(s)
                s.depth -= 1
            end
        end
    end

    return NDTree(Tree, perm)

end

function solve(tree::NDTree, x::Vector{Float64})
    y = x[tree.p]
    maxLevel = length(tree.t)
    # Forward
    for l = 1:maxLevel
        for s = 1:2^(maxLevel-l)
            Lpp = tril(tree.t[l][s].ALow[1,1])
            p   = tree.t[l][s].NbrRanges[1]
            # Pivot
            y[p] = Lpp \ y[p]
            # Panel
            for (Lnp, n) in zip(tree.t[l][s].ALow[2:end,1], tree.t[l][s].NbrRanges[2:end])
                y[n] -= Lnp * y[p]
            end
        end
    end
    # Backward
    for l = maxLevel:-1:1
        for s = 1:2^(maxLevel-l)
            Lpp = tril(tree.t[l][s].ALow[1,1])
            p   = tree.t[l][s].NbrRanges[1]
            # Panel
            for (Lnp, n) in zip(tree.t[l][s].ALow[2:end,1], tree.t[l][s].NbrRanges[2:end])
                y[p] -= Lnp' * y[n]
            end
            # Pivot
            y[p] = Lpp' \ y[p]
        end
    end
    return y[invperm(tree.p)]
end

