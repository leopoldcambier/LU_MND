Cluster = Tuple{Int64,Int64,Int64}

function lt_cl(x::Cluster, y::Cluster)
    return     (x[1] < y[1]) ||
             ( (x[1] == y[1]) && (x[2] < y[2]) ) ||
             ( (x[1] == y[1]) && (x[2] == y[2]) && (x[3] < y[3]) )
end

mutable struct NDSep # A nested dissection separator, "s"

    ## Hierarchy & dofs data
	ptr::Array{Int64,1} 			# A pointer array to the end of each cluster within s at the current stage of the elimination
    hrch::IT                        # A hierarchy, i.e., a tree of integer, monotone left -> right
    start::Int64                    # Where does the separator start
    depth::Int64                    # Where the clustering stand in terms of depth

    ## Data
	Piv::Matrix{Float64}			# Ass: Diagonal
	Low::Matrix{Float64}            # Ans: Lower-part 
	Upp::Matrix{Float64}			# Asn: Upper-part

    ## Block information, basically pointers to the data above
    APiv::Matrix{StridedMatrix{Float64}} # APiv[i,j] = ith partition of s, jth partition of s - (self lvl, self sep, i)
    ALow::Matrix{StridedMatrix{Float64}} # ALow[n,i] = ith partition of s, nth neighbor
    AUpp::Matrix{StridedMatrix{Float64}} # AUpp[i,n] = ith partition of s, nth neighbor
    Nbr::Vector{Cluster}                 # Nbr[n]    = nth neighbor (lvl, sep, part)
    # When fully merged, size(APiv) == (1,1), size(ALow,2) == size(AUpp,1) == 1

    ## Range information, for easier backward/forward pass
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

function Base.iterate(iter::BinTreeIt, state=(0, 2^iter.maxLevel))
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

function factorize(A::SparseMatrixCSC{Float64,Int64}, maxLevel::Int64, verbose::Bool=false)
   
    N = size(A, 1)
    @assert size(A, 2) == N
    @assert maxLevel >= 1

    # Symmetric partition
    (subseps, hrch, dofs) = mnd(A+A', maxLevel, verbose=verbose);

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

    ## Symbolic computation: Compute sparsity pattern, allocate & assemble
	if(verbose) @printf("Assembling...\n") end
    for (lvl, sep) in BinTreeIt(maxLevel)
        h = Tree[lvl][sep].hrch
        s = length(dofs[lvl][sep])
        np = length(h.nit[end])
        # Get the original neighbors
        (I, t1, t2) = findnz(A[:,get_dofs(Tree[lvl][sep])])
        (t1, J, t2) = findnz(A[get_dofs(Tree[lvl][sep]),:]) # FIXME: probably slow since in CSC
        nbrs = sort(unique(vcat(I,J)))
        Nbr = Set(dofs2cl[nbrs])
        # Add the children
        if lvl > 1
            union!(Nbr, Tree[lvl-1][2*sep-1].Nbr)
            union!(Nbr, Tree[lvl-1][2*sep  ].Nbr)
        end
        # Remove too small
        filter!(x -> x[1] > lvl, Nbr)
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
        nsizes = Vector{Int64}(map(x -> Tree[x[1]][x[2]].ptr[x[3]+1]-Tree[x[1]][x[2]].ptr[x[3]], Nbr))
        n = sum(nsizes)
        Low = zeros(n, s)
        Upp = zeros(s, n)
        Piv = zeros(s, s)
        Tree[lvl][sep].Low = Low
        Tree[lvl][sep].Upp = Upp
        Tree[lvl][sep].Piv = Piv
        # Get subblocks
        ALow = Matrix{StridedMatrix{Float64}}(undef, (length(Nbr), np))
        AUpp = Matrix{StridedMatrix{Float64}}(undef, (np, length(Nbr)))
        APiv = Matrix{StridedMatrix{Float64}}(undef, (np, np))
        sptr = Tree[lvl][sep].ptr
        nptr = sizes_to_ptr(nsizes)
        # Fill diag
        for i in 1:np
            idofs = get_dofs(Tree[lvl][sep], i)
            iids  = sptr[i]:(sptr[i+1]-1)
            # Diagonal
            for j in 1:np
                jdofs = get_dofs(Tree[lvl][sep], j)
                jids  = sptr[j]:(sptr[j+1]-1)
                APiv[i,j] = view(Piv, iids, jids)
                APiv[i,j][:,:] = A[idofs,jdofs]
            end
            # Lower/Upper
            for n in 1:length(Nbr)
				(lvln,sepn,in) = Nbr[n]
                ndofs = get_dofs(Tree[lvln][sepn], in)
                nids  = nptr[n]:(nptr[n+1]-1)
                ALow[n,i] = view(Low, nids, iids);  
				ALow[n,i][:,:] = A[ndofs, idofs]
                AUpp[i,n] = view(Upp, iids, nids)
                AUpp[i,n][:,:] = A[idofs, ndofs]
            end
        end
        Tree[lvl][sep].ALow = ALow
        Tree[lvl][sep].AUpp = AUpp
        Tree[lvl][sep].APiv = APiv
    end

    ## Numerical factorization
    for lvl = 1:maxLevel

		if(verbose) @printf("Lvl: %d, eliminating\n", lvl) end
        
        ## Eliminate
        for sep = 1:length(dofs[lvl])
            s = Tree[lvl][sep]
            # Check it's fully merged
            @assert all(size(s.APiv) .== (1,1))
            # GETRF, factor pivot (no pivoting so far)
            Ass = s.APiv[1,1]
            lu!(Ass, Val(false)); # Ass -> L U
            # TRSM, solve panel in place
            for n = 1:length(s.Nbr)
                (ln,sn,in) = s.Nbr[n]
                push!(s.NbrRanges, get_dofs(Tree[ln][sn],in))
                # Lower
                Ans = s.ALow[n,1]
                BLAS.trsm!('R','U','N','N', 1.0, Ass, Ans) # Ans -> Ans U^-1
                # Upper
                Asn = s.AUpp[1,n]
                BLAS.trsm!('L','L','N','U', 1.0, Ass, Asn) # Asn -> L^-1 Asn
            end
            # GEMM, update schur complement
            for ni = 1:length(s.Nbr)
                for nj = 1:length(s.Nbr)
                    Ais = s.ALow[ni,1]
                    Asj = s.AUpp[1,nj]
                    (li,si,ii) = s.Nbr[ni]
                    (lj,sj,jj) = s.Nbr[nj]
                    if (li,si) == (lj,sj)
                        Aij = Tree[li][si].APiv[ii,jj]
                    else
                        @assert li != lj # Can't be any edge between different separators at a given level
                        if li < lj
                            n = Tree[li][si]
                            id = searchsortedfirst(n.Nbr, s.Nbr[nj], lt=lt_cl)
                            Aij = n.AUpp[ii,id]
                        else
                            n = Tree[lj][sj]
                            id = searchsortedfirst(n.Nbr, s.Nbr[ni], lt=lt_cl)
                            Aij = n.ALow[id,jj]
                        end
                    end
                    BLAS.gemm!('N', 'N', - 1.0, Ais, Asj, 1.0, Aij) # Aij -> Aij - Ais Asj
                end
            end
        end

        ## Merge, i.e., redefine all the views
        for lvl2 = lvl+1:maxLevel
            for sep2 = 1:length(dofs[lvl2])
                s = Tree[lvl2][sep2]
                ptr2 = ptr_merged(s)
                # Pivot
                for i = 1:length(ptr2)-1
                    for j = 1:length(ptr2)-1
                        s.APiv[i,j] = view(s.Piv, ptr2[i]:(ptr2[i+1]-1), ptr2[j]:(ptr2[j+1]-1))
                    end
                end
                s.APiv = s.APiv[1:length(ptr2)-1, 1:length(ptr2)-1]
                # Lower/Upper
                n_i1, n_b1, n_merged = 1, 1, 1
                while n_b1 <= length(s.Nbr)
                    # What are n's siblings ?
                    (l_, s_, p_) = s.Nbr[n_b1]
                    (n_siblg, parent) = siblings(Tree[l_][s_], p_)
                    # All the siblings should follow
                    n_i2, n_b2 = n_i1, n_b1
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
                        s.AUpp[c, n_merged] = view(s.Upp, ptr2[c]:(ptr2[c+1]-1), n_i1:(n_i2-1))
                    end
                    # Next
                    n_merged += 1
                    n_i1, n_b1 = n_i2, n_b2
                end
                # Resize
                s.ALow = s.ALow[1:n_merged-1, 1:length(ptr2)-1]
                s.AUpp = s.AUpp[1:length(ptr2)-1, 1:n_merged-1]
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
            Lss = tree.t[l][s].APiv[1,1]
            p   = get_dofs(tree.t[l][s])
            # Pivot
            BLAS.trsv!('L', 'N', 'U', Lss, view(y, p)) # y[p] <- Lss^-1 y[p] 
            # Panel
            for (Lnp, n) in zip(tree.t[l][s].ALow[:,1], tree.t[l][s].NbrRanges)
                y[n] -= Lnp * y[p]
            end
        end
    end
    # Backward
    for l = maxLevel:-1:1
        for s = 1:2^(maxLevel-l)
            Uss = tree.t[l][s].APiv[1,1]
            p   = get_dofs(tree.t[l][s])
            # Panel
            for (Upn, n) in zip(tree.t[l][s].AUpp[1,:], tree.t[l][s].NbrRanges)
                y[p] -= Upn * y[n]
            end
            # Pivot
            BLAS.trsv!('U', 'N', 'N', Uss, view(y, p)) # y[p] <- Uss^-1 y[p] 
        end
    end
    return y[invperm(tree.p)]
end

