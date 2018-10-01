using Printf
using LinearAlgebra
using SparseArrays

function tensor_grid(grids)
    n = 1
    for g in grids
        n *= length(g)
    end
    X = zeros(length(grids), n)
    for (i,x) in enumerate(Base.product(grids...))
        X[:,i] = [x...]
    end
    return X
end

function show_array(arr, name=nothing)
    if name != nothing
        @printf "Array %s\n" name
    else
        @printf "Array\n"
    end
    show(IOContext(stdout, :compact => true), "text/plain", arr)
    @printf "\n"
end



# Test matrix generation
An = n -> spdiagm(-1 => -ones(n-1), 0 => 2*ones(n), 1 => -ones(n-1))

# order-2 Laplacian in 1, 2, or 3D
function Ad(n, d)
    I = sparse(UniformScaling{Float64}(1.0), n, n)
    if d == 1
        A = An(n)
    elseif d == 2
        A = kron(An(n), I) + kron(I, An(n))
    elseif d == 3
        A = kron(An(n), I, I) + kron(I, An(n), I) + kron(I, I, An(n))
    end
    return A
end

# order-2/4/6/8 Laplacian in 1, 2, or 3D
function Ad_o(n, d, o)
    order2 = ([1, -2, 1], [1, 0, 1])
    order4 = ([-1/12, 4/3, -5/2, 4/3, -1/12], [2, 1, 0, 1, 2])
    order6 = ([1/90, -3/20, 3/2, -49/18, 3/2, -3/20, 1/90], [3, 2, 1, 0, 1, 2, 3])
    order8 = ([-1/560, 8/315, -1/5, 8/5, -205/72, 8/5, -1/5, 8/315, -1/560], [4, 3, 2, 1, 0, 1, 2, 3, 4])
    e = ones(n)
    I = sparse(UniformScaling{Float64}(1.0), n, n)
    if o == 2
        Ao = spdiagm(([c*e[1:n-s] for (c,s) in zip(order2...)]...), (-1, 0, 1), n, n)
    elseif o == 4
        Ao = spdiagm(([c*e[1:n-s] for (c,s) in zip(order4...)]...), (-2, -1, 0, 1, 2), n, n)
    elseif o == 6
        Ao = spdiagm(([c*e[1:n-s] for (c,s) in zip(order6...)]...), (-3, -2, -1, 0, 1, 2, 3), n, n)
    elseif o == 8
        Ao = spdiagm(([c*e[1:n-s] for (c,s) in zip(order8...)]...), (-4, -3, -2, -1, 0, 1, 2, 3, 4), n, n)
    end
    if d == 1
        A = - Ao
    elseif d == 2
        A = kron(- Ao, I) + kron(I, - Ao)
    elseif d == 3
        A = kron(- Ao, I, I) + kron(I, - Ao, I) + kron(I, I, - Ao)
    end
    return A
end

function count_sorted(v::AbstractArray{V,1}) where {V}
    if length(v) == 0
        return 0
    end
    lastv = v[1]
    n = 1
    for i = 2:length(v)
        if v[i] != lastv
            n += 1
            lastv = v[i]
        end
    end
    return n
end

function sorted_to_ptr(v::AbstractArray{V}) where {V}
    l = count_sorted(v)
    r = Array{Int64,1}(undef, l+1)
    r[1] = 1
    r[l+1] = length(v)+1
    if length(v) == 0 || length(v) == 1
        return r
    end
    lastv = v[1]
    n = 1
    for i = 2:length(v)
        if v[i] != lastv
            r[n+1] = i
            n += 1
            lastv = v[i]
        end
    end
    return r
end

# [1, 2, 1] -> [1, 2, 4, 5]
function sizes_to_ptr(v::AbstractArray)
    if length(v) == 0
        return [1]
    else
        return vcat([1], 1 .+ cumsum(v))
    end
end

# Fill Adst with Asrc[n,s]
function fill_dense_block!(Asrc::SparseMatrixCSC{Tv,Ti}, Adst::StridedMatrix{Tv}, n::UnitRange{Ti}, s::UnitRange{Ti},) where {Tv, Ti <: Integer}
    for j = 1:length(s)
        k = A.colptr[s[j]]:(A.colptr[s[j]+1]-1)
        I = A.rowval[k]
        V = A.nzval[k]
        ids = searchsortedfirst(I, n[1]):searchsortedlast(I, n[end])
        for i in ids
            Adst[I[i] - n[1] + 1, j] = V[i]
        end
    end
end

# In a binary tree, returns true if (l1, s1) and (l2, s2) are on same path to the root
#			(3,1)
#	      /       \
#     (2,1)      (2,2)
#     /  \       /   \
# (1,1) (1,2) (1,3) (1,4)
function are_connected(l1, s1, l2, s2, maxLevel)
    @assert l1 > 0 && l1 <= maxLevel && s1 > 0 && s1 <= 2^(maxLevel-l1)
    @assert l2 > 0 && l2 <= maxLevel && s2 > 0 && s2 <= 2^(maxLevel-l2)
    if l1 == l2
        return s1 == s2
    end
	L = [l1, l2]
    S = [s1, s2]
    (lmax, imax) = findmax(L)
    (lmin, imin) = findmin(L)
    smax = S[imax]
    smin = S[imin]
    while lmin < lmax
        lmin += 1
        smin = div(smin+1,2)
    end
    smin == smax
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

