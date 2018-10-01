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
                BLAS.gemv!('N', -1.0, Lnp, view(y, p), 1.0, view(y, n)) # y[n] <- y[n] - Lnp * y[p]
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
                BLAS.gemv!('N', -1.0, Upn, view(y, n), 1.0, view(y, p)) # y[p] <- y[p] - Upn * y[n]
            end
            # Pivot
            BLAS.trsv!('U', 'N', 'N', Uss, view(y, p)) # y[p] <- Uss^-1 y[p] 
        end
    end
    return y[invperm(tree.p)]
end

