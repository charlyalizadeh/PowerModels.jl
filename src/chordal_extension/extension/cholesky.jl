"""
    cadj, p = _chordal_extension(adj, extension_alg)

Return a cholesky chordal extension of `adj`.
The permutation used by the cholesky factorisation is `extension_alg.perm`.
If `extension_alg.perm` is `nothing` use an approximate minimum degree permutation.
"""
function _chordal_extension(adj::SparseMatrixCSC, extension_alg::CholeskyExtension)
    nb = size(adj, 1)
    diag_el = sum(adj, dims=1)[:]
    W = Hermitian(-adj + spdiagm(0 => diag_el .+ 1))

    F = cholesky(W; shift = extension_alg.shift, perm = extension_alg.perm)
    L = sparse(F.L)
    p = F.p
    q = invperm(p)

    Rchol = L - spdiagm(0 => diag(L))
    f_idx, t_idx, V = findnz(Rchol)
    cadj = sparse([f_idx;t_idx], [t_idx;f_idx], ones(2*length(f_idx)), nb, nb)
    cadj = cadj[q, q] # revert to original bus ordering (invert cholfact permutation)
    return cadj, p
end


"""
    cadj, p = _chordal_extension(adj, extension_alg)

Return a cholesky chordale extension of `adj`.
The permutation used by the cholesky factorisation is a minimum degree ordering
"""
function _chordal_extension(adj::SparseMatrixCSC, extension_alg::MinimumDegreeExtension)
    perm = _perfect_minimum_degree_ordering(adj)
    alg = CholeskyExtension(0.0, perm)
    cadj, p = _chordal_extension(adj, alg)
    return cadj, p
end
