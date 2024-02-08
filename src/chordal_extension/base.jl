import LinearAlgebra: Hermitian, cholesky, Symmetric, diag, I
import SparseArrays: SparseMatrixCSC, sparse, spdiagm, findnz, spzeros, nonzeros

abstract type AbstractChordalExtension end
mutable struct CholeskyExtension <: AbstractChordalExtension
    shift::Float64
    perm::Union{Nothing, Vector{Int}}
    CholeskyExtension(shift = 0.0, perm = nothing) = new(shift, perm)
end
mutable struct MinimumDegreeExtension <: AbstractChordalExtension end


abstract type AbstractMerge end
mutable struct MolzahnMerge <: AbstractMerge
    L::Float64
    MolzahnMerge(L::Float64=0.1) = new(L)
end


function _filter_flipped_pairs!(pairs)
    for (i, j) in pairs
        if i != j && (j, i) in pairs
            filter!(x -> x != (j, i), pairs)
        end
    end
end


"""
    idx_a, idx_b = _overlap_indices(A, B)
Given two arrays (sizes need not match) that share some values, return:

- linear index of shared values in A
- linear index of shared values in B

Thus, A[idx_a] == B[idx_b].
"""
function _overlap_indices(A::Array, B::Array, symmetric=true)
    overlap = intersect(A, B)
    symmetric && _filter_flipped_pairs!(overlap)
    idx_a = [something(findfirst(isequal(o), A), 0) for o in overlap]
    idx_b = [something(findfirst(isequal(o), B), 0) for o in overlap]
    return idx_a, idx_b
end


"""
    ps = _problem_size(groups)
Returns the sum of variables and linking constraints corresponding to the
semidefinite constraint decomposition given by `groups`. This function is
not necessary for the operation of clique merge, since `merge_cost`
computes the change in problem size for a proposed group merge.
"""
function _problem_size(groups)
    nvars(n::Integer) = n*(2*n + 1)
    A = _prim(_overlap_graph(groups))
    return sum(nvars.(Int.(nonzeros(A)))) + sum(nvars.(length.(groups)))
end


"""
    cadj, lookup_index, ordering = _chordal_extension(pm, nw, extension_alg)
Return:
- a sparse adjacency matrix corresponding to a chordal extension
of the power grid graph. The algorithm used is defined by `extension_alg`.
- `lookup_index` s.t. `lookup_index[bus_id]` returns the integer index
of the bus with `bus_id` in the adjacency matrix.
- the graph ordering that may be used to reconstruct the chordal extension
"""
function _chordal_extension(pm::AbstractPowerModel, nw::Int, extension_alg::AbstractChordalExtension)
    adj, lookup_index = _adjacency_matrix(pm, nw)
    cadj, p = _chordal_extension(adj, extension_alg)
    return cadj, lookup_index, p
end
