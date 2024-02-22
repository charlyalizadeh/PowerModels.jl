"""
    adj, lookup_index = _adjacency_matrix(pm, nw)
Return:
- a sparse adjacency matrix
- `lookup_index` s.t. `lookup_index[bus_id]` returns the integer index
of the bus with `bus_id` in the adjacency matrix.
"""
function _adjacency_matrix(pm::AbstractPowerModel, nw::Int=nw_id_default)
    bus_ids = ids(pm, nw, :bus)
    buspairs = ref(pm, nw, :buspairs)

    nb = length(bus_ids)
    nl = length(buspairs)

    lookup_index = Dict((bi, i) for (i, bi) in enumerate(bus_ids))
    f = [lookup_index[bp[1]] for bp in keys(buspairs)]
    t = [lookup_index[bp[2]] for bp in keys(buspairs)]

    return sparse([f;t], [t;f], ones(2nl), nb, nb), lookup_index
end


"""
    mc = _maximal_cliques(cadj, peo)
Given a chordal graph adjacency matrix and perfect elimination
ordering, return the set of maximal cliques.
"""
function _maximal_cliques(cadj::SparseMatrixCSC, peo::Vector{Int})
    nb = size(cadj, 1)

    # use peo to obtain one clique for each vertex
    cliques = Vector(undef, nb)
    for (i, v) in enumerate(peo)
        Nv = findall(x->x!=0, cadj[:, v])
        cliques[i] = union(v, intersect(Nv, peo[i+1:end]))
    end

    # now remove cliques that are strict subsets of other cliques
    mc = Vector()
    for c1 in cliques
        # declare clique maximal if it is a subset only of itself
        if sum([issubset(c1, c2) for c2 in cliques]) == 1
            push!(mc, c1)
        end
    end
    # sort node labels within each clique
    mc = [sort(c) for c in mc]
    return mc
end
_maximal_cliques(cadj::SparseMatrixCSC) = _maximal_cliques(cadj, _mcs(cadj))


"""
    peo = _mcs(A)
Maximum cardinality search for graph adjacency matrix A.
Returns a perfect elimination ordering for chordal graphs.
"""
function _mcs(A)
    n = size(A, 1)
    w = zeros(Int, n)
    peo = zeros(Int, n)
    unnumbered = collect(1:n)

    for i = n:-1:1
        z = unnumbered[argmax(w[unnumbered])]
        filter!(x -> x != z, unnumbered)
        peo[i] = z

        Nz = findall(x->x!=0, A[:, z])
        for y in intersect(Nz, unnumbered)
            w[y] += 1
        end
    end
    return peo
end


"""
    T = _prim(A, minweight=false)
Return minimum spanning tree adjacency matrix, given adjacency matrix.
If minweight == false, return the *maximum* weight spanning tree.

Convention: start with node 1.
"""
function _prim(A, minweight=false)
    n = size(A, 1)
    candidate_edges = []
    unvisited = collect(1:n)
    next_node = 1 # convention
    T = spzeros(Int, n, n)

    while length(unvisited) > 1
        current_node = next_node
        filter!(node -> node != current_node, unvisited)

        neighbors = intersect(findall(x->x!=0, A[:, current_node]), unvisited)
        current_node_edges = [(current_node, i) for i in neighbors]
        append!(candidate_edges, current_node_edges)
        filter!(edge -> length(intersect(edge, unvisited)) == 1, candidate_edges)
        weights = [A[edge...] for edge in candidate_edges]
        next_edge = minweight ? candidate_edges[indmin(weights)] : candidate_edges[argmax(weights)]
        filter!(edge -> edge != next_edge, candidate_edges)
        weight = minweight ? minimum(weights) : maximum(weights)
        T[next_edge[1], next_edge[2]] = weight
        T[next_edge[2], next_edge[1]] = weight
        next_node = intersect(next_edge, unvisited)[1]
    end
    return T
end


"""
    A = _overlap_graph(groups)
Return adjacency matrix for overlap graph associated with `groups`.
I.e. if `A[i, j] = k`, then `groups[i]` and `groups[j]` share `k` elements.
"""
function _overlap_graph(groups)
    n = length(groups)
    I = Vector{Int}()
    J = Vector{Int}()
    V = Vector{Int}()
    for (i, gi) in enumerate(groups)
        for (j, gj) in enumerate(groups)
            if gi != gj
                overlap = length(intersect(gi, gj))
                if overlap > 0
                    push!(I, i)
                    push!(J, j)
                    push!(V, overlap)
                end
            end
        end
    end
    return sparse(I, J, V, n, n)
end

"""
    n = _neighbors(adj, v; exclude)

Return the neighbors of the vertex `v` in `adj` (undirected).
"""
function _neighbors(adj::SparseMatrixCSC, v::Int; exclude=[v])
	return [i for i in 1:adj.n if (adj[v, i] != 0 || adj[i, v] != 0) && i ∉ exclude]
end


"""
    is_complete = _check_neighboor_complete(adj, v, v_sub)

Return `true` if the neighbors of `v` in the subgraph of `adj` formed by the vertices in `v_sub` is a complete graph.
"""
function _check_neighboor_complete(adj::SparseMatrixCSC, v::Int, v_sub::Vector{Int}=1:adj.n)
    neighbors = filter(n -> n in v_sub, _neighbors(adj, v))
    for i in 1:length(neighbors)-1
        for j in i+1:length(neighbors)
            if adj[neighbors[i], neighbors[j]] == 0
                return false
            end
        end
    end
    return true
end


"""
    is_chordal = _check_chordal(adj)

Return `true` if the adjacency matrix `adj` corresponds to a chordal graph.
"""
function _check_chordal(adj::SparseMatrixCSC)
    peo = _mcs(adj)
    for i in 1:length(peo)-1
        if !_check_neighboor_complete(adj, peo[i], peo[i:end])
            return false
        end
    end
    return true
end


"""
    _make_neighborhood_complete!(adj, v)

Make the neighborhood of `v` complete in `adj`.
"""
function _make_neighborhood_complete!(adj::SparseMatrixCSC, v::Int)
    neighbors = _neighbors(adj, v)
    for i in 1:length(neighbors) - 1
        for j in i+1:length(neighbors)
            adj[neighbors[i], neighbors[j]] = 1.0
            adj[neighbors[j], neighbors[i]] = 1.0
        end
    end
end

"""
    pmdo = _perfect_minimum_degree_ordering(adj)

Return a perfect minimum degree ordering of the graph represented by `adj`.
"""
function _perfect_minimum_degree_ordering(adj::SparseMatrixCSC)
    vertex_map = 1:adj.n
    perm = Int[]
    for i in 1:adj.n
        vmindeg = argmin([sum(adj[1:end, j]) for j in 1:adj.n])
        push!(perm, vertex_map[vmindeg])
        _make_neighborhood_complete!(adj, vmindeg)
        vertex_map = [vertex_map[vmindeg] <= vertex_map[j] ? vertex_map[j + 1] : vertex_map[j] for j in 1:adj.n - 1]
        adj = adj[1:end .!= vmindeg, 1:end .!= vmindeg]
    end
    return perm
end


"""
    pmdo = _perfect_minimum_degree_ordering(adj)

Return a perfect minimum degree ordering of the graph build from the power model `pm`
"""
function _perfect_minimum_degree_ordering(pm::AbstractSparseSDPWRMModel, nw::Int=nw_id_default)
    adj, lookup_index = _adjacency_matrix(pm, nw)
    return _perfect_minimum_degree_ordering(adj), lookup_index
end

"""
    _make_subgraph_complete!(adj, vertices)

Make the subgraph of `adj` composed of the vertices in `vertices` complete.
"""
function _make_subgraph_complete!(adj::SparseMatrixCSC, vertices::Vector{Int})
    for i in 1:length(vertices)-1
        for j in i+1:length(vertices)
            adj[vertices[i], vertices[j]] = 1.0
            adj[vertices[j], vertices[i]] = 1.0
        end
    end
end
