using PowerModels
import JuMP
import SCS
using SparseArrays
import SparseArrays: spdiagm, sparse, findnz
import LinearAlgebra: cholesky, Hermitian, diag
import InfrastructureModels
import PowerModels: parse_file, instantiate_model, build_opf, _adjacency_matrix, _chordal_extension, MolzahnMerge, _maximal_cliques, _overlap_graph, _prim_max, _merge_cliques!, _check_chordal
using Ipopt

using Graphs
using GraphPlot
using Compose
using ImageView
using FileIO
using Mosek
using MosekTools
using JSON

function ne(adj::SparseMatrixCSC)
    nedges = 0
    for i in 1:adj.n-1
        for j in i+1:adj.n
            nedges += (adj[i, j] != 0)
        end
    end
    return nedges
end

function edges_adj(adj::SparseMatrixCSC)
	edge_list = []
	for i in 1:adj.n-1
		for j in i+1:adj.n
			if adj[i, j] != 0
				push!(edge_list, [i, j])
			end
		end
	end
	return edge_list
end

function build_adj_from_cliques(cliques)
	nb_vertex = length(union(cliques...))
	adj = sparse(zeros(nb_vertex, nb_vertex))
	for c in cliques
		PowerModels._make_subgraph_complete!(adj, c)
	end
	return adj
end

function display_cliquegraph(cliques)
	g = SimpleGraph(length(cliques))
	for i in 1:length(cliques)-1
		for j in i+1:length(cliques)
			if length(intersect(cliques[i], cliques[j])) != 0
				add_edge!(g, i, j)
			end
		end
	end
	nodelabel = map(string, cliques)
	edgelabel = [length(intersect(cliques[e.src], cliques[e.dst])) for e in edges(g)]
	gimg = gplot(g; nodelabel=nodelabel, edgelabel=edgelabel)
	draw(PNG("cliquegraph.png", 16cm, 16cm), gimg)
end

function display_cliquetree(cliques, cliquetree)
	g = SimpleGraph(length(cliques))
	for edge in cliquetree
		add_edge!(g, edge[1], edge[2])
	end
	nodelabel = map(string, cliques)
	edgelabel = [length(intersect(cliques[e.src], cliques[e.dst])) for e in edges(g)]
	gimg = gplot(g; nodelabel=nodelabel, edgelabel=edgelabel)
	draw(PNG("cliquetree.png", 16cm, 16cm), gimg)
end

function display_graph(adj)
	g = SimpleGraph(adj.n)
	for i in 1:adj.n-1
		for j in i+1:adj.n
			if adj[i, j] != 0
				add_edge!(g, i, j)
			end
		end
	end
	nodelabel = 1:adj.n
	gimg = gplot(g; nodelabel=nodelabel)
	draw(PNG("graph.png", 16cm, 16cm), gimg)
end


function display_cliquetree_weight(cliques, cliquetree)
	weight = 0
	for edge in cliquetree
		weight += length(intersect(cliques[edge[1]], cliques[edge[2]]))
	end
	println("WEIGHT CLIQUETREE = $weight")
end

function test_merge()
    data = PowerModels.parse_file("data/MATPOWER/case145.m")
    #data = PowerModels.parse_file("data/MATPOWER/case1354pegase.m")
    pm = InfrastructureModels.InitializeInfrastructureModel(SparseSDPWRMPowerModel, data, PowerModels._pm_global_keys, PowerModels.pm_it_sym)
    PowerModels.ref_add_core!(pm.ref)
    nw = collect(nw_ids(pm))[1]

	adj, lookup_index = PowerModels._adjacency_matrix(pm)
	lookup_index_inverse = Dict(v => k for (k, v) in lookup_index)
	order = sort(1:adj.n, by=x->lookup_index_inverse[x])
	#println(order)
	display(sparse(adj))
	adj = adj[order, order]
    diag_el = sum(adj, dims=1)[:]
    W = Hermitian(-adj + spdiagm(0 => diag_el .+ 1))
	display(sparse(W))
	display_graph(adj)
	println("EDGES BEFORE CH: $(ne(adj))")


    cadj, lookup_index, perm = PowerModels._chordal_extension(pm, nw, PowerModels.CholeskyExtension())
    maximal_cliques = PowerModels._maximal_cliques(cadj)
	cadj_bis = build_adj_from_cliques(maximal_cliques)
    clique_tree = _prim_max(_overlap_graph(maximal_cliques))

	println("EDGES CH: $(ne(cadj))")
	println("Number of cliques CHOLESKY: $(length(maximal_cliques))")

	display_cliquegraph(maximal_cliques)
	display_cliquetree(maximal_cliques, edges_adj(clique_tree))
	display_cliquetree_weight(maximal_cliques, edges_adj(clique_tree))

    merge_alg = PowerModels.MolzahnMerge(0.5)
    PowerModels._merge_cliques!(cadj, maximal_cliques, clique_tree, merge_alg)
    clique_tree = _prim_max(_overlap_graph(maximal_cliques))
	cadj_bis = build_adj_from_cliques(maximal_cliques)

	println("EDGES MERGED: $(ne(cadj))")
	println("GRAPH MERGED CHORDAL: $(_check_chordal(cadj))")
	println("Number of cliques MERGE: $(length(maximal_cliques))")
	display_cliquetree_weight(maximal_cliques, edges_adj(clique_tree))
end

function test_solve()
    sdp_solver = Mosek.Optimizer
    data = PowerModels.parse_file("data/MATPOWER/case9.m")
    #data = PowerModels.parse_file("data/MATPOWER/case14.m")
    pm = InfrastructureModels.InitializeInfrastructureModel(SparseSDPWRMPowerModel, data, PowerModels._pm_global_keys, PowerModels.pm_it_sym)
    #pm_merged = InfrastructureModels.InitializeInfrastructureModel(SparseSDPWRMPowerModel, data, PowerModels._pm_global_keys, PowerModels.pm_it_sym)
    PowerModels.ref_add_core!(pm.ref)
    #PowerModels.ref_add_core!(pm_merged.ref)

    nw = collect(nw_ids(pm))[1]

    ## Chordal Extension
    cadj, lookup_index, perm = PowerModels._chordal_extension(pm, nw, PowerModels.CholeskyExtension())
    #println(lookup_index)

    cliques = PowerModels._maximal_cliques(cadj)
    lookup_bus_index = Dict((reverse(p) for p = pairs(lookup_index)))
    #println(lookup_bus_index)
    groups = [[lookup_bus_index[gi] for gi in g] for g in cliques]

    pm.ext[:SDconstraintDecomposition] = PowerModels._SDconstraintDecomposition(groups, lookup_index, perm)
    PowerModels.build_opf(pm)
	#println(pm.model)
	#print(pm.model)
    result = optimize_model!(pm, optimizer=sdp_solver)
	#println(json(result, 4))
    println("Time: $(result["solve_time"])\nObjective: $(result["objective"])")
end

test_solve()
#solve_ac_opf("data/MATPOWER/case9.m", Ipopt.Optimizer)
#test_solve()
