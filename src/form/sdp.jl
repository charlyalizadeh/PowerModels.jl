### sdp relaxations in the rectangular W-space
import LinearAlgebra: Hermitian, cholesky, Symmetric, diag, I
import SparseArrays: SparseMatrixCSC, sparse, spdiagm, findnz, spzeros, nonzeros


struct _SDconstraintDecomposition
    "Each sub-vector consists of bus IDs corresponding to a clique grouping"
    decomp::Vector{Vector{Int}}
    "`lookup_index[bus_id] --> idx` for mapping between 1:n and bus indices"
    lookup_index::Dict
    "A chordal extension and maximal cliques are uniquely determined by a graph ordering"
    ordering::Vector{Int}
end
import Base: ==
function ==(d1::_SDconstraintDecomposition, d2::_SDconstraintDecomposition)
    eq = true
    for f in fieldnames(_SDconstraintDecmposition)
        eq = eq && (getfield(d1, f) == getfield(d2, f))
    end
    return eq
end


function variable_bus_voltage(pm::AbstractSparseSDPWRMModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true,
                              extension_alg::AbstractChordalExtension=CholeskyExtension())

    if haskey(pm.ext, :SDconstraintDecomposition)
        decomp = pm.ext[:SDconstraintDecomposition]
        groups = decomp.decomp
        lookup_index = decomp.lookup_index
        lookup_bus_index = Dict((reverse(p) for p = pairs(lookup_index)))
    else
        cadj, lookup_index, ordering = _chordal_extension(pm, nw, extension_alg)
        groups = _maximal_cliques(cadj)
        lookup_bus_index = Dict((reverse(p) for p = pairs(lookup_index)))
        groups = [[lookup_bus_index[gi] for gi in g] for g in groups]
        pm.ext[:SDconstraintDecomposition] = _SDconstraintDecomposition(groups, lookup_index, ordering)
    end

    voltage_product_groups =
        var(pm, nw)[:voltage_product_groups] =
        Vector{Dict{Symbol, Array{JuMP.VariableRef,2}}}(undef, length(groups))

    for (gidx, group) in enumerate(groups)
        n = length(group)
        wr_start = zeros(n, n) + I
        voltage_product_groups[gidx] = Dict()
        WR = voltage_product_groups[gidx][:WR] =
            var(pm, nw)[:voltage_product_groups][gidx][:WR] =
            JuMP.@variable(pm.model, [i=1:n, j=1:n], Symmetric,
                base_name="$(nw)_$(gidx)_WR", start=wr_start[i,j])
        if report
            sol(pm, nw, :w_group, gidx)[:WR] = WR
        end

        WI = voltage_product_groups[gidx][:WI] =
            var(pm, nw)[:voltage_product_groups][gidx][:WI] =
            JuMP.@variable(pm.model, [1:n, 1:n],
                base_name="$(nw)_$(gidx)_WI", start=0.0)
        if report
            sol(pm, nw, :w_group, gidx)[:WI] = WI
        end
    end

    # voltage product bounds
    visited_buses = []
    visited_buspairs = []
    var(pm, nw)[:w] = Dict{Int,Any}()
    var(pm, nw)[:wr] = Dict{Tuple{Int,Int},Any}()
    var(pm, nw)[:wi] = Dict{Tuple{Int,Int},Any}()
    wr_min, wr_max, wi_min, wi_max = ref_calc_voltage_product_bounds(ref(pm, nw, :buspairs))
    for (gidx, voltage_product_group) in enumerate(voltage_product_groups)
        WR, WI = voltage_product_group[:WR], voltage_product_group[:WI]
        group = groups[gidx]
        ng = length(group)

        # diagonal bounds
        for (group_idx, bus_id) in enumerate(group)
            # group_idx indexes into group
            # bus_id indexes into ref(pm, nw, :bus)
            bus = ref(pm, nw, :bus, bus_id)

            wr_ii = WR[group_idx, group_idx]

            if bounded
                JuMP.set_upper_bound(wr_ii, (bus["vmax"])^2)
                JuMP.set_lower_bound(wr_ii, (bus["vmin"])^2)
            else
                JuMP.set_lower_bound(wr_ii, 0)
            end

            # for non-semidefinite constraints
            if !(bus_id in visited_buses)
                push!(visited_buses, bus_id)
                var(pm, nw, :w)[bus_id] = wr_ii
            end
        end

        # off-diagonal bounds
        offdiag_indices = [(i, j) for i in 1:ng, j in 1:ng if i != j]
        for (i, j) in offdiag_indices
            i_bus, j_bus = group[i], group[j]
            if (i_bus, j_bus) in ids(pm, nw, :buspairs)
                if bounded
                    JuMP.set_upper_bound(WR[i, j], wr_max[i_bus, j_bus])
                    JuMP.set_lower_bound(WR[i, j], wr_min[i_bus, j_bus])

                    JuMP.set_upper_bound(WI[i, j], wi_max[i_bus, j_bus])
                    JuMP.set_lower_bound(WI[i, j], wi_min[i_bus, j_bus])
                end

                # for non-semidefinite constraints
                if !((i_bus, j_bus) in visited_buspairs)
                    push!(visited_buspairs, (i_bus, j_bus))
                    var(pm, nw, :wr)[(i_bus, j_bus)] = WR[i, j]
                    var(pm, nw, :wi)[(i_bus, j_bus)] = WI[i, j]
                end
            end
        end
    end

    report && sol_component_value(pm, nw, :bus, :w, ids(pm, nw, :bus), var(pm, nw)[:w])
    report && sol_component_value_buspair(pm, nw, :buspairs, :wr, ids(pm, nw, :buspairs), var(pm, nw)[:wr])
    report && sol_component_value_buspair(pm, nw, :buspairs, :wi, ids(pm, nw, :buspairs), var(pm, nw)[:wi])
end


function constraint_model_voltage(pm::AbstractSparseSDPWRMModel, n::Int)
    _check_missing_keys(var(pm, n), [:voltage_product_groups], typeof(pm))

    pair_matrix(group) = [(i, j) for i in group, j in group]

    decomp = pm.ext[:SDconstraintDecomposition]
    groups = decomp.decomp
    voltage_product_groups = var(pm, n)[:voltage_product_groups]

    # semidefinite constraint for each group in clique grouping
    for (gidx, voltage_product_group) in enumerate(voltage_product_groups)
        _check_missing_keys(voltage_product_group, [:WR,:WI], typeof(pm))

        group = groups[gidx]
        ng = length(group)
        WR = voltage_product_group[:WR]
        WI = voltage_product_group[:WI]

        # Lower-dimensional SOC constraint equiv. to SDP for 2-vertex
        # clique
        if ng == 2
            wr_ii = WR[1, 1]
            wr_jj = WR[2, 2]
            wr_ij = WR[1, 2]
            wi_ij = WI[1, 2]
            wi_ji = WI[2, 1]

            # standard SOC form (Mosek doesn't like rotated form)
            JuMP.@constraint(pm.model, [(wr_ii + wr_jj), (wr_ii - wr_jj), 2*wr_ij, 2*wi_ij] in JuMP.SecondOrderCone())
            JuMP.@constraint(pm.model, wi_ij == -wi_ji)
        else
            JuMP.@constraint(pm.model, [WR WI; -WI WR] in JuMP.PSDCone())
        end
    end

    # linking constraints
    tree = _prim(_overlap_graph(groups))
    overlapping_pairs = [Tuple(CartesianIndices(tree)[i]) for i in (LinearIndices(tree))[findall(x->x!=0, tree)]]
    for (i, j) in overlapping_pairs
        gi, gj = groups[i], groups[j]
        var_i, var_j = voltage_product_groups[i], voltage_product_groups[j]

        Gi, Gj = pair_matrix(gi), pair_matrix(gj)
        overlap_i, overlap_j = _overlap_indices(Gi, Gj)
        indices = zip(overlap_i, overlap_j)
        for (idx_i, idx_j) in indices
            JuMP.@constraint(pm.model, var_i[:WR][idx_i] == var_j[:WR][idx_j])
            JuMP.@constraint(pm.model, var_i[:WI][idx_i] == var_j[:WI][idx_j])
        end
    end
end
