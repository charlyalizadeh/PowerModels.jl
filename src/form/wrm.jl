### sdp relaxations in the rectangular W-space
import LinearAlgebra: Hermitian, cholesky, Symmetric, diag, I
import SparseArrays: SparseMatrixCSC, sparse, spdiagm, findnz, spzeros, nonzeros


""
function constraint_current_limit_from(pm::AbstractWRMModel, n::Int, f_idx, c_rating_a)
    l,i,j = f_idx

    w_fr = var(pm, n, :w, i)

    p_fr = var(pm, n, :p, f_idx)
    q_fr = var(pm, n, :q, f_idx)
    JuMP.@constraint(pm.model, [w_fr*c_rating_a^2+1, 2*p_fr, 2*q_fr, w_fr*c_rating_a^2-1] in JuMP.SecondOrderCone())
end

""
function constraint_current_limit_to(pm::AbstractWRMModel, n::Int, t_idx, c_rating_a)
    l,j,i = t_idx

    w_to = var(pm, n, :w, j)

    p_to = var(pm, n, :p, t_idx)
    q_to = var(pm, n, :q, t_idx)
    JuMP.@constraint(pm.model, [w_to*c_rating_a^2+1, 2*p_to, 2*q_to, w_to*c_rating_a^2-1] in JuMP.SecondOrderCone())
end




""
function constraint_model_voltage(pm::AbstractWRMModel, n::Int)
    _check_missing_keys(var(pm, n), [:WR,:WI], typeof(pm))

    WR = var(pm, n)[:WR]
    WI = var(pm, n)[:WI]

    JuMP.@constraint(pm.model, [WR WI; -WI WR] in JuMP.PSDCone())
end


""
function variable_bus_voltage(pm::AbstractWRMModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    wr_min, wr_max, wi_min, wi_max = ref_calc_voltage_product_bounds(ref(pm, nw, :buspairs))
    bus_ids = ids(pm, nw, :bus)

    w_index = 1:length(bus_ids)
    lookup_w_index = Dict((bi,i) for (i,bi) in enumerate(bus_ids))

    WR_start = zeros(length(bus_ids), length(bus_ids)) + I

    WR = var(pm, nw)[:WR] = JuMP.@variable(pm.model,
        [i=1:length(bus_ids), j=1:length(bus_ids)], Symmetric, base_name="$(nw)_WR", start=WR_start[i,j]
    )
    if report
        sol(pm, nw)[:WR] = WR
    end

    WI = var(pm, nw)[:WI] = JuMP.@variable(pm.model,
        [1:length(bus_ids), 1:length(bus_ids)], base_name="$(nw)_WI", start=0.0
    )
    if report
        sol(pm, nw)[:WI] = WI
    end

    # bounds on diagonal
    for (i, bus) in ref(pm, nw, :bus)
        w_idx = lookup_w_index[i]
        wr_ii = WR[w_idx,w_idx]
        wi_ii = WR[w_idx,w_idx]

        if bounded
            JuMP.set_lower_bound(wr_ii, (bus["vmin"])^2)
            JuMP.set_upper_bound(wr_ii, (bus["vmax"])^2)

            #this breaks SCS on the 3 bus exmple
            #JuMP.set_lower_bound(wi_ii, 0)
            #JuMP.set_upper_bound(wi_ii, 0)
        else
            JuMP.set_lower_bound(wr_ii, 0)
        end
    end

    # bounds on off-diagonal
    for (i,j) in ids(pm, nw, :buspairs)
        wi_idx = lookup_w_index[i]
        wj_idx = lookup_w_index[j]

        if bounded
            JuMP.set_upper_bound(WR[wi_idx, wj_idx], wr_max[(i,j)])
            JuMP.set_lower_bound(WR[wi_idx, wj_idx], wr_min[(i,j)])

            JuMP.set_upper_bound(WI[wi_idx, wj_idx], wi_max[(i,j)])
            JuMP.set_lower_bound(WI[wi_idx, wj_idx], wi_min[(i,j)])
        end
    end

    var(pm, nw)[:w] = Dict{Int,Any}()
    for (i, bus) in ref(pm, nw, :bus)
        w_idx = lookup_w_index[i]
        var(pm, nw, :w)[i] = WR[w_idx,w_idx]
    end
    report && sol_component_value(pm, nw, :bus, :w, ids(pm, nw, :bus), var(pm, nw)[:w])

    var(pm, nw)[:wr] = Dict{Tuple{Int,Int},Any}()
    var(pm, nw)[:wi] = Dict{Tuple{Int,Int},Any}()
    for (i,j) in ids(pm, nw, :buspairs)
        w_fr_index = lookup_w_index[i]
        w_to_index = lookup_w_index[j]

        var(pm, nw, :wr)[(i,j)] = WR[w_fr_index, w_to_index]
        var(pm, nw, :wi)[(i,j)] = WI[w_fr_index, w_to_index]
    end
    report && sol_component_value_buspair(pm, nw, :buspairs, :wr, ids(pm, nw, :buspairs), var(pm, nw)[:wr])
    report && sol_component_value_buspair(pm, nw, :buspairs, :wi, ids(pm, nw, :buspairs), var(pm, nw)[:wi])
end
