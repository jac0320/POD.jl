function init_acpf_package(m::PODNonlinearModel, acpf=Dict())

    up = m.user_parameters

    N = acpf[:bus_count] = length(up.data["bus"])
    E = acpf[:branch_count] = length(up.data["branch"])
    G = acpf[:gen_count] = length(up.data["gen"])

    acpf[:vr] = Dict(i[1]=>up.var[:nw][0][:vr][i[1]].col for i in keys(up.var[:nw][0][:vr]))
    acpf[:vi] = Dict(i[1]=>up.var[:nw][0][:vi][i[1]].col for i in keys(up.var[:nw][0][:vi]))
    acpf[:pg] = Dict(i[1]=>up.var[:nw][0][:pg][i[1]].col for i in keys(up.var[:nw][0][:pg]))
    acpf[:qg] = Dict(i[1]=>up.var[:nw][0][:vr][i[1]].col for i in keys(up.var[:nw][0][:vr]))
    acpf[:p] = Dict(i[1]=>up.var[:nw][0][:p][i[1]].col for i in keys(up.var[:nw][0][:p]))
    acpf[:q] = Dict(i[1]=>up.var[:nw][0][:q][i[1]].col for i in keys(up.var[:nw][0][:q]))
    acpf[:warmstarter] = []

    @assert 2*N+2*G+4*E == m.num_var_orig

    return acpf
end

function acpf_algo_measurements(m::PODNonlinearModel; sol=nothing)

    sol == nothing ? sol = m.best_bound_sol : sol = sol

    # println("@@@@@@@@@@@@@@@@ ACPF @@@@@@@@@@@@@@@@@@")
    # measure_relaxed_deviation(m, sol=sol)           # Experimental code
    # # measure_binding_constraints(m, sol=sol)         # Experimental code
    # # acpf_get_v_val(m, sol)
    # # acpf_get_g_val(m, sol)
    # # acpf_get_flow_val(m, sol)
    # println("@@@@@@@@@@@@@@@@ ACPF @@@@@@@@@@@@@@@@@@")

    return
end

function acpf_pre_partition_construction(m::PODNonlinearModel; sol=nothing)
    println("@@@@@@@@@@@@@@@@ ACPF @@@@@@@@@@@@@@@@@@")
    sol == nothing ? sol = m.best_bound_sol : sol = sol
    measure_relaxed_deviation(m, sol=sol)           # Experimental cod
    acpf_disc_vars_heuristics(m, sol)
    println("@@@@@@@@@@@@@@@@ ACPF @@@@@@@@@@@@@@@@@@")
    return
end

function acpf_position_bounding_model(m::PODNonlinearModel)

    acpf_examine_partition_status(m)
    # adjust_branch_priority(m, acpf_build_priority(m))

    return
end

function acpf_examine_partition_status(m::PODNonlinearModel)

    acpf_initialize_partition_status(m)
    for k in keys(m.nonconvex_terms)
        va, vb = m.nonconvex_terms[k][:var_idxs]
        va_pcnt, vb_pcnt = length(m.discretization[va])-1, length(m.discretization[vb])-1
        y_idx = m.nonconvex_terms[k][:y_idx]
        y_part = (m.l_var_tight[y_idx], m.u_var_tight[y_idx])
        for i in 1:va_pcnt
            if m.disc_status[va][i]
                va_part = [m.discretization[va][i], m.discretization[va][i+1]]
                vb_part = [m.discretization[vb][1], m.discretization[vb][end]]
                pd_part = [minimum(va_part*vb_part'), maximum(va_part*vb_part')]
                if pd_part[2] < y_part[1] || pd_part[1] > y_part[2]
                    println("VAR$(va) PART $(i) (X)")
                    m.disc_status[va][i] = false
                    m.model_mip.colUpper[m.convhull_binary_links[va][i]] = 0.0
                end
            end
        end
        for i in 1:vb_pcnt
            if m.disc_status[vb][i]
                va_part = [m.discretization[va][1], m.discretization[va][end]]
                vb_part = [m.discretization[vb][i], m.discretization[vb][i+1]]
                pd_part = [minimum(va_part*vb_part'), maximum(va_part*vb_part')]
                if pd_part[2] < y_part[1] || pd_part[1] > y_part[2]
                    println("VAR$(vb) PART $(i) (X)")
                    m.disc_status[vb][i] = false
                    m.model_mip.colUpper[m.convhull_binary_links[vb][i]] = 0.0
                end
            end
        end
    end

    return
end

function acpf_initialize_partition_status(m::PODNonlinearModel)

    m.disc_status = Dict(i=>[true for j in 1:length(m.discretization[i])-1] for i in m.candidate_disc_vars)

    return
end


function acpf_relaxation_heuristic(m::PODNonlinearModel)

    m.extension[:warmstarter] = []
    mip_solver_verbosity(m, 0)
    println("~~~~~~~~~~~~~~~~ NEW CODE ~~~~~~~~~~~~~~~~~")
    # Set up an bounding model based on current discretization (with 1 extra partition)
    create_bounding_mip(m)

    neg_vr_part = Dict()
    for k in keys(m.extension[:vr])
        i = m.extension[:vr][k]
        neg_vr_part[i] = find_local_partition_idx(m.discretization[i], -m.best_bound_sol[i])
    end

    for k in keys(m.extension[:vr])
        i = m.extension[:vr][k]
        setlowerbound(Variable(m.model_mip, m.convhull_binary_links[i][neg_vr_part[i]]), 1.0)
    end

    heuristic_status = solve(m.model_mip)
    if heuristic_status in [:Optimal, :UserObjLimit, :UserLimit, :Suboptimal]
        m.loglevel > 99 && println("Found feasible bounding solution OBJ=$(m.model_mip.objVal)")
        m.extension[:warmstarter] = [round.(getvalue(Variable(m.model_mip, i)), 6) for i in 1:(m.num_var_orig+m.num_var_linear_mip+m.num_var_nonlinear_mip)]
    end
    println("~~~~~~~~~~~~~~~~ NEW CODE ~~~~~~~~~~~~~~~~~")

    mip_solver_verbosity(m, 1)

    return
end

function acpf_warmstart_bounding_model(m::PODNonlinearModel)

    
    if !isempty(m.extension[:warmstarter])
        m.model_mip.colVal = copy(m.extension[:warmstarter])
    end

    return
end

function acpf_build_priority(m::PODNonlinearModel)

    priority = zeros(Int, len(m.model_mip.colVal))

    for (cnt, v) in enumerate(keys(m.convhull_binary_links))
        for i in m.convhull_binary_links[v]
            priority[i] = cnt
        end
    end

    return priority
end

function acpf_get_v_val(m::PODNonlinearModel, sol::Vector)

    bus_idxs = sort([i for i in keys(m.extension[:vr])])
    for i in bus_idxs
        vr = m.extension[:vr][i]
        vi = m.extension[:vi][i]
        println("BUS$(i) VR=$(sol[vr]) [$(m.l_var_tight[vr]),$(m.u_var_tight[vr])] | VI=$(sol[vi]) [$(m.l_var_tight[vi]), $(m.u_var_tight[vi])]")
    end

    return
end

function acpf_get_g_val(m::PODNonlinearModel, sol::Vector)

    gen_idxs = sort([i for i in keys(m.extension[:pg])])
    for i in gen_idxs
        pg = m.extension[:pg][i]
        qg = m.extension[:qg][i]
        println("GEN$(i) VR=$(sol[pg]) [$(m.l_var_tight[pg]),$(m.u_var_tight[pg])] | VI=$(sol[qg]) [$(m.l_var_tight[qg]), $(m.u_var_tight[qg])]")
    end

    return
end

function acpf_get_flow_val(m::PODNonlinearModel, sol::Vector)

    branch_idxs = sort([i for i in keys(m.extension[:p])], by=x->x[1])
    binding_cnt = 0
    for i in branch_idxs
        p = m.extension[:p][i]
        q = m.extension[:q][i]
        rate_a = m.user_parameters.data["branch"][string(i[1])]["rate_a"]
        if abs(sol[p]^2+sol[q]^2-rate_a^2) <= 1e-1
            binding_cnt += 1
        end
        println("BINDING FLOW$(i[1]) $(i[2])->$(i[3]) P=$(sol[p]) | Q=$(sol[q]) | P²+Q²=$(sol[p]^2+sol[q]^2) | RATE-A² = $(rate_a^2) | RA-DIST $(abs(sol[p]^2+sol[q]^2-rate_a^2))")
    end
    println("BINDING TOTAL $(binding_cnt)")
    return
end

function acpf_disc_vars_heuristics(m::PODNonlinearModel, sol::Vector)

    branch_idxs = sort([i for i in keys(m.extension[:p])], by=x->x[1])
    congested_cnt = 0

    congestions = []
    for i in branch_idxs
        p = m.extension[:p][i]
        q = m.extension[:q][i]
        rate_a = m.user_parameters.data["branch"][string(i[1])]["rate_a"]
        if abs(sol[p]^2+sol[q]^2-rate_a^2) < 1e-3
            congested_cnt += 1
        end
        bind_diff = abs(sol[p]^2+sol[q]^2-rate_a^2)
        y_diff = 0.0
        for j in acpf_search_v_terms(m, acpf_convhull_info_index(m), acpf_get_v_idxs(m, i), i)
            y_diff += m.convhull_sol_info[j][2]
        end
        println("FLOW $(i[1]) : $(i[2]) -> $(i[3]) | RA-DIST $(bind_diff) | TOTAL-Y-DIFF $(y_diff) ")
        push!(congestions, (i, bind_diff, y_diff))
    end
    sort!(congestions, by=x->x[3], rev=true)
    println("--- TOTAL congested lines (1e-3) = $(congested_cnt)")

    m.disc_vars = acpf_reselect_disc_vars(m, congestions, congested_cnt >= 0.1 * length(branch_idxs))
    return
end

function acpf_reselect_disc_vars(m::PODNonlinearModel, congestions::Any, consider_congested_line=true)

    disc_vars = Set()

    println("--- Focused lines $(congestions[1][1]): R-DIST=$(congestions[1][2]) Y-DIST=$(congestions[1][3])")
    for i in 1:2
        for v in acpf_get_v_idxs(m, congestions[i][1])
            println("--- LINE $(congestions[i][1][1]):$(congestions[i][1][2])->$(congestions[i][1][3]) | VAR $(congestions[i])")
            push!(disc_vars, v)
        end
    end

    if consider_congested_line
        for b in congestions
            if b[2] <= 1e-4
                vars = acpf_get_v_idxs(m, b[1])
                for v in vars
                    push!(disc_vars, v)
                end
            end
        end
    end

    println("--- Reselecting $(length(disc_vars)) Variables")

    return [v for v in disc_vars]
end

function acpf_search_v_terms(m::PODNonlinearModel, info_dict::Dict, v_vars::Vector, branch::Tuple)

    vr_from_idx = m.user_parameters.var[:nw][0][:vr][branch[2]].col
    vi_from_idx = m.user_parameters.var[:nw][0][:vi][branch[2]].col
    vr_to_idx = m.user_parameters.var[:nw][0][:vr][branch[3]].col
    vi_to_idx = m.user_parameters.var[:nw][0][:vi][branch[3]].col

    entry_idx = [info_dict[Set(vr_from_idx)], info_dict[Set(vi_from_idx)]]
    push!(entry_idx, info_dict[Set(v for v in [vr_from_idx, vr_to_idx])])
    push!(entry_idx, info_dict[Set(v for v in [vi_from_idx, vi_to_idx])])
    push!(entry_idx, info_dict[Set(v for v in [vr_to_idx, vi_from_idx])])
    push!(entry_idx, info_dict[Set(v for v in [vr_from_idx, vi_to_idx])])

    return entry_idx
end

function acpf_get_v_idxs(m::PODNonlinearModel, branch::Tuple)

    return [m.user_parameters.var[:nw][0][:vr][branch[2]].col,
            m.user_parameters.var[:nw][0][:vr][branch[3]].col,
            m.user_parameters.var[:nw][0][:vi][branch[2]].col,
            m.user_parameters.var[:nw][0][:vi][branch[3]].col;]

end

function acpf_set_branching_priority(m::PODNonlinearModel, congestions::Any)

    return
end

function acpf_convhull_info_index(m::PODNonlinearModel)
    return index_info = Dict(i[5]=>cnt for (cnt, i) in enumerate(m.convhull_sol_info))
end
