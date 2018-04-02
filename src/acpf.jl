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

    @assert 2*N+2*G+4*E == m.num_var_orig

    return acpf
end

function acpf_algo_measurements(m::PODNonlinearModel; sol=nothing)

    sol == nothing ? sol = m.best_bound_sol : sol = sol

    println("@@@@@@@@@@@@@@@@ ACPF @@@@@@@@@@@@@@@@@@")
    # measure_relaxed_deviation(m, sol=sol)           # Experimental code
    # measure_binding_constraints(m, sol=sol)         # Experimental code
    # acpf_get_v_val(m, sol)
    # acpf_get_g_val(m, sol)
    # acpf_get_flow_val(m, sol)
    acpf_focus_on_congested_line(m, sol)

    println("@@@@@@@@@@@@@@@@ ACPF @@@@@@@@@@@@@@@@@@")

    return
end

function acpf_algo_position(m::PODNonlinearModel)

    adjust_branch_priority(m, acpf_build_priority(m))

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

function acpf_focus_on_congested_line(m::PODNonlinearModel, sol::Vector)

    branch_idxs = sort([i for i in keys(m.extension[:p])], by=x->x[1])
    congested_cnt = 0

    congestions = []
    for i in branch_idxs
        p = m.extension[:p][i]
        q = m.extension[:q][i]
        rate_a = m.user_parameters.data["branch"][string(i[1])]["rate_a"]
        push!(congestions, (i, abs(sol[p]^2+sol[q]^2-rate_a^2)))
        if abs(sol[p]^2+sol[q]^2-rate_a^2) < 1e-3
            congested_cnt += 1
        end
        println("FLOW $(i[1]) : $(i[2]) -> $(i[3]) P=$(sol[p]) | Q=$(sol[q]) | RA-DIST $(abs(sol[p]^2+sol[q]^2-rate_a^2))")
    end
    sort!(congestions, by=x->x[2])
    println("--- TOTAL congested lines (1e-3) = $(congested_cnt)")
    # m.disc_vars = acpf_reselect_disc_vars(m, congestions)
    return
end

function acpf_reselect_disc_vars(m::PODNonlinearModel, congestions::Any)

    disc_vars = Set()
    println("--- Most congested lines $(congestions[1][1]): DIST=$(congestions[1][2])")
    for v in acpf_get_v_idxs(m, congestions[1][1])
        push!(disc_vars, v)
    end

    for b in congestions
        if b[2] <= 1e-3
            vars = acpf_get_v_idxs(m, b[1])
            for v in vars
                push!(disc_vars, v)
            end
        end
    end

    println("--- Reselecting $(length(disc_vars)) Variables")

    return [v for v in disc_vars]
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
