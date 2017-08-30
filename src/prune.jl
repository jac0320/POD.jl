"""
    Algorithm preparation with info.
    THe partition information every iteration.
"""
function build_graph(m::PODNonlinearModel)

    node_graph = Dict()
    for i in m.var_discretization_mip
        node_graph[i] = [true for j in 1:(length(m.discretization[i])-1)]
    end

    return node_graph
end

function arc_consistency(m::PODNonlinearModel)

end

function process_init_graph(m::PODNonlinearModel)

    init_g = build_graph(m)
    create_bounding_mip(m)
    m.sense_orig == :Min ? @constraint(m.model_mip, m.model_mip.obj >= m.best_obj) : @constraint(m.model_mip, m.model_mip.obj <= m.best_obj)
    @objective(m.model_mip, Min, 0.0)
    for var in keys(init_g)
        for i in 1:(length(m.discretization[var])-1)
            println("Limiting VAR $(var) with bounds $([m.discretization[var][i], m.discretization[var][i+1]])")
            kept_bounds = switch_bound(m, var, [m.discretization[var][i], m.discretization[var][i+1]])
            # print(m.model_mip)
            status = solve(m.model_mip)
            if status == :Infeasible
                init_g[var][i] = false
            else
                tighter_relax_obj = getobjectivevalue(m.model_mip)
                (tighter_relax_obj > m.best_obj) && println("PRUNING!!! $(tighter_relax_obj)")
            end
            restore_bound(m, var, kept_bounds)
        end
        println(init_g[var])
    end
    # error("STOP")
    return
end

function switch_bound(m::PODNonlinearModel, var_ind::Int, bound::Vector)

    @assert bound[1] <= bound[2]
    kept_bounds = [m.l_var_tight[var_ind], m.u_var_tight[var_ind]]

    setlowerbound(Variable(m.model_mip, var_ind), bound[1])
    setupperbound(Variable(m.model_mip, var_ind), bound[2])

    return kept_bounds
end

function restore_bound(m::PODNonlinearModel, var_ind::Int, kept_bound::Vector=[])

    isempty(kept_bound) && (kept_bound=[m.l_var_tight[var_ind], m.u_var_tight[var_ind]])

    @assert kept_bound[1] <= kept_bound[2]
    setlowerbound(Variable(m.model_mip, var_ind), kept_bound[1])
    setupperbound(Variable(m.model_mip, var_ind), kept_bound[2])

    return
end

"""
    Eliminate partitions by adding constraints
"""
function prune_partitions(m::PODNonlinearModel, var_ind::Int;bounds::Vector=[], )

    return
end
