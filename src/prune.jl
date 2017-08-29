"""
    Algorithm preparation with info.
    THe partition information every iteration.
"""
function build_graph(m::PODNonlinearModel)

    graph = Dict()


    return graph
end




"""
    Main algorithm
"""
function check_consistency()

    return
end

function switch_bound(m::PODNonlinearModel, var_ind::Int, bound::Vector)

    @assert bound[1] <= bound[2]
    kept_bounds = [m.l_var_tight[var_ind], m.u_var_tight[var_ind]]

    setlowerbound(Variable(m.model_mip, var_ind), bound[1])
    setupperbound(Variable(m.model_mip, var_ind), bound[2])

    return kept_bound
end

function restore_bound(m::PODNonlinearModel, var_ind::Int; kept_bound::Vector=[])

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
