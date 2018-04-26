function initialize_solution_pool(m::PODNonlinearModel, cnt::Int)

    s = Dict()

    s[:cnt] = cnt

    # Column dimension changing variable
    s[:len] = m.num_var_orig+m.num_var_linear_mip+m.num_var_nonlinear_mip

    # !! Be careful with the :vars when utilizing the dynamic discretization variable selection !!
    s[:vars] = [i for i in m.candidate_disc_vars if length(m.discretization[i]) > 2]

    s[:sol] = Vector{Vector}(cnt)                   # Solution value
    s[:obj] = Vector{Float64}(cnt)                  # Objecitve value
    s[:disc] = Vector{Dict}(cnt)                    # Discretization
    s[:stat] = [:Alive for i in 1:cnt]              # Solution status
    s[:iter] = [m.logs[:n_iter] for i in 1:cnt]     # Iteration collected
    s[:ubstart] = [false for i in 1:cnt]            # Solution used for ub multistart

    return s
end

"""
    Collect LB solutions
    Don't test this function
"""
function collect_lb_pool(m::PODNonlinearModel)

    # Always stick to the structural .discretization for algorithm consideration info
    # If in need, the scheme need to be refreshed with customized discretization info

    if m.mip_solver_id != "Gurobi" || m.obj_structure == :convex || isempty([i for i in m.model_mip.colCat if i in [:Int, :Bin]])
        warn("Skipping collecting solution pool procedure", once=true) # Only feaible with Gurobi solver
        return
    end

    # Construct a new solution pool for just this new iteration
    s = initialize_solution_pool(m, Gurobi.get_intattr(m.model_mip.internalModel.inner, "SolCount"))

    # Collect Solution and corresponding objective values
    for i in 1:s[:cnt]
        if m.mip_solver_id == "Gurobi"
            Gurobi.set_int_param!(m.model_mip.internalModel.inner, "SolutionNumber", i-1)
            s[:sol][i] = Gurobi.get_dblattrarray(m.model_mip.internalModel.inner, "Xn", 1, s[:len])
            s[:obj][i] = Gurobi.get_dblattr(m.model_mip.internalModel.inner, "PoolObjVal")
            s[:obj][i] = MathProgBase.eval_f(m.d_orig, s[:sol][i][1:m.num_var_orig])
        elseif m.mip_solver_id == "CPLEX"
            error("No implementation for CPLEX")
        end
        s[:disc][i] = Dict(j=>get_active_partition_idx(m.discretization, s[:sol][i][j],j) for j in s[:vars])
    end

    merge_solution_pool(m, s)

    return
end

"""
    Merge collected solution pools
"""
function merge_solution_pool(m::PODNonlinearModel, s::Dict)

    # Always stick to the structural .discretization for algorithm consideration info
    # If in need, the scheme need to be refreshed with customized discretization info

    # Update a few dimensional parameter
    var_idxs = s[:vars]

    lbv2p = Dict()  # Bookeeping the partition idx between iterations
    for v in var_idxs
        vpcnt = length(m.discretization[v]) - 1
        chosen_p = track_new_partition_idx(m.discretization, v, m.best_bound_sol[v])
        lbv2p[v] = [i for i in 1:vpcnt if !(i in chosen_p)]
    end

    for i in 1:m.bound_sol_pool[:cnt] # First update the existing solution pool
        # First update the discretization idx with the existing
        m.bound_sol_pool[:disc][i] = Dict(j => get_active_partition_idx(m.discretization, m.bound_sol_pool[:sol][i][j],j) for j in var_idxs)
        act = true # Then check if the update pool solution active partition idex is within the deactivated region
        for v in var_idxs
            (m.bound_sol_pool[:disc][i][v] in lbv2p[v]) || (act = false)
            act || (m.bound_sol_pool[:stat][i] = :Dead)  # No re-activation
            act || break
        end
    end

    for i in 1:s[:cnt] # Now perform the merge
        act = true # Then check if the update pool solution active partition index is within the deactivated region
        for v in var_idxs
            (s[:disc][i][v] in lbv2p[v]) || (act = false)
            act || (s[:stat][i] = :Dead)
            act || break
        end
        # Reject solutions that is around best bound to avoid traps
        if isapprox(s[:obj][i], m.best_bound;atol=m.tol)
            s[:stat][i] = :Dead
        end
        push!(m.bound_sol_pool[:sol], s[:sol][i])
        push!(m.bound_sol_pool[:obj], s[:obj][i])
        push!(m.bound_sol_pool[:disc], s[:disc][i])
        push!(m.bound_sol_pool[:stat], s[:stat][i])
        push!(m.bound_sol_pool[:iter], s[:iter][i])
        push!(m.bound_sol_pool[:ubstart], s[:ubstart][i])
    end

    # Update dimensional parameters
    m.bound_sol_pool[:cnt] = length(m.bound_sol_pool[:sol])
    m.bound_sol_pool[:vars] = var_idxs

    # Show the summary
    m.loglevel > 199 && println("POOL size = $(length([i for i in 1:m.bound_sol_pool[:cnt] if m.bound_sol_pool[:stat][i] != :Dead])) / $(m.bound_sol_pool[:cnt]) ")
    for i in 1:m.bound_sol_pool[:cnt]
        m.loglevel > 99 && m.bound_sol_pool[:stat][i] != :Dead && println("ITER $(m.bound_sol_pool[:iter][i]) | SOL $(i) | POOL solution obj = $(m.bound_sol_pool[:obj][i])")
    end

    return
end

"""
    check_solution_history(m::PODNonlinearModel, ind::Int)

Check if the solution is alwasy the same within the last disc_consecutive_forbid iterations. Return true if suolution in invariant.
"""
function check_solution_history(m::PODNonlinearModel, ind::Int)

    m.disc_consecutive_forbid == 0 && return false
    (m.logs[:n_iter] < m.disc_consecutive_forbid) && return false

    sol_val = m.bound_sol_history[mod(m.logs[:n_iter]-1, m.disc_consecutive_forbid)+1][ind]
    for i in 1:(m.disc_consecutive_forbid-1)
        search_pos = mod(m.logs[:n_iter]-1-i, m.disc_consecutive_forbid)+1
        !isapprox(sol_val, m.bound_sol_history[search_pos][ind]; atol=m.disc_rel_width_tol) && return false
    end

    m.loglevel > 99 && println("Consecutive bounding solution on VAR$(ind) obtained. Diverting...")
    return true
end
