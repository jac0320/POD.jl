function measure_relaxed_deviation(m::PODNonlinearModel;sol=nothing)

    sol == nothing ? sol = m.best_bound_sol : sol = sol

    isempty(sol) && return

    dev = []
    for k in keys(m.nonconvex_terms)
        y_idx = m.nonconvex_terms[k][:y_idx]
        y_hat = sol[y_idx]
        y_val = m.nonconvex_terms[k][:evaluator](m.nonconvex_terms[k], sol)
        y_ext = measure_extreme_points(m, k, m.discretization, sol)
        push!(dev, (y_idx, abs(y_hat-y_val), y_hat, y_val, m.nonconvex_terms[k][:var_idxs], y_ext, [sol[j] for j in m.nonconvex_terms[k][:var_idxs]]))
    end

    sort!(dev, by=x->x[1])
    println("====================================================================")
    for i in 1:m.num_var_orig
        local_interval = find_local_partition(m.discretization[i], sol[i])
        local_binding = false
        for j in local_interval
            if isapprox(sol[i], j;atol=1e-6)
                local_binding = true
            end
        end
        m.loglevel > 99 && println("X-VAR$(i): X-val=$(sol[i]) || $(local_interval) || BIND=$(local_binding)")
    end
    for i in dev
        m.loglevel > 99 && println("Y-VAR$(i[1]): DIFF=$(i[2]) || Y-hat=$(i[3]), Y-val=$(i[4]) || EXT=$(i[6]) || COMP $(i[5]) || X $(i[7])")
    end
    println("====================================================================")

    return
end

function measure_extreme_points(m::PODNonlinearModel, nlk::Any, d::Dict, sol::Vector)

    lifted_idx = m.nonconvex_terms[nlk][:lifted_var_ref].args[2]

    cnt = 0
    bound = []
    for var in nlk
        cnt += 1
        var_idx = var.args[2]
        var_bounds = find_local_partition(d[var_idx], sol[var_idx])
        if cnt == 1
            bound = copy(var_bounds)
        elseif cnt == 2
            bound = bound * var_bounds'
        else
            bound = diag(bound) * var_bounds'
        end
    end

    return bound
end

function measure_binding_constraints(m::PODNonlinearModel;sol=nothing)

    sol == nothing ? sol = m.best_bound_sol : sol = sol

    isempty(sol) && return

    bound_dev = []
    for i in 1:m.num_constr_orig
        constr = m.bounding_constr_mip[i]
        lhs = 0.0
        for j in 1:constr[:cnt]
            lhs += constr[:coefs][j] * sol[constr[:vars][j].args[2]]
        end
        rhs = constr[:rhs]
        push!(bound_dev, (i, abs(lhs-rhs), isapprox(lhs, rhs; atol=1e-8)))
    end

    orig_dev = []
    evaluated_lhs = eval_lhs_orig(m, sol)
    for i in 1:m.num_constr_orig
        evaluated = evaluated_lhs[i]
        original = nothing
        violation = nothing
        if m.constr_type_orig[i] == :(==)
            original = m.l_constr_orig[i]
        elseif m.constr_type_orig[i] == :(>=)
            original = m.l_constr_orig[i]
        elseif m.constr_type_orig[i] == :(<=)
            original = m.u_constr_orig[i]
        end
        push!(orig_dev, (i, abs(evaluated-original), isapprox(evaluated, original; atol=1e-4)))
    end

    nlvar_in_constraint = find_nlvar_in_constr(m)

    println("====================================================================")
    for i in 1:m.num_constr_orig
        b = bound_dev[i]
        o = orig_dev[i]
        m.loglevel > 99 && println("CONS-$(i) | CONVEX $(m.constr_structure[i]==:convex) | NLP DIFF=$(round(o[2],9)) BINDING=$(o[3]) || MILP DIFF=$(round(b[2],7)) BINDING=$(b[3]) || $(nlvar_in_constraint[i])")
    end
    println("====================================================================")

    return
end

function eval_lhs_orig(m::PODNonlinearModel, sol::Vector)

    eval_lhs = zeros(m.num_constr_orig)
    interface_eval_g(m.d_orig, eval_lhs, sol[1:m.num_var_orig])

    return eval_lhs
end


"""
    Evaluate a solution feasibility: Solution bust be in the feasible category and evaluated rhs must be feasible
"""
function eval_feasibility(m::PODNonlinearModel, sol::Vector)

    length(sol) == m.num_var_orig || error("Candidate solution length mismatch.")

    for i in 1:m.num_var_orig
        # Check solution category and bounds
        if m.var_type[i] == :Bin
            isapprox(sol[i], 1.0;atol=m.tol) || isapprox(sol[i], 0.0;atol=m.tol) || return false
        elseif m.var_type[i] == :Int
            isapprox(mod(sol[i], 1.0), 0.0;atol=m.tol) || return false
        end
        # Check solution bounds (with tight bounds)
        sol[i] <= m.l_var_tight[i] - m.tol || return false
        sol[i] >= m.u_var_tight[i] + m.tol || return false
    end

    # Check constraint violation
    eval_rhs = eval_rhs_orig(m, rounded_sol)
    feasible = true
    for i in 1:m.num_constr_orig
        if m.constr_type_orig[i] == :(==)
            if !isapprox(eval_rhs[i], m.l_constr_orig[i]; atol=m.tol)
                feasible = false
                m.loglevel >= 100 && println("[BETA] Violation on CONSTR $(i) :: EVAL $(eval_rhs[i]) != RHS $(m.l_constr_orig[i])")
                m.loglevel >= 100 && println("[BETA] CONSTR $(i) :: $(m.bounding_constr_expr_mip[i])")
                return false
            end
        elseif m.constr_type_orig[i] == :(>=)
            if !(eval_rhs[i] >= m.l_constr_orig[i] - m.tol)
                m.loglevel >= 100 && println("[BETA] Violation on CONSTR $(i) :: EVAL $(eval_rhs[i]) !>= RHS $(m.l_constr_orig[i])")
                m.loglevel >= 100 && println("[BETA] CONSTR $(i) :: $(m.bounding_constr_expr_mip[i])")
                return false
            end
        elseif m.constr_type_orig[i] == :(<=)
            if !(eval_rhs[i] <= m.u_constr_orig[i] + m.tol)
                m.loglevel >= 100 && println("[BETA] Violation on CONSTR $(i) :: EVAL $(eval_rhs[i]) !<= RHS $(m.u_constr_orig[i])")
                m.loglevel >= 100 && println("[BETA] CONSTR $(i) :: $(m.bounding_constr_expr_mip[i])")
                return false
            end
        end
    end

    return feasible
end
