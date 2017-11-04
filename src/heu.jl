"""
    A simple heuristic that round realxation solutions to seek for feasible solutions
"""
function heu_relaxation_rounding(m::PODNonlinearModel)

    !(:Bin in m.var_type_orig || :Int in m.var_type_orig) && return # No need to solve the continuous problem
    m.nlp_local_solver == UnsetSolver() && return # Always require a local nonlinear solver

    local_solve_nlp_model = MathProgBase.NonlinearModel(m.nlp_local_solver)
    l_var, u_var = fix_domains(m)
    MathProgBase.loadproblem!(local_solve_nlp_model,
                              m.num_var_orig,
                              m.num_constr_orig,
                              l_var,
                              u_var,
                              m.l_constr_orig,
                              m.u_constr_orig,
                              m.sense_orig,
                              m.d_orig)

    (!m.d_orig.want_hess) && MathProgBase.initialize(m.d_orig, [:Grad,:Jac,:Hess,:HessVec,:ExprGraph]) # Safety scheme for sub-solvers re-initializing the NLPEvaluator
    MathProgBase.setwarmstart!(local_solve_nlp_model, m.best_sol[1:m.num_var_orig])

    start_heuristic_solve = time()
    MathProgBase.optimize!(local_solve_nlp_model)
    m.logs[:total_time] += cputime_local_solve
    m.logs[:time_left] = max(0.0, m.timeout - m.logs[:total_time])
    status_pass = [:Optimal, :Suboptimal, :UserLimit, :LocalOptimal]
    status_reroute = [:Infeasible]

    local_solve_nlp_status = MathProgBase.status(local_solve_nlp_model)
    if local_solve_nlp_status in status_pass
        relaxed_sol = MathProgBase.getsolution(local_solve_nlp_model)
        rounded_sol = [m.var_type_orig[i] in [:Bin, :Int] ? round(relaxed_sol[i]) : round(relaxed_sol[i], 6) for i in 1:num_var_orig]
        eval_rhs = zeros(m.num_var_orig)
        MathProgBase.eval_g(m.d_orig, eval_rhs, rounded_sol)
        feasible = true
        for i in 1:m.num_constr_orig
            if m.constr_type_orig[i] == :(==)
                !isapprox(eval_rhs[i], m.constr_lb_orig[i];atol=m.tol_fea) && (feasible = false) && break
            elseif m.constr_type_orig[i] == :(>=)
                !(eval_rhs[i] >= m.constr_lb_orig[i] - m.tol_fea) && (feasible = false) && break
            elseif m.constr_type_orig[i] == :(<=)
                !(eval_rhs[i] <= m.constr_ub_orig[i] + m.tol_fea) && (feasible = false) && break
            end
        end
        if feasible
            candidate_obj = MathProgBase.eval_f(m.d_orig, rounded_sol) # Re-obtain the objective value
            push!(m.logs[:obj], candidate_obj)
            if eval(convertor[m.sense_orig])(candidate_obj, m.best_obj + 1e-5)
                m.best_obj = candidate_obj
                m.best_sol = MathProgBase.getsolution(local_solve_nlp_model)
                m.best_sol = round.(MathProgBase.getsolution(local_solve_nlp_model), 5)
                m.status[:feasible_solution] = :Detected
            end
            m.status[:local_solve] = local_solve_nlp_status
            return :UserLimit
        end
        return
    elseif local_solve_nlp_status in status_reroute
        return :Infeasible
    elseif local_solve_nlp_status == :Unbounded
        warn("[ROUNDING HEURISTIC] NLP local solve is unbounded.")
        return :Unbounded
    else
        warn("[ROUNDING HEURISTIC] NLP local solve failure.")
        return :Error
    end

    return :Infeasible
end
