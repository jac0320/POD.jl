"""
    A simple heuristic that round realxation solutions to seek for feasible solutions
"""
function heu_relaxation_rounding(m::PODNonlinearModel)

    println("[BETA] Running rounding heuristic to retrive a feasible solution...")
    !(:Bin in m.var_type_orig || :Int in m.var_type_orig) && return # No need to solve the continuous problem
    m.nlp_local_solver == UnsetSolver() && return # Always require a local nonlinear solver

    start_heuristic_solve = time()
    local_solve_nlp_model = MathProgBase.NonlinearModel(m.nlp_local_solver)
    l_var, u_var = fix_domains(m, relaxed=true)
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
    if isempty(m.sol_incumb_lb)
        MathProgBase.setwarmstart!(local_solve_nlp_model, m.best_sol[1:m.num_var_orig])
    else
        MathProgBase.setwarmstart!(local_solve_nlp_model, m.sol_incumb_lb[1:m.num_var_orig])
    end
    println("[BETA] Resolving relaxed local NLP ...")
    MathProgBase.optimize!(local_solve_nlp_model)
    heuristic_solve_time = time() - start_heuristic_solve
    m.logs[:total_time] += heuristic_solve_time
    m.logs[:time_left] = max(0.0, m.timeout - m.logs[:total_time])
    status_pass = [:Optimal, :Suboptimal, :UserLimit, :LocalOptimal]
    status_reroute = [:Infeasible]

    local_solve_nlp_status = MathProgBase.status(local_solve_nlp_model)
    println("[BETA] Resolve status = $(local_solve_nlp_status)")
    if local_solve_nlp_status in status_pass
        relaxed_sol = MathProgBase.getsolution(local_solve_nlp_model)
        rounded_sol = [m.var_type_orig[i] in [:Bin, :Int] ? round(relaxed_sol[i]) : relaxed_sol[i] for i in 1:m.num_var_orig]
        eval_rhs = zeros(m.num_constr_orig)
        MathProgBase.eval_g(m.d_orig, eval_rhs, rounded_sol)
        feasible = true
        for i in 1:m.num_constr_orig
            if m.constr_type_orig[i] == :(==)
                if !isapprox(eval_rhs[i], m.l_constr_orig[i]; atol=m.tol_fea)
                    feasible = false
                    println("[BETA] Violation on CONSTR $(i) :: EVAL $(eval_rhs[i]) != RHS $(m.l_constr_orig[i])")
                    println("[BETA] CONSTR $(i) :: $(m.bounding_constr_expr_mip[i])")
                    break
                end
            elseif m.constr_type_orig[i] == :(>=)
                if !(eval_rhs[i] >= m.l_constr_orig[i] - m.tol_fea)
                    feasible = false
                    println("[BETA] Violation on CONSTR $(i) :: EVAL $(eval_rhs[i]) !>= RHS $(m.l_constr_orig[i])")
                    println("[BETA] CONSTR $(i) :: $(m.bounding_constr_expr_mip[i])")
                    break
                end
            elseif m.constr_type_orig[i] == :(<=)
                if !(eval_rhs[i] <= m.u_constr_orig[i] + m.tol_fea)
                    feasible = false
                    println("[BETA] Violation on CONSTR $(i) :: EVAL $(eval_rhs[i]) !<= RHS $(m.u_constr_orig[i])")
                    println("[BETA] CONSTR $(i) :: $(m.bounding_constr_expr_mip[i])")
                    break
                end
            end
        end

        if feasible
            println("[BETA] Rounded solution FEASIBLE. Recording...")
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
        println("[BETA] Rounded solution INfeasible. Recording...")
        push!(m.logs[:obj], "Infea")
        return :Infeasible
    elseif local_solve_nlp_status in status_reroute
        push!(m.logs[:obj], "Infea")
        return :Infeasible
    elseif local_solve_nlp_status == :Unbounded
        warn("[ROUNDING HEURISTIC] NLP local solve is unbounded.")
        push!(m.logs[:obj], "Ubd")
        return :Unbounded
    else
        warn("[ROUNDING HEURISTIC] NLP local solve failure.")
        push!(m.logs[:obj], "Err")
        return :Error
    end

end
