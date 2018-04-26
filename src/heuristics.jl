"""
    One-time rounding heuristic to obtain a feasible solution
    For integer solutions
"""
function heu_basic_rounding(m::PODNonlinearModel, local_model)

    println("Basic Rounding Heuristic Avtivated...")

    convertor = Dict(:Max=>:>, :Min=>:<)

    rounded_sol = round_sol(m, nlp_model=local_model)
    l_var, u_var = fix_domains(m, discrete_sol = rounded_sol)

    heuristic_model = interface_init_nonlinear_model(m.nlp_solver)
    interface_load_nonlinear_model(m, heuristic_model, l_var, u_var)
    interface_optimize(heuristic_model)
    heuristic_model_status = interface_get_status(heuristic_model)

    if heuristic_model_status in [:Infeasible, :Error]
        m.loglevel > 0 && println("Rounding obtained an Infeasible point.")
        push!(m.logs[:obj], "INF")
        return :Infeasibles
    elseif heuristic_model_status in [:Optimal, :Suboptimal, :LocalOptimal]
        candidate_obj = interface_get_objval(heuristic_model)
        candidate_sol = round.(interface_get_solution(heuristic_model), 5)
        update_incumb_objective(m, candidate_obj, candidate_sol)
        m.loglevel > 0 && println("Rounding obtained a feasible solution OBJ = $(m.best_obj)")
        return :LocalOptimal
    else
        error("[EXCEPTION] Unknown NLP solver status.")
    end

    return
end

"""
    Use all lower bound solution pools as starting points
"""
function heu_pool_multistart(m::PODNonlinearModel)

    convertor = Dict(:Max=>:>, :Min=>:<)
    m.sense_orig == :Min ? incumb_obj = Inf : incumb_obj = -Inf
    incumb_sol = []
    found_feasible = false

    for i in 1:m.bound_sol_pool[:cnt]
        if !m.bound_sol_pool[:ubstart][i]
            rounded_sol = round_sol(m, nlp_sol=m.bound_sol_pool[:sol][i])
            l_var, u_var = fix_domains(m, discrete_sol=rounded_sol, use_orig=true)
            heuristic_model = interface_init_nonlinear_model(m.nlp_solver)
            interface_load_nonlinear_model(m, heuristic_model, l_var, u_var)
            interface_optimize(heuristic_model)
            heuristic_model_status = interface_get_status(heuristic_model)
            if heuristic_model_status in [:Optimal, :Suboptimal, :LocalOptimal]
                candidate_obj = interface_get_objval(heuristic_model)
                if eval(convertor[m.sense_orig])(candidate_obj, incumb_obj)
                    incumb_obj = candidate_obj
                    incumb_sol = round.(interface_get_solution(heuristic_model), 5)
                    m.loglevel > 0 && println("Feasible solution obtained using lower bound solution pool [SOL:$(i)] [OBJ=$(incumb_obj)]")
                end
                found_feasible = true
            else
                m.loglevel > 199 && println("Multi-start heuristic returns $(heuristic_model_status) [SOL:$(i)]")
            end
            m.bound_sol_pool[:ubstart][i] = true
        end
    end

    if found_feasible
        update_incumb_objective(m, incumb_obj, incumb_sol)
        return :LocalOptimal
    end

    return :Infeasible
end
