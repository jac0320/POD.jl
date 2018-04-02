"""
    update_rel_gap(m::PODNonlinearModel)

Update POD model relative & absolute optimality gap.

The relative gap calculation is

```math
    \\textbf{Gap} = \\frac{|UB-LB|}{ϵ+|UB|}
```

The absolute gap calculation is
```
    |UB-LB|
```
"""
function update_opt_gap(m::PODNonlinearModel)

    if m.best_obj in [Inf, -Inf]
        m.best_rel_gap = Inf
        return
    else
        p = convert(Int, round(abs(log(10,m.relgap))))
        n = round(abs(m.best_obj-m.best_bound), p)
        dn = round(abs(1e-12+abs(m.best_obj)), p)
        if isapprox(n, 0.0;atol=m.tol) && isapprox(m.best_obj,0.0;atol=m.tol)
            m.best_rel_gap = 0.0
            return
        end
        if m.gapref == "ub"
            m.best_rel_gap = abs(m.best_obj - m.best_bound)/(m.tol+abs(m.best_obj))
        else
            m.best_rel_gap = abs(m.best_obj - m.best_bound)/(m.tol+abs(m.best_bound))
        end
    end

    m.best_abs_gap = abs(m.best_obj - m.best_bound)
    return
end

"""
    discretization_to_bounds(d::Dict, l::Int)

    Same as [`update_var_bounds`](@ref)
"""
discretization_to_bounds(d::Dict, l::Int) = update_var_bounds(d, len=l)

"""
    Update the data structure with feasible solution and its associated objective (if better)
"""
function update_incumb_objective(m::PODNonlinearModel, objval::Float64, sol::Vector)

    convertor = Dict(:Max=>:>, :Min=>:<)
    push!(m.logs[:obj], objval)
    if eval(convertor[m.sense_orig])(objval, m.best_obj) #&& !eval(convertor[m.sense_orig])(objval, m.best_bound)
        m.best_obj = objval
        m.best_sol = sol
        m.status[:feasible_solution] = :Detected
    end

    return
end

"""
    Utility function for debugging.
"""
function show_solution(m::JuMP.Model)
    for i in 1:length(m.colNames)
        println("$(m.colNames[i])=$(m.colVal[i])")
    end
    return
end

function find_local_partition(part::Vector, sol::Float64)

    P = length(part) - 1
    for i in 1:P
        if sol >= part[i] - 1e-6 && sol <= part[i+1] + 1e-6
            return [part[i], part[i+1];]
        end
    end

    error("$(sol) $(part) Solution didn't showed in any partition")

    return
end


"""
    @docstring
"""
function insert_timeleft_symbol(options, val::Any, keywords::Symbol, timeout; options_string_type=1)

    for i in 1:length(options)
        if options_string_type == 1
            if keywords in collect(options[i])
                deleteat!(options, i)
                break
            end
        elseif options_string_type == 2
            if keywords == split(options[i],"=")[1]
                deleteat!(options, i)
                break
            end
        end
    end

    if options_string_type == 1
        (val != Inf) && push!(options, (keywords, val))
    elseif options_string_type == 2
        (val != Inf) && push!(options, "$(keywords)=$(val)")
    end

    return
end

"""
    fetch_boundstop_symbol(m::PODNonlinearModel)

An utility function used to recongize different sub-solvers and return the bound stop option key words
"""
function update_boundstop_options(m::PODNonlinearModel)

    if m.mip_solver_id == "Gurobi"
        # Calculation of the bound
        if m.sense_orig == :Min
            m.gapref == "ub" ? stopbound=(1-m.relgap+m.tol)*abs(m.best_obj) : stopbound=(1-m.relgap+m.tol)*abs(m.best_bound)
        elseif m.sense_orig == :Max
            m.gapref == "ub" ? stopbound=(1+m.relgap-m.tol)*abs(m.best_obj) : stopbound=(1+m.relgap-m.tol)*abs(m.best_bound)
        end

        for i in 1:length(m.mip_solver.options)
            if m.mip_solver.options[i][1] == :BestBdStop
                deleteat!(m.mip_solver.options, i)
                if m.mip_solver_id == "Gurobi"
                    push!(m.mip_solver.options, (:BestBdStop, stopbound))
                else
                    return
                end
            end
        end
    end

    return
end


"""
    is_fully_convexified(m::PODNonlinearModel)
"""
function is_fully_convexified(m::PODNonlinearModel)

    # Other more advanced convexification check goes here
    for term in keys(m.nonconvex_terms)
        if !m.nonconvex_terms[term][:convexified]
            warn("Detected terms that is not convexified $(term[:lifted_constr_ref]), bounding model solver may report a error due to this")
            return
        else
            m.nonconvex_terms[term][:convexified] = false    # Reset status for next iteration
        end
    end

    return
end



"""
    Need to be careful with the diverted point.
"""
function track_new_partition_idx(d::Dict, idx::Int, val::Float64, acp::Int=-1)

    acp > 0 ? acp = acp : acp = get_active_partition_idx(d,val,idx)

    pcnt = length(d[idx]) - 1
    newpidx = []                # Tracks the newly constructed partition idxes
    pcnt == 1 && return [1;]    # Still keep non-discretizing variables
    if acp == 1
        return newpidx = [1,2;]
    elseif acp == pcnt
        return newpidx = [pcnt-1,pcnt;]
    else
        tlb = d[idx][acp-1]
        tub = d[idx][acp+1]
        if abs(val-tlb) == abs(val-tub)
            return [acp-1, acp, acp+1;]
        elseif abs(val-tlb) > abs(val-tub)
            return [acp-1, acp;]
        elseif abs(val-tlb) < abs(val-tub)
            return [acp, acp+1;]
        end
    end

    return
end

"""
    Collect active partition idx
    Need to be careful with the diverted point
"""
function get_active_partition_idx(discretization::Dict, val::Float64, idx::Int; tol=1e-6)

    for j in 1:length(discretization[idx])-1
        if val > discretization[idx][j] - tol && val < discretization[idx][j+1] + tol
            return j
        end
    end

    warn("Activate parition not found [VAR$(idx)]. Returning default partition 1.")
    return 1
end


function eval_objective(m::PODNonlinearModel; svec::Vector=[])

    isempty(svec) ? svec = m.best_bound_sol : svec = svec
    m.sense_orig == :Min ? obj = Inf : obj=-Inf

    if m.obj_structure == :affine
        obj = m.bounding_obj_mip[:rhs]
        for i in 1:m.bounding_obj_mip[:cnt]
            obj += m.bounding_obj_mip[:coefs][i]*svec[m.bounding_obj_mip[:vars][i].args[2]]
        end
    elseif m.structural_obj == :convex
        error("need implementation for local objective function evaluation for convex form")
    else
        error("Unknown structural obj type $(m.structural_obj)")
    end

    return obj
end

function adjust_branch_priority(m::PODNonlinearModel, priority=nothing)

    priority == nothing ? priority = m.convhull_branch_priority : priority = priority

    isempty(priority) && return
    m.mip_solver_id != "Gurobi" && return
    !m.model_mip.internalModelLoaded && JuMP.build(m.model_mip)

    len = length(m.model_mip.colVal)
    Gurobi.set_intattrarray!(m.model_mip.internalModel.inner, "BranchPriority", 1, len, [priority[i] for i in 1:len])

    return
end

"""
    Special funtion for debugging bounding models
"""
function print_iis_gurobi(m::JuMP.Model)

    grb = interface_get_rawsolver(m)
    Gurobi.computeIIS(grb)
    numconstr = Gurobi.num_constrs(grb)
    numvar = Gurobi.num_vars(grb)

    iisconstr = Gurobi.get_intattrarray(grb, "IISConstr", 1, numconstr)
    iislb = Gurobi.get_intattrarray(grb, "IISLB", 1, numvar)
    iisub = Gurobi.get_intattrarray(grb, "IISUB", 1, numvar)

    info("Irreducible Inconsistent Subsystem (IIS)")
    info("Variable bounds:")
    for i in 1:numvar
        v = Variable(m, i)
        if iislb[i] != 0 && iisub[i] != 0
            println(getlowerbound(v), " <= ", getname(v), " <= ", getupperbound(v))
        elseif iislb[i] != 0
            println(getname(v), " >= ", getlowerbound(v))
        elseif iisub[i] != 0
            println(getname(v), " <= ", getupperbound(v))
        end
    end

    info("Constraints:")
    for i in 1:numconstr
        if iisconstr[i] != 0
            println(m.linconstr[i])
        end
    end

    return
end

function round_sol(m::PODNonlinearModel;nlp_model=nothing, nlp_sol=[])

    if nlp_model != nothing
        relaxed_sol = interface_get_solution(nlp_model)
    end

    if !isempty(nlp_sol)
        relaxed_sol = nlp_sol
    end

    if nlp_model != nothing && !isempty(nlp_sol)
        error("In function collision. Special usage")
    end

    rounded_sol = copy(relaxed_sol)
    for i in 1:m.num_var_orig
        if m.var_type_orig[i] == :Bin
            relaxed_sol[i] >= 0.5 ? rounded_sol[i] = 1 : rounded_sol[i] = 0
        elseif m.var_type_orig[i] == :Int
            rounded_sol[i] = round(relaxed_sol[i])
        else
            rounded_sol[i] = relaxed_sol[i]
        end
    end

    return rounded_sol
end

function fetch_mip_solver_identifier(m::PODNonlinearModel;override="")

    isempty(override) ? solverstring = string(m.mip_solver) : solverstring=override

    # Higher-level solvers: that can use sub-solvers
    if contains(solverstring,"Pajarito")
        m.mip_solver_id = "Pajarito"
        return
    end

    # Lower level solvers
    if contains(solverstring,"Gurobi")
        m.mip_solver_id = "Gurobi"
    elseif contains(solverstring,"CPLEX")
        m.mip_solver_id = "CPLEX"
    elseif contains(solverstring,"Cbc")
        m.mip_solver_id = "Cbc"
    elseif contains(solverstring,"GLPK")
        m.mip_solver_id = "GLPK"
    else
        error("Unsupported mip solver name. Using blank")
    end

    return
end

function fetch_nlp_solver_identifier(m::PODNonlinearModel;override="")

    isempty(override) ? solverstring = string(m.nlp_solver) : solverstring=override

    # Higher-level solver
    if contains(solverstring, "Pajarito")
        m.nlp_solver_id = "Pajarito"
        return
    end

    # Lower-level solver
    if contains(solverstring, "Ipopt")
        m.nlp_solver_id = "Ipopt"
    elseif contains(solverstring, "AmplNL") && contains(solverstring, "bonmin")
        m.nlp_solver_id = "Bonmin"
    elseif contains(solverstring, "KNITRO")
        m.nlp_solver_id = "Knitro"
    elseif contains(solverstring, "NLopt")
        m.nlp_solver_id = "NLopt"
    else
        error("Unsupported nlp solver name. Using blank")
    end

    return
end

function fetch_minlp_solver_identifier(m::PODNonlinearModel;override="")

    (m.minlp_solver == UnsetSolver()) && return

    isempty(override) ? solverstring = string(m.minlp_solver) : solverstring=override

    # Higher-level solver
    if contains(solverstring, "Pajarito")
        m.minlp_solver_id = "Pajarito"
        return
    end

    # Lower-level Solver
    if contains(solverstring, "AmplNL") && contains(solverstring, "bonmin")
        m.minlp_solver_id = "Bonmin"
    elseif contains(solverstring, "KNITRO")
        m.minlp_solver_id = "Knitro"
    elseif contains(solverstring, "NLopt")
        m.minlp_solver_id = "NLopt"
    elseif contains(solverstring, "CoinOptServices.OsilSolver(\"bonmin\"")
        m.minlp_solver_id = "Bonmin"
    else
        error("Unsupported nlp solver name. Using blank")
    end

    return
end

"""

    update_mip_time_limit(m::PODNonlinearModel)

An utility function used to dynamically regulate MILP solver time limits to fit POD solver time limits.
"""
function update_mip_time_limit(m::PODNonlinearModel; kwargs...)

    options = Dict(kwargs)
    haskey(options, :timelimit) ? timelimit = options[:timelimit] : timelimit = max(0.0, m.timeout-m.logs[:total_time])

    if m.mip_solver_id == "CPLEX"
        insert_timeleft_symbol(m.mip_solver.options,timelimit,:CPX_PARAM_TILIM,m.timeout)
    elseif m.mip_solver_id == "Gurobi"
        insert_timeleft_symbol(m.mip_solver.options,timelimit,:TimeLimit,m.timeout)
    elseif m.mip_solver_id == "Cbc"
        insert_timeleft_symbol(m.mip_solver.options,timelimit,:seconds,m.timeout)
    elseif m.mip_solver_id == "GLPK"
        insert_timeleft_symbol(m.mip_solver.opts, timelimit,:tm_lim,m.timeout)
    elseif m.mip_solver_id == "Pajarito"
        (timelimit < Inf) && (m.mip_solver.timeout = timelimit)
    else
        error("Needs support for this MIP solver")
    end

    return
end

function mip_solver_verbosity(m::PODNonlinearModel, verbosity::Int)

    if m.mip_solver_id == "CPLEX"
        insert_timeleft_symbol(m.mip_solver.options,verbosity,:CPX_PARAM_SCRIND,m.timeout)
    elseif m.mip_solver_id == "Gurobi"
        insert_timeleft_symbol(m.mip_solver.options,verbosity,:OutputFlag,m.timeout)
    else
        error("Needs support for this MIP solver")
    end

    return
end

"""

    update_mip_time_limit(m::PODNonlinearModel)

An utility function used to dynamically regulate MILP solver time limits to fit POD solver time limits.
"""
function update_nlp_time_limit(m::PODNonlinearModel; kwargs...)

    options = Dict(kwargs)
    haskey(options, :timelimit) ? timelimit = options[:timelimit] : timelimit = max(0.0, m.timeout-m.logs[:total_time])

    if m.nlp_solver_id == "Ipopt"
        insert_timeleft_symbol(m.nlp_solver.options,timelimit,:CPX_PARAM_TILIM,m.timeout)
    elseif m.nlp_solver_id == "Pajarito"
        (timelimit < Inf) && (m.nlp_solver.timeout = timelimit)
    elseif m.nlp_solver_id == "AmplNL"
        insert_timeleft_symbol(m.nlp_solver.options,timelimit,:seconds,m.timeout, options_string_type=2)
    elseif m.nlp_solver_id == "Knitro"
        error("You never tell me anything about knitro. Probably because they have a very short trail length.")
    elseif m.nlp_solver_id == "NLopt"
        m.nlp_solver.maxtime = timelimit
    else
        error("Needs support for this MIP solver")
    end

    return
end

"""

    update_mip_time_limit(m::PODNonlinearModel)

    An utility function used to dynamically regulate MILP solver time limits to fit POD solver time limits.
"""
function update_minlp_time_limit(m::PODNonlinearModel; kwargs...)

    options = Dict(kwargs)
    haskey(options, :timelimit) ? timelimit = options[:timelimit] : timelimit = max(0.0, m.timeout-m.logs[:total_time])

    if m.minlp_solver_id == "Pajarito"
        (timelimit < Inf) && (m.minlp_solver.timeout = timelimit)
    elseif m.minlp_solver_id == "AmplNL"
        insert_timeleft_symbol(m.minlp_solver.options,timelimit,:seconds,m.timeout,options_string_type=2)
    elseif m.minlp_solver_id == "Knitro"
        error("You never tell me anything about knitro. Probably because they charge everything they own.")
    elseif m.minlp_solver_id == "NLopt"
        m.minlp_solver.maxtime = timelimit
    else
        error("Needs support for this MIP solver")
    end

    return
end
