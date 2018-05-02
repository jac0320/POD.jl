function create_bounding_mip(m::PODNonlinearModel; use_disc=nothing, warmstart=true, cuts=true)

    use_disc == nothing ? discretization = m.discretization : discretization = use_disc

    m.model_mip = Model(solver=m.mip_solver) # Construct JuMP Model
    start_build = time()
    # ------- Model Construction ------ #
    amp_post_vars(m)                                      # Post original and lifted variables
    amp_post_lifted_constraints(m)                        # Post lifted constraints
    amp_post_lifted_objective(m)                          # Post objective
    amp_post_convexification(m, use_disc=discretization, warmstart=warmstart, cuts=cuts)
    # --------------------------------- #
    cputime_build = time() - start_build
    m.logs[:total_time] += cputime_build
    m.logs[:time_left] = max(0.0, m.timeout - m.logs[:total_time])

    return
end

function create_bounding_slackness_mip(m::PODNonlinearModel; use_disc=nothing, warmstart=true, cuts=true)

    use_disc == nothing ? discretization = m.discretization : discretization = use_disc

    m.model_mip = Model(solver=m.mip_solver) # Construct JuMP Model
    start_build = time()
    # ------- Model Construction ------ #
    amp_post_vars(m, enable_slack=true)                            # Post original and lifted variables
    amp_post_lifted_constraints(m, enable_slack=true)              # Post lifted constraints
    amp_post_convexification(m, use_disc=discretization, warmstart=warmstart, cuts=cuts)
                                                                   # Convexify problem
    amp_post_slackness_objective(m)                                # Post objective
    # --------------------------------- #
    cputime_build = time() - start_build
    m.logs[:total_time] += cputime_build
    m.logs[:time_left] = max(0.0, m.timeout - m.logs[:total_time])

    return
end

function create_bounding_feasibility_mip(m::PODNonlinearModel; use_disc=nothing, warmstart=true, cuts=true)

    use_disc == nothing ? discretization = m.discretization : discretization = use_disc

    m.model_mip = Model(solver=m.mip_solver) # Construct JuMP Model
    start_build = time()
    # ------- Model Construction ------ #
    amp_post_vars(m, enable_slack=true)                            # Post original and lifted variables
    amp_post_lifted_constraints(m, enable_slack=true)              # Post lifted constraints
    amp_post_convexification(m, use_disc=discretization, warmstart=warmstart, cuts=cuts)           # Convexify problem
    amp_post_zero_objective(m)                                     # Post objective
    # --------------------------------- #
    cputime_build = time() - start_build
    m.logs[:total_time] += cputime_build
    m.logs[:time_left] = max(0.0, m.timeout - m.logs[:total_time])

end

"""
    amp_post_convexification(m::PODNonlinearModel; kwargs...)

warpper function to convexify the problem for a bounding model. This function talks to nonconvex_terms and convexification methods
to finish the last step required during the construction of bounding model.
"""
function amp_post_convexification(m::PODNonlinearModel; use_disc=nothing, warmstart=true, cuts=true)

    use_disc == nothing ? discretization = m.discretization : discretization = use_disc

    for i in 1:length(m.method_convexification)             # Additional user-defined convexification method
        eval(m.method_convexification[i])(m)
    end

    amp_post_mccormick(m, use_disc=discretization)              # handles all bi-linear and monomial convexificaitons
    amp_post_convhull(m, use_disc=discretization, warmstart=warmstart, cuts=cuts)    # convex hull representation

    is_fully_convexified(m) # Exam to see if all non-linear terms have been convexificed

    return
end

function amp_post_vars(m::PODNonlinearModel;enable_slack=false, kwargs...)

    options = Dict(kwargs)

    if haskey(options, :use_disc)
        l_var = [options[:use_disc][i][1]   for i in 1:(m.num_var_orig+m.num_var_linear_mip+m.num_var_nonlinear_mip)]
        u_var = [options[:use_disc][i][end] for i in 1:(m.num_var_orig+m.num_var_linear_mip+m.num_var_nonlinear_mip)]
    else
        l_var = m.l_var_tight
        u_var = m.u_var_tight
    end

    @variable(m.model_mip, x[i=1:(m.num_var_orig+m.num_var_linear_mip+m.num_var_nonlinear_mip)])

    enable_slack && @variable(m.model_mip, s[i=1:m.num_slack_vars]>=0)

    for i in 1:(m.num_var_orig+m.num_var_linear_mip+m.num_var_nonlinear_mip)

        (i <= m.num_var_orig) && setcategory(x[i], m.var_type_orig[i])  # This is a tricky step, not enforcing category of lifted variables is able to improve performance
        l_var[i] > -Inf && setlowerbound(x[i], l_var[i])    # Changed to tight bound, if no bound tightening is performed, will be just .l_var_orig
        u_var[i] < Inf && setupperbound(x[i], u_var[i])    # Changed to tight bound, if no bound tightening is performed, will be just .u_var_orig

        #TODO experimental code for integer problems, make sure it is properly re-structured
        m.int_enable && m.var_type[i] == :Int && setcategory(x[i], :Cont)
        if m.int_enable && m.var_type[i] == :Int && i in m.disc_vars
            setlowerbound(x[i], floor(m.l_var_tight[i]) - 0.5)
        end
        if m.int_enable && m.var_type[i] == :Int && i in m.disc_vars
            setupperbound(x[i], ceil(m.u_var_tight[i]) + 0.5)
        end
    end

    return
end

function amp_post_lifted_constraints(m::PODNonlinearModel;enable_slack=false)

    for i in 1:m.num_constr_orig
        if enable_slack
            if m.constr_structure[i] == :affine
                amp_post_affine_slack_constraint(m, i)
            elseif m.constr_structure[i] == :convex
                warn("Support for convex constrint in feasibility mode is highly experimental.", once=true)
                amp_post_convex_constraint(m.model_mip, m.bounding_constr_mip[i])
            else
                error("Unknown constr_structure type $(m.constr_structure[i]). No support for convex constraints under feasibility_mode")
            end
        else
            if m.constr_structure[i] == :affine
                amp_post_affine_constraint(m.model_mip, m.bounding_constr_mip[i])
            elseif m.constr_structure[i] == :convex
                amp_post_convex_constraint(m.model_mip, m.bounding_constr_mip[i])
            else
                error("Unknown constr_structure type $(m.constr_structure[i])")
            end
        end
    end

    for i in keys(m.linear_terms)
        amp_post_linear_lift_constraints(m.model_mip, m.linear_terms[i])
    end

    return
end

function amp_post_affine_constraint(model_mip::JuMP.Model, affine::Dict; constr_id::Int=0)

    if affine[:sense] == :(>=)
        @constraint(model_mip, sum(affine[:coefs][j]*Variable(model_mip, affine[:vars][j].args[2]) for j in 1:affine[:cnt]) >= affine[:rhs])
    elseif affine[:sense] == :(<=)
        @constraint(model_mip, sum(affine[:coefs][j]*Variable(model_mip, affine[:vars][j].args[2]) for j in 1:affine[:cnt]) <= affine[:rhs])
    elseif affine[:sense] == :(==)
        @constraint(model_mip, sum(affine[:coefs][j]*Variable(model_mip, affine[:vars][j].args[2]) for j in 1:affine[:cnt]) == affine[:rhs])
    else
        error("Unkown sense.")
    end

    return
end

function amp_post_convex_constraint(model_mip::JuMP.Model, convex::Dict)

    !prod([i == 2 for i in convex[:powers]]) && error("No relaxation implementation for convex constraints $(convex[:expr])")

    if convex[:sense] == :(<=)
        @constraint(model_mip,
            sum(convex[:coefs][j]*Variable(model_mip, convex[:vars][j].args[2])^2 for j in 1:convex[:cnt]) <= convex[:rhs])
    elseif convex[:sense] == :(>=)
        @constraint(model_mip,
            sum(convex[:coefs][j]*Variable(model_mip, convex[:vars][j].args[2])^2 for j in 1:convex[:cnt]) >= convex[:rhs])
    else
        error("No equality constraints should be recognized as supported convex constriants")
    end

    return
end

function amp_post_linear_lift_constraints(model_mip::JuMP.Model, l::Dict)

    @assert l[:ref][:sign] == :+
    @constraint(model_mip, Variable(model_mip, l[:y_idx]) == sum(i[1]*Variable(model_mip, i[2]) for i in l[:ref][:coef_var]) + l[:ref][:scalar])
    return
end

function amp_post_lifted_objective(m::PODNonlinearModel)

    if m.obj_structure == :affine
        @objective(m.model_mip, m.sense_orig, m.bounding_obj_mip[:rhs] + sum(m.bounding_obj_mip[:coefs][i]*Variable(m.model_mip, m.bounding_obj_mip[:vars][i].args[2]) for i in 1:m.bounding_obj_mip[:cnt]))
    elseif m.obj_structure == :convex
        @objective(m.model_mip, m.sense_orig, m.bounding_obj_mip[:rhs] + sum(m.bounding_obj_mip[:coefs][i]*Variable(m.model_mip, m.bounding_obj_mip[:vars][i].args[2])^2 for i in 1:m.bounding_obj_mip[:cnt]))
    else
        error("Unknown structural obj type $(m.obj_structure)")
    end

    return
end

function add_partition(m::PODNonlinearModel; use_disc=nothing, use_solution=nothing, use_ratio=nothing, method=nothing)

    use_disc == nothing ? discretization = m.discretization : discretization = use_disc
    use_solution == nothing ? point_vec = m.best_bound_sol : point_vec = use_solution
    method == nothing ? method = m.disc_add_partition_method : method = method
    use_ratio == nothing ? ratio = m.disc_ratio : ratio = use_ratio

    if isa(method, Function)
        m.discretization = eval(method)(m, use_disc=discretization, use_solution=point_vec)
    elseif method == "adaptive"
        m.discretization = add_adaptive_partition(m, use_disc=discretization, use_solution=point_vec, use_ratio=use_ratio)
    elseif method == "separative"
        m.discretization = add_separative_partition(m, use_disc=discretization, user_solution=point_vec)
    elseif method == "uniform"
        m.discretization = add_uniform_partition(m, use_disc=discretization)
    else
        error("Unknown input on how to add partitions.")
    end

    return
end

function add_separative_partition(m::PODNonlinearModel;kwargs...)

    options = Dict(kwargs)

    haskey(options, :use_disc) ? discretization = options[:use_disc] : discretization = m.discretization
    haskey(options, :use_solution) ? point_vec = copy(options[:use_solution]) : point_vec = copy(m.best_bound_sol)
    haskey(options, :use_ratio) ? ratio = options[:use_ratio] : ratio = m.disc_ratio
    haskey(options, :branching) ? branching = options[:branching] : branching = false

    if length(point_vec) < m.num_var_orig + m.num_var_linear_mip + m.num_var_nonlinear_mip
        point_vec = resolve_lifted_var_value(m, point_vec)  # Update the solution vector for lifted variable
    end

    if branching
        discretization = deepcopy(discretization)
    end

    processed = Set{Int}()

    # ? Perform discretization base on type of nonlinear terms ? #
    for i in m.disc_vars
        point = point_vec[i]                # Original Variable
        λCnt = length(discretization[i])
        # Built-in method based-on variable type
        if m.var_type[i] == :Cont
            i in processed && continue
            point = correct_point(m, discretization[i], point, i)
            for j in 1:λCnt
                partvec = discretization[i]
                if point >= partvec[j] && point <= partvec[j+1]  # Locating the right location
                    insert_separative_partition(m, i, j, point, partvec)
                    push!(processed, i)
                    break
                end
            end
        else
            error("Unexpected variable types during injecting partitions")
        end
    end

    return discretization
end

function insert_separative_partition(m::PODNonlinearModel, var::Int, partidx::Int, point::Number, partvec::Vector)

    abstol, reltol = m.disc_abs_width_tol, m.disc_rel_width_tol

    lb_local, ub_local = partvec[partidx], partvec[partidx+1]
    ub_touch, lb_touch = true, true

    if point < ub_local && !isapprox(point, ub_local; atol=abstol) && abs(point-ub_local)/(1e-8+abs(ub_local)) > reltol # Insert new UB-based partition
        ub_touch = false
    end

    if point > lb_local && !isapprox(point, lb_local; atol=abstol) && abs(point-lb_local)/(1e-8+abs(lb_local)) > reltol # Insert new LB-based partition
        lb_touch = false
    end

    if ub_touch && lb_touch
        distvec = [(j, partvec[j+1]-partvec[j]) for j in 1:length(partvec)-1]
        sort!(distvec, by=x->x[2])
        point_orig = point
        pos = distvec[end][1]
        lb_local = partvec[pos]
        ub_local = partvec[pos+1]
        isapprox(lb_local, ub_local;atol=m.tol) && return
        chunk = (ub_local - lb_local) / m.disc_divert_chunks
        point = lb_local + (ub_local - lb_local) / m.disc_divert_chunks
        for i in 2:m.disc_divert_chunks
            insert!(partvec, pos+1, lb_local + chunk * (m.disc_divert_chunks-(i-1)))
        end
        (m.loglevel > 99) && println("[DEBUG] !D! VAR$(var): SOL=$(round(point_orig,4))=>$(point) |$(round(lb_local,4)) | $(m.disc_divert_chunks) SEGMENTS | $(round(ub_local,4))|")
    elseif ub_touch || lb_touch
        lb_local = partvec[partidx]
        ub_local = partvec[partidx+1]
        chunk = (ub_local - lb_local) / m.disc_divert_chunks
        point = lb_local + (ub_local - lb_local) / m.disc_divert_chunks
        for i in 2:m.disc_divert_chunks
            insert!(partvec, partidx+1, lb_local + chunk * (m.disc_divert_chunks-(i-1)))
        end
        (m.loglevel > 99) && println("[DEBUG] VAR$(var): SOL=$(round(point,4)) PARTITIONS=$(length(partvec)-1) |$(round(lb_local,4)) <- | $(m.disc_divert_chunks) SEGMENTS | -> $(round(ub_local,4))|")
        return
    else
        (m.loglevel > 99) && println("[DEBUG] VAR$(var): SOL=$(round(point,4)), PARTITIONS=$(length(partvec)-1) |$(round(lb_local,4)) <- | $(round(point,4)) | -> $(round(ub_local,4))|")
        insert!(partvec, partidx+1, point)
        return
    end

    return
end


"""
    add_discretization(m::PODNonlinearModel; use_disc::Dict, use_solution::Vector)

Basic built-in method used to add a new partition on feasible domains of discretizing variables.
This method make modification in .discretization

Consider original partition [0, 3, 7, 9], where LB/any solution is 4.
Use ^ as the new partition, "|" as the original partition

A case when discretize ratio = 4
| -------- | - ^ -- * -- ^ ---- | -------- |
0          3  3.5   4   4.5     7          9

A special case when discretize ratio = 2
| -------- | ---- * ---- ^ ---- | -------- |
0          3      4      5      7          9

There are two options for this function,

    * `use_disc(default=m.discretization)`:: to regulate which is the base to add new partitions on
    * `use_solution(default=m.best_bound_sol)`:: to regulate which solution to use when adding new partitions on

TODO: also need to document the speical diverted cases when new partition touches both sides

This function belongs to the hackable group, which means it can be replaced by the user to change the behvaior of the solver.
"""
function add_adaptive_partition(m::PODNonlinearModel;use_ratio=nothing, branching=nothing, use_disc=nothing, use_solution=nothing)

    use_disc == nothing ? discretization = m.discretization : discretization = use_disc
    use_solution == nothing ? point_vec = copy(use_solution) : point_vec = copy(use_solution)
    use_ratio == nothing ? ratio = m.disc_ratio : ratio = use_ratio
    branching == nothing ? branching = false : branching = branching

    if length(point_vec) < m.num_var_orig + m.num_var_linear_mip + m.num_var_nonlinear_mip
        point_vec = resolve_lifted_var_value(m, point_vec)  # Update the solution vector for lifted variable
    end

    if branching
        discretization = deepcopy(discretization)
    end

    processed = Set{Int}()

    # ? Perform discretization base on type of nonlinear terms ? #
    for i in m.disc_vars
        point = point_vec[i]                # Original Variable
        λCnt = length(discretization[i])

        # Built-in method based-on variable type
        if m.var_type[i] == :Cont
            i in processed && continue
            point = correct_point(m, discretization[i], point, i)
            for j in 1:(λCnt-1)
                if point >= discretization[i][j] && point <= discretization[i][j+1]  # Locating the right location
                    @show m.suggest_ratio
                    if m.suggest_ratio
                        haskey(m.disc_suggest,i) ? ratio = m.disc_suggest[i] : ratio = 8
                    end
                    radius = calculate_radius(discretization[i], j, ratio)
                    insert_partition(m, i, j, point, radius, discretization[i])
                    push!(processed, i)
                    break
                end
            end
        elseif m.var_type[i] == :Bin # This should never happen
            warn("Binary variable in m.disc_vars. Check out what is wrong...")
            continue  # No partition should be added to binary variable unless user specified
        elseif m.var_type[i] == :Int
            i in processed && continue
            point = discover_int_point(m, i, discretization[i], point)
            point == nothing && continue
            for j in 1:λCnt
                if point >= discretization[i][j] && point <= discretization[i][j+1]
                    insert_partition(m, i, j, point, 0.5, discretization[i])
                    push!(processed, i)
                    break
                end
            end
        else
            error("Unexpected variable types during injecting partitions")
        end
    end

    return discretization
end

function correct_point(m::PODNonlinearModel, partvec::Vector, point::Float64, var::Int)

    if point < partvec[1] - m.tol || point > partvec[end] + m.tol
        warn("VAR$(var) SOL=$(point) out of discretization [$(partvec[1]),$(partvec[end])]. Taking middle point...")
        return 0.5*(partvec[1] + partvec[end]) # Should choose the longest range
    end

    isapprox(point, partvec[1];atol=m.tol) && return partvec[1]
    isapprox(point, partvec[end];atol=m.tol) && return partvec[end]

    return point
end

function calculate_radius(partvec::Vector, part::Int, ratio::Any)

    lb_local = partvec[part]
    ub_local = partvec[part+1]

    distance = ub_local - lb_local
    if isa(ratio, Float64) || isa(ratio, Int)
        radius = distance / ratio
    elseif isa(ratio, Function)
        radius = distance / ratio(m)
    else
        error("Undetermined discretization ratio")
    end

    return radius
end

function insert_partition(m::PODNonlinearModel, var::Int, partidx::Int, point::Number, radius::Float64, partvec::Vector)

    abstol, reltol = m.disc_abs_width_tol, m.disc_rel_width_tol

    lb_local, ub_local = partvec[partidx], partvec[partidx+1]
    ub_touch, lb_touch = true, true
    lb_new, ub_new = max(point - radius, lb_local), min(point + radius, ub_local)

    if ub_new < ub_local && !isapprox(ub_new, ub_local; atol=abstol) && abs(ub_new-ub_local)/(1e-8+abs(ub_local)) > reltol # Insert new UB-based partition
        insert!(partvec, partidx+1, ub_new)
        ub_touch = false
    end

    if lb_new > lb_local && !isapprox(lb_new, lb_local; atol=abstol) && abs(lb_new-lb_local)/(1e-8+abs(lb_local)) > reltol # Insert new LB-based partition
        insert!(partvec, partidx+1, lb_new)
        lb_touch = false
    end

    if (ub_touch && lb_touch) || check_solution_history(m, var)
        distvec = [(j, partvec[j+1]-partvec[j]) for j in 1:length(partvec)-1]
        sort!(distvec, by=x->x[2])
        point_orig = point
        pos = distvec[end][1]
        lb_local = partvec[pos]
        ub_local = partvec[pos+1]
        isapprox(lb_local, ub_local;atol=m.tol) && return
        chunk = (ub_local - lb_local) / m.disc_divert_chunks
        point = lb_local + (ub_local - lb_local) / m.disc_divert_chunks
        for i in 2:m.disc_divert_chunks
            insert!(partvec, pos+1, lb_local + chunk * (m.disc_divert_chunks-(i-1)))
        end
        (m.loglevel > 99) && println("[DEBUG] !D! VAR$(var): SOL=$(round(point_orig,4))=>$(point) |$(round(lb_local,4)) | $(m.disc_divert_chunks) SEGMENTS | $(round(ub_local,4))|")
    else
        (m.loglevel > 99) && println("[DEBUG] VAR$(var): SOL=$(round(point,4)) RADIUS=$(radius), PARTITIONS=$(length(partvec)-1) |$(round(lb_local,4)) |$(round(lb_new,6)) <- * -> $(round(ub_new,6))| $(round(ub_local,4))|")
    end

    return
end

function add_uniform_partition(m::PODNonlinearModel; kwargs...)

    options = Dict(kwargs)
    haskey(options, :use_disc) ? discretization = options[:use_disc] : discretization = m.discretization

    for i in m.disc_vars  # Only construct when discretized
        lb_local = discretization[i][1]
        ub_local = discretization[i][end]
        distance = ub_local - lb_local
        chunk = distance / ((m.logs[:n_iter]+1)*m.disc_uniform_rate)
        discretization[i] = [lb_local+chunk*(j-1) for j in 1:(m.logs[:n_iter]+1)*m.disc_uniform_rate]
        push!(discretization[i], ub_local)   # Safety Scheme
        (m.loglevel > 99) && println("[DEBUG] VAR$(i): RATE=$(m.disc_uniform_rate), PARTITIONS=$(length(discretization[i]))  |$(round(lb_local,4)) | $(m.disc_uniform_rate*(1+m.logs[:n_iter])) SEGMENTS | $(round(ub_local,4))|")
    end

    return discretization
end

function update_disc_ratio(m::PODNonlinearModel, presolve=false)

    m.logs[:n_iter] > 2 && return m.disc_ratio # Stop branching after the second iterations

    ratio_pool = [8:2:20;]  # Built-in try range
    convertor = Dict(:Max=>:<, :Min=>:>)
    revconvertor = Dict(:Max=>:>, :Min=>:<)

    incumb_ratio = ratio_pool[1]
    m.sense_orig == :Min ? incumb_res = -Inf : incumb_res = Inf
    res_collector = Float64[]

    for r in ratio_pool
        st = time()
        if !isempty(m.best_sol)
            branch_disc = add_adaptive_partition(m, use_disc=m.discretization, branching=true, use_ratio=r, use_solution=m.best_sol)
        else
            branch_disc = add_adaptive_partition(m, use_disc=m.discretization, branching=true, use_ratio=r)
        end
        create_bounding_mip(m, use_disc=branch_disc)
        res = disc_branch_solve(m)
        push!(res_collector, res)
        if eval(convertor[m.sense_orig])(res, incumb_res) # && abs(abs(collector[end]-res)/collector[end]) > 1e-1  # %1 of difference
            incumb_res = res
            incumb_ratio = r
        end
        et = time() - st
        if et > 300  # 5 minutes limitation
            println("Expensive disc branching pass... Fixed at 8")
            return 8
        end
        m.loglevel > 0 && println("BRANCH RATIO = $(r), METRIC = $(res) || TIME = $(time()-st)")
    end

    if std(res_collector) >= 1e-2    # Detect if all solution are similar to each other
        m.loglevel > 0 && println("RATIO BRANCHING OFF due to solution variance test passed.")
        m.disc_ratio_branch = false # If an incumbent ratio is selected, then stop the branching scheme
    end

    if !isempty(m.best_sol)
        m.discretization = add_adaptive_partition(m, use_disc=m.discretization, branching=true, use_ratio=incumb_ratio, use_solution=m.best_sol)
    else
        m.discretization = add_adaptive_partition(m, use_disc=m.discretization, branching=true, use_ratio=incumb_ratio)
    end

    m.loglevel > 0 && println("INCUMB_RATIO = $(incumb_ratio)")

    return incumb_ratio
end

function disc_branch_solve(m::PODNonlinearModel)

    # ================= Solve Start ================ #
    update_mip_time_limit(m)
    start_bounding_solve = time()
    status = solve(m.model_mip, suppress_warnings=true)
    cputime_branch_bounding_solve = time() - start_bounding_solve
    m.logs[:total_time] += cputime_branch_bounding_solve
    m.logs[:time_left] = max(0.0, m.timeout - m.logs[:total_time])
    # ================= Solve End ================ #

    if status in [:Optimal, :Suboptimal, :UserLimit]
        return m.model_mip.objBound
    else
        warn("Unexpected solving condition $(status) during disc branching.")
    end

    # Safety scheme
    if m.sense_orig == :Min
        return -Inf
    else
        return Inf
    end
end
