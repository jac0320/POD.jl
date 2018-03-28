"""
    pick_disc_vars(m::PODNonlinearModel)

This function helps pick the variables for discretization. The method chosen depends on user-inputs.
In case when `indices::Int` is provided, the method is chosen as built-in method. Currently,
there exist two built-in method:

    * `max-cover(m.disc_var_pick=0, default)`: pick all variables involved in the non-linear term for discretization
    * `min-vertex-cover(m.disc_var_pick=1)`: pick a minimum vertex cover for variables involved in non-linear terms so that each non-linear term is at least convexified

For advance usage, `m.disc_var_pick` allows `::Function` inputs. User is required to perform flexible methods in choosing the non-linear variable.
For more information, read more details at [Hacking Solver](@ref).

"""
function pick_disc_vars(m::PODNonlinearModel)

    if isa(m.disc_var_pick, Function)
        eval(m.disc_var_pick)(m)
        length(m.disc_vars) == 0 && length(m.nonconvex_terms) > 0 && error("[USER FUNCTION] must select at least one variable to perform discretization for convexificiation purpose")
    elseif isa(m.disc_var_pick, Int)
        if m.disc_var_pick == 0
            ncvar_collect_nodes(m)
        elseif m.disc_var_pick == 1
            min_vertex_cover(m)
        elseif m.disc_var_pick == 2
            (length(m.candidate_disc_vars) > 15) ? min_vertex_cover(m) : ncvar_collect_nodes(m)
        elseif m.disc_var_pick == 3 # Initial
            (length(m.candidate_disc_vars) > 15) ? min_vertex_cover(m) : ncvar_collect_nodes(m)
        else
            error("Unsupported default indicator for picking variables for discretization")
        end
    else
        error("Input for parameter :disc_var_pick is illegal. Should be either a Int for default methods indexes or functional inputs.")
    end

    return
end

"""
    Collects distance of each variable's bounding solution to best feasible solution and select the ones that is the furthest
    Currently don't support recursively convexification
"""
function update_disc_cont_var(m::PODNonlinearModel)

    length(m.candidate_disc_vars) <= 15 && return   # Algorithm Separation Point

    # If no feasible solution found, do NOT update
    if m.status[:feasible_solution] != :Detected
        println("no feasible solution detected. No update disc var selection.")
        return
    end

	var_idxs = copy(m.candidate_disc_vars)
    var_diffs = Vector{Float64}(m.num_var_orig+length(keys(m.linear_terms))+length(keys(m.nonconvex_terms)))

    for i in 1:m.num_var_orig       # Original Variables
        var_diffs[i] = abs(m.best_sol[i]-m.best_bound_sol[i])
    end

    for i in 1:length(keys(m.linear_terms))
        for j in keys(m.linear_terms)   # sequential evaluation to avoid dependency issue
            if m.linear_terms[j][:id] == i
                var_diffs[m.linear_terms[j][:lifted_var_ref].args[2]] = m.linear_terms[j][:evaluator](m.linear_terms[j], var_diffs)
            end
        end
    end

    for i in 1:length(keys(m.nonconvex_terms))
        for j in keys(m.nonconvex_terms)    # sequential evaluation to avoid dependency issue
            if m.nonconvex_terms[j][:id] == i
                var_diffs[m.nonconvex_terms[j][:lifted_var_ref].args[2]] = m.nonconvex_terms[j][:evaluator](m.nonconvex_terms[j], var_diffs)
            end
        end
    end

    distance = Dict(zip(var_idxs,var_diffs))
    weighted_min_vertex_cover(m, distance)

    (m.loglevel > 100) && println("updated partition var selection => $(m.disc_vars)")
    return
end

function update_disc_int_var(m::PODNonlinearModel)

    length(m.candidate_disc_vars) <= 15 && return   # Algorithm Separation Point

    # If all discretized integer solution have been fully discovered, then include all other integer variables
    checker = [is_fully_discovered_integer(i, m.discretization[i]) for i in m.disc_vars if m.var_type[i] == :Int]
    if prod(checker)
        for i in m.candidate_disc_vars
            if m.var_type[i] == :Int && !(i in m.disc_vars)
                m.loglevel > 99 && println("PUMPing Integer VAR$(i) into discretization var set.")
                push!(m.disc_vars, i)
            end
        end
    end

    return
end

"""

    ncvar_collect_nodes(m:PODNonlinearModel)

A built-in method for selecting variables for discretization. It selects all variables in the nonlinear terms.

"""
function ncvar_collect_nodes(m::PODNonlinearModel;getoutput=false)

    # Pick variables that is bound width more than tolerance length
    if getoutput
        return [i for i in m.candidate_disc_vars]
    else
        m.disc_vars = [i for i in m.candidate_disc_vars]
        m.num_var_disc_mip = length(m.disc_vars)
    end

    return
end

"""
    Reconsideration required
"""
function ncvar_collect_arcs(m::PODNonlinearModel, nodes::Vector)

    arcs = Set()

    for k in keys(m.nonconvex_terms)
        if m.nonconvex_terms[k][:nonlinear_type] == :BILINEAR
            arc = [i.args[2] for i in k]
            length(arc) == 2 && push!(arcs, sort(arc))
        elseif m.nonconvex_terms[k][:nonlinear_type] == :MONOMIAL
            @assert isa(m.nonconvex_terms[k][:var_idxs][1], Int)
            varidx = m.nonconvex_terms[k][:var_idxs][1]
            push!(arcs, [varidx, varidx;])
        elseif m.nonconvex_terms[k][:nonlinear_type] == :MULTILINEAR
            varidxs = m.nonconvex_terms[k][:var_idxs]
            for i in 1:length(varidxs)
                for j in 1:length(varidxs)
                    if i != j
                        push!(arcs, sort([varidxs[i], varidxs[j];]))
                    end
                end
            end
            if length(varidxs) == 1
                push!(arcs, sort([varidxs[1], varidxs[1];]))
            end
        elseif m.nonconvex_terms[k][:nonlinear_type] == :INTLIN
            var_idxs = copy(m.nonconvex_terms[k][:var_idxs])
            push!(arcs, sort(var_idxs))
        elseif m.nonconvex_terms[k][:nonlinear_type] == :INTPROD
            var_idxs = m.nonconvex_terms[k][:var_idxs]
            for i in 1:length(var_idxs)
                for j in 1:length(var_idxs)
                    i != j && push!(arcs, sort([var_idxs[i], var_idxs[j];]))
                end
            end
        elseif m.nonconvex_terms[k][:nonlinear_type] in [:cos, :sin]
            @assert length(m.nonconvex_terms[k][:var_idxs]) == 1
            var_idx = m.nonconvex_terms[k][:var_idxs][1]
            push!(arcs, [var_idx, var_idx;])
        elseif m.nonconvex_terms[k][:nonlinear_type] in [:BININT, :BINLIN, :BINPROD]
            continue
        else
            error("[EXCEPTION] Unexpected nonlinear term when building discvar graph.")
        end
    end

    return arcs
end

"""
    Tell what would be the variable type of a lifted term.
    This function is with limited functionality
    @docstring TODO
"""
function resolve_lifted_var_type(var_types::Vector{Symbol}, operator::Symbol)

    if operator == :+
        detector = [i in [:Bin, :Int] ? true : false for i in var_types]
        if length(detector) == 1 && detector[1] # Special case
            if :Bin in var_types
                return :Bin
            else
                return :Int
            end
        end
        prod(detector) && return :Int
        # o/w continous variables
    elseif operator == :*
        detector = [i == :Bin ? true : false for i in var_types]
        prod(detector) && return :Bin
        detector = [i in [:Bin, :Int] ? true : false for i in var_types]
        prod(detector) && return :Int
        # o/w continous variables
    end

    return :Cont
end

"""
    TODO can be improved
"""
function build_discvar_graph(m::PODNonlinearModel)

    # Collect the information of nonlinear terms in terms of arcs and nodes
    nodes = ncvar_collect_nodes(m, getoutput=true)
    arcs = ncvar_collect_arcs(m, nodes)

    # Collect integer variables
    for i in 1:m.num_var_orig
        if !(i in nodes) && m.var_type[i] == :Int
            push!(nodes, i)
            push!(arcs, [i,i;])
        end
    end

    nodes = collect(nodes)
    arcs = collect(arcs)

    return nodes, arcs
end

function min_vertex_cover(m::PODNonlinearModel)

    nodes, arcs = build_discvar_graph(m)

    # Set up minimum vertex cover problem
    update_mip_time_limit(m, timelimit=60.0)  # Set a timer to avoid waste of time in proving optimality
    minvertex = Model(solver=m.mip_solver)
    @variable(minvertex, x[nodes], Bin)
    @constraint(minvertex, [a in arcs], x[a[1]] + x[a[2]] >= 1)
    @objective(minvertex, Min, sum(x))

    status = solve(minvertex, suppress_warnings=true)
    xVal = getvalue(x)

    # Collecting required information
    m.num_var_disc_mip = Int(sum(xVal))
    m.disc_vars = [i for i in nodes if xVal[i] > m.tol && abs(m.u_var_tight[i]-m.l_var_tight[i]) >= m.tol]

    return
end

function weighted_min_vertex_cover(m::PODNonlinearModel, distance::Dict)

    # Collect the graph information
    nodes, arcs = build_discvar_graph(m)

    # A little bit redundency before
    disvec = [distance[i] for i in keys(distance) if i in m.candidate_disc_vars]
    disvec = abs.(disvec[disvec .> 0.0])
    isempty(disvec) ? heavy = 1.0 : heavy = 1/minimum(disvec)
    weights = Dict()
    for i in m.candidate_disc_vars
        isapprox(distance[i], 0.0; atol=1e-6) ? weights[i] = heavy : (weights[i]=(1/distance[i]))
        (m.loglevel > 100) && println("VAR$(i) WEIGHT -> $(weights[i]) ||| DISTANCE -> $(distance[i])")
    end

    # Set up minimum vertex cover problem
    update_mip_time_limit(m, timelimit=60.0)  # Set a timer to avoid waste of time in proving optimality
    minvertex = Model(solver=m.mip_solver)
    @variable(minvertex, x[nodes], Bin)
    for arc in arcs
        @constraint(minvertex, x[arc[1]] + x[arc[2]] >= 1)
    end
    @objective(minvertex, Min, sum(weights[i]*x[i] for i in nodes))

    # Solve the minimum vertex cover
    status = solve(minvertex, suppress_warnings=true)

    xVal = getvalue(x)
    m.num_var_disc_mip = Int(sum(xVal))
    m.disc_vars = [i for i in nodes if xVal[i] > 0 && abs(m.u_var_tight[i]-m.l_var_tight[i]) >= m.tol]
    m.loglevel >= 99 && println("UPDATED DISC-VAR COUNT = $(length(m.disc_vars)) : $(m.disc_vars)")

    return
end

"""
    Follow the definition of terms to calculate the value of lifted terms
"""
function resolve_lifted_var_value(m::PODNonlinearModel, sol_vec::Array)

    @assert length(sol_vec) == m.num_var_orig
    sol_vec = [sol_vec; fill(NaN, m.num_var_linear_mip+m.num_var_nonlinear_mip)]

    for i in 1:length(m.term_seq)
        k = m.term_seq[i]
        if haskey(m.nonconvex_terms, k)
            lvar_idx = m.nonconvex_terms[k][:y_idx]
            sol_vec[lvar_idx] = m.nonconvex_terms[k][:evaluator](m.nonconvex_terms[k], sol_vec)
        elseif haskey(m.linear_terms, k)
            lvar_idx = m.linear_terms[k][:y_idx]
            sol_vec[lvar_idx] = m.linear_terms[k][:evaluator](m.linear_terms[k], sol_vec)
        else
            error("[RARE] Found homeless term key $(k) during bound resolution.")
        end
    end

    return sol_vec
end
