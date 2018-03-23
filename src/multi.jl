function amp_post_convhull(m::PODNonlinearModel; use_disc=nothing, warmstart=true, cuts=true)

    use_disc == nothing ? d = m.discretization : d = use_disc

    # Variable holders
    λ = Dict()  # Extreme points and multipliers
    α = Dict()  # Partitioning Variables
    β = Dict()  # Lifted variables for exact formulation

    # Convexification Treatment for Complex Non-Convex Terms
    for k in keys(m.nonconvex_terms)
        nl_type = m.nonconvex_terms[k][:nonlinear_type]
        if ((nl_type == :MULTILINEAR) || (nl_type == :BILINEAR)) && (m.nonconvex_terms[k][:convexified] == false)
            λ, α = amp_convexify_multilinear(m, k, λ, α, d)
        elseif nl_type == :MONOMIAL && !m.nonconvex_terms[k][:convexified]
            λ, α = amp_convexify_monomial(m, k, λ, α, d)
        elseif nl_type == :BINLIN && !m.nonconvex_terms[k][:convexified]
            β = amp_convexify_binlin(m, k, β)
        elseif nl_type == :BINPROD && !m.nonconvex_terms[k][:convexified]
            β = amp_convexify_binprod(m, k, β)
        elseif nl_type == :BININT && !m.nonconvex_terms[k][:convexified]
            β = amp_convexify_binint(m, k, β)
        elseif nl_type == :INTPROD && !m.nonconvex_terms[k][:convexified]
            m.int_enable || error("Integer features are OFF. No support for INTPROD at this condition.")
            λ, α = amp_convexify_intprod(m, k, λ, α, d)
        elseif nl_type == :INTLIN && !m.nonconvex_terms[k][:convexified]
            m.int_enable || error("Integer features are OFF. No support for INTLIN at this condition.")
            λ, α = amp_convexify_intlin(m, k, λ, α, d)
        elseif nl_type in [:sin, :cos] && !m.nonconvex_terms[k][:convexified]
            λ, α = amp_convexify_sincos(m, k, λ, α, d)
        end
    end

    # Experimental Code for Integer Problems
    if m.int_enable
        for i in 1:m.num_var_orig
            m.var_type[i] == :Int && amp_convexify_integer(m, i, λ, α, d)
        end
    end

    # Experimental code for Warm starting
    warmstart & m.convhull_warmstart && !m.convhull_ebd && !isempty(m.best_bound_sol) && amp_warmstart_α(m, α)
    cuts && amp_post_ac_cuts(m, α)

    return
end

function amp_convexify_integer(m::PODNonlinearModel, idx::Int, λ::Dict, α::Dict, d::Dict)

    # Check if VAR idx exist in any Complex Nonlinear Terms
    # If so, VAR idx has already been discretized/relaxed, in its linear term, the value will be relaxed as well
    for k in keys(m.nonconvex_terms)
        idx in m.nonconvex_terms[k][:var_idxs] && return
    end

    # Easy way to check whether integer variable idx has been discretized already
    haskey(α, idx) && return

    # println("Resolving SOLE integer VAR$(idx) ..")
    # Main Convexification Steps
    indices, dim, ext_cnt = amp_convhull_prepare(m, d, idx)
    λ = amp_convhull_λ(m, idx, λ, ext_cnt, dim)
    λ = populate_convhull_extreme_values(m, d, idx, λ)
    α = amp_convhull_α(m, indices, α, dim, d)
    amp_post_convhull_constrs(m, λ, α, indices, dim, ext_cnt, d)

    # Additional Relaxation : Safety
    setcategory(Variable(m.model_mip, idx), :Cont)

    return
end

function amp_convexify_sincos(m::PODNonlinearModel, k::Any, λ::Dict, α::Dict, d::Dict)

    # Be careful with the single dimension here
    # TODO build test cases for this problem
    m.nonconvex_terms[k][:convexified] = true  # Bookeeping the convexified terms
    trifunc_index, dim, ext_cnt = amp_convhull_prepare(m, d, k)   # convert key to easy read mode
    # haskey(λ, trifunc_index) && error("Overlapping but not sure if will cause bug...")
    λ = amp_convhull_λ(m, k, trifunc_index, λ, ext_cnt, dim)
    λ = populate_sincos_extreme_values(m, d, trifunc_index, m.nonconvex_terms[k][:nonlinear_type], λ)
    α = amp_convhull_α(m, trifunc_index, α, dim, d)
    amp_post_constrs_sincos(m, λ, α, trifunc_index, dim, ext_cnt, d, m.nonconvex_terms[k][:nonlinear_type])

    # post dedicated cutting planes

    return λ, α
end

function amp_convexify_multilinear(m::PODNonlinearModel, k::Any, λ::Dict, α::Dict, discretization::Dict)

    m.nonconvex_terms[k][:convexified] = true  # Bookeeping the convexified terms

    ml_indices, dim, extreme_point_cnt = amp_convhull_prepare(m, discretization, k)   # convert key to easy read mode
    λ = amp_convhull_λ(m, k, ml_indices, λ, extreme_point_cnt, dim)
    λ = populate_convhull_extreme_values(m, discretization, ml_indices, λ, dim, ones(Int,length(dim)))
    α = amp_convhull_α(m, ml_indices, α, dim, discretization)
    amp_post_convhull_constrs(m, λ, α, ml_indices, dim, extreme_point_cnt, discretization)

    return λ, α
end

function amp_convexify_intprod(m::PODNonlinearModel, k::Any, λ::Dict, α::Dict, d::Dict)

    m.nonconvex_terms[k][:convexified] = true  # Bookeeping the convexified terms

    intprod_indices, dim, ext_cnt = amp_convhull_prepare(m, d, k)   # convert key to easy read mode
    λ = amp_convhull_λ(m, k, intprod_indices, λ, ext_cnt, dim)
    λ = populate_convhull_extreme_values(m, d, intprod_indices, λ, dim, ones(Int,length(dim)))
    α = amp_convhull_α(m, intprod_indices, α, dim, d)
    amp_post_convhull_constrs(m, λ, α, intprod_indices, dim, ext_cnt, d)
    amp_post_special_λ_ub(m, intprod_indices, dim, λ, d)

    return λ, α
end

function amp_convexify_intlin(m::PODNonlinearModel, k::Any, λ::Dict, α::Dict, d::Dict)

    m.nonconvex_terms[k][:convexified] = true  # Bookeeping the convexified terms

    intlin_indices, dim, ext_cnt = amp_convhull_prepare(m, d, k)   # convert key to easy read mode
    @assert length(intlin_indices) == 2

    λ = amp_convhull_λ(m, k, intlin_indices, λ, ext_cnt, dim)
    λ = populate_convhull_extreme_values(m, d, intlin_indices, λ, dim, ones(Int,length(dim)))
    α = amp_convhull_α(m, intlin_indices, α, dim, d)
    amp_post_convhull_constrs(m, λ, α, intlin_indices, dim, ext_cnt, d)
    amp_post_special_λ_lb(m, intlin_indices, dim, λ, α, d)

    return λ, α
end

function amp_convexify_monomial(m::PODNonlinearModel, k::Any, λ::Dict, α::Dict, discretization::Dict)

    m.nonconvex_terms[k][:convexified] = true  # Bookeeping the convexified terms

    monomial_index, dim, extreme_point_cnt = amp_convhull_prepare(m, discretization, k, monomial=true)
    λ = amp_convhull_λ(m, k, monomial_index, λ, extreme_point_cnt, dim)
    λ = populate_convhull_extreme_values(m, discretization, monomial_index, λ, 2)
    α = amp_convhull_α(m, [monomial_index], α, dim, discretization)
    amp_post_convhull_constrs(m, λ, α, monomial_index, dim, discretization)

    return λ, α
end

function amp_convexify_binlin(m::PODNonlinearModel, k::Any, β::Dict)

    m.nonconvex_terms[k][:convexified] = true  # Bookeeping the convexified terms

    @assert length(m.nonconvex_terms[k][:var_idxs]) == 2

    lift_idx = m.nonconvex_terms[k][:y_idx]

    if haskey(β, lift_idx)
        return β
    else
        β[lift_idx] = Variable(m.model_mip, lift_idx)
    end

    bin_idx = [i for i in m.nonconvex_terms[k][:var_idxs] if m.var_type[i] == :Bin]
    cont_idx = [i for i in m.nonconvex_terms[k][:var_idxs] if m.var_type[i] != :Bin]

    @assert length(bin_idx) == length(cont_idx) == 1

    bin_idx = bin_idx[1]
    cont_idx = cont_idx[1]

    mccormick_binlin(m.model_mip, Variable(m.model_mip, lift_idx),
        Variable(m.model_mip, bin_idx), Variable(m.model_mip, cont_idx),
        m.l_var_tight[cont_idx], m.u_var_tight[cont_idx])

    return β
end

amp_convexify_binint(m::PODNonlinearModel, k::Any, β::Dict) = amp_convexify_binlin(m, k, β)

function amp_convexify_binprod(m::PODNonlinearModel, k::Any, β::Dict)

    m.nonconvex_terms[k][:convexified] = true  # Bookeeping the convexified terms

    lift_idx = m.nonconvex_terms[k][:y_idx]
    if haskey(β, lift_idx)
        return β    # Already constructed
    else
        β[lift_idx] = Variable(m.model_mip, lift_idx)
    end

    z = Variable(m.model_mip, m.nonconvex_terms[k][:y_idx])
    x = [Variable(m.model_mip, i) for i in m.nonconvex_terms[k][:var_idxs]]
    for i in x
        @constraint(m.model_mip, z <= i)
    end
    @constraint(m.model_mip, z >= sum(x) - (length(x)-1))

    return β
end

"""
    Method for general nonlinear terms
"""
function amp_convhull_prepare(m::PODNonlinearModel, d::Dict, nonlinear_key::Any; monomial=false)

    counted_var = []                # Keep both vector and set for collection sake
    id = Set()                      # Coverting the nonlinear indices into a set

    if isa(nonlinear_key, Vector)
        for var in nonlinear_key        # This output regulates the sequence of how composing variable should be arranged
            @assert isa(var, Expr)
            m.var_type[var.args[2]] in [:Cont, :Int] && push!(id, var.args[2])
            m.var_type[var.args[2]] in [:Cont, :Int] && push!(counted_var, var.args[2])
        end
    elseif isa(nonlinear_key, Dict)
        for var in nonlinear_key[:vars]
            @assert isa(var, Int)
            m.var_type[var] in [:Cont, :Int] && push!(id, var)
            m.var_type[var] in [:Cont, :Int] && push!(counted_var, var)
        end
    end

    if length(id) < length(counted_var) # Got repeating terms, now the sequence matters
        id = []
        for var in nonlinear_key
            m.var_type[var.args[2]] in [:Cont, :Int] && push!(id, var.args[2])
        end
    end

    dim = [length(d[i]) for i in id]

    monomial && return id[1], tuple(dim[1]), dim[1]   # One less dimension is required
    return id, tuple([i for i in dim]...), prod(dim)
end

"""
    Method for integers
"""
function amp_convhull_prepare(m::PODNonlinearModel, d::Dict, idx::Int)
    return [idx], tuple(length(d[idx])), length(d[idx])
end

"""
    Method for general nonlinear terms
"""
function amp_convhull_λ(m::PODNonlinearModel, nonlinear_key::Any, indices::Any, λ::Dict, ext_cnt::Int, dim::Tuple)

    y_idx = m.nonconvex_terms[nonlinear_key][:y_idx]

    @assert !(y_idx in keys(λ))
    λ[indices] = Dict(:dim=>dim,
                     :lifted_var_idx=>y_idx,
                     :indices=>reshape([1:ext_cnt;], dim),
                     :vars=>@variable(m.model_mip, [1:ext_cnt], lowerbound=0, basename="L$(y_idx)"),
                     :vals=>ones(dim))

    return λ
end

"""
    Method for integers
"""
function amp_convhull_λ(m::PODNonlinearModel, idx::Int, λ::Dict, ext_cnt::Int, dim::Tuple)

    @assert !([idx] in keys(λ))
    λ[[idx]] = Dict(:dim=>dim,
                      :lifted_var_idx=>idx,
                      :indices=>[1:ext_cnt;],
                      :vars=>@variable(m.model_mip, [1:ext_cnt], lowerbound=0, basename="L$(idx)"),
                      :vals=>ones(dim))

    return λ
end

"""
    Method for integer variables
"""
function populate_convhull_extreme_values(m::PODNonlinearModel, d::Dict, int_index::Int, λ::Dict)
    λ[[int_index]][:vals] = [d[int_index][i] for i in 1:length(d[int_index])]
    return λ
end

"""
    Method for power terms
"""
function populate_convhull_extreme_values(m::PODNonlinearModel, d::Dict, mono_idx::Int, λ::Dict, p::Int)
    λ[mono_idx][:vals] = [d[mono_idx][i]^p for i in 1:length(d[mono_idx])]
    return λ
end

"""
    Method for sin/cos... and potential some more terms
"""
function populate_sincos_extreme_values(m::PODNonlinearModel, d::Dict, λ_Key::Any, operator::Symbol, λ::Dict)

    @assert isa(λ_Key, Vector) || isa(λ_Key, Set)
    @assert length(λ_Key) == 1

    var = pop!(λ_Key)
    push!(λ_Key, var)
    λCnt = length(d[var])
    λ[λ_Key][:vals] = [eval(operator)(d[var][i]) for i in 1:λCnt]

    return λ
end

"""
    Method for regular muiltilinear terms
"""
function populate_convhull_extreme_values(m::PODNonlinearModel, discretization::Dict, indices::Any, λ::Dict, dim::Tuple, locator::Array, level::Int=1)

    if level > length(dim)
        @assert length(indices) == length(dim)
        @assert length(indices) == length(locator)
        val = 1.0
        k = 0
        for i in indices
            k += 1
            val *= discretization[i][locator[k]]                      # Calculate extreme point z-value
        end
        λ[indices][:vals][CartesianIndex(tuple([i for i in locator]...))] = val  # Value assignment
        return λ                                                         # finished with last dimension
    else
        for i in 1:dim[level]
            locator[level] = i
            λ = populate_convhull_extreme_values(m, discretization, indices, λ, dim, locator, level+1)
        end
    end

    return λ
end

"""
    General Method for all term
"""
function amp_convhull_α(m::PODNonlinearModel, indices::Any, α::Dict, dim::Tuple, discretization::Dict)

    for i in indices
        if !(i in keys(α))
            lambda_cnt = length(discretization[i])
            partition_cnt = length(discretization[i]) - 1
            if m.convhull_ebd && m.var_type[i] == :Cont && partition_cnt > 2
                αCnt = Int(ceil(log(2,partition_cnt)))
                α[i] = @variable(m.model_mip, [1:αCnt], Bin, basename=string("YL",i))
            else
                α[i] = @variable(m.model_mip, [1:partition_cnt], Bin, basename="A$(i)")
                @constraint(m.model_mip, sum(α[i]) == 1)
                @constraint(m.model_mip, Variable(m.model_mip, i) >= sum(α[i][j]*discretization[i][j] for j in 1:lambda_cnt-1)) # Add x = f(α) for regulating the domains
                @constraint(m.model_mip, Variable(m.model_mip, i) <= sum(α[i][j-1]*discretization[i][j] for j in 2:lambda_cnt))
            end
        end
    end

    return α
end

amp_convhull_α(m::PODNonlinearModel, idx::Int, α::Dict, dim, d::Dict) = amp_convhull_α(m, [idx], α, dim, d)

function amp_no_good_cut_α(m::PODNonlinearModel, α::Dict)

    println("Global Incumbent solution objective = $(m.best_obj)")

    for i in 1:m.bound_sol_pool[:cnt]
        (m.bound_sol_pool[:stat][i] == :Cutoff) && (m.bound_sol_pool[:stat][i] = :Alive)
        if m.best_obj < m.bound_sol_pool[:obj][i] && m.bound_sol_pool[:stat][i] == :Alive
            no_good_idxs = keys(m.bound_sol_pool[:disc][i])
            no_good_size = length(no_good_idxs) - 1
            @constraint(m.model_mip, sum(α[v][m.bound_sol_pool[:disc][i][v]] for v in no_good_idxs) <= no_good_size)
            m.loglevel > 0 && println("!! GLOBAL cuts off POOL_SOL-$(i) POOL_OBJ=$(m.bound_sol_pool[:obj][i])!")
            m.bound_sol_pool[:stat][i] = :Cutoff
        end
    end

    return
end

function amp_warmstart_α(m::PODNonlinearModel, α::Dict)

    d = m.discretization

    if m.bound_sol_pool[:cnt] >= 2 # can only warm-start the problem when pool is large enough
        ws_idx = -1
        m.sense_orig == :Min ? ws_obj = Inf : ws_obj = -Inf
        comp_opr = Dict(:Min=>:<, :Max=>:>)

        # Search for the pool for incumbent warm starter
        for i in 1:m.bound_sol_pool[:cnt]
            m.bound_sol_pool[:stat][i] == :Warmstarter && (m.bound_sol_pool[:stat][i] = :Alive)   # reset the status if not dead
            if m.bound_sol_pool[:stat][i] != :Dead && eval(comp_opr[m.sense_orig])(m.bound_sol_pool[:obj][i], ws_obj)
                ws_idx = i
                ws_obj = m.bound_sol_pool[:obj][i]
            end
        end

        if ws_idx > 0 # If a warm starter is found
            for v in m.bound_sol_pool[:vars]
                partition_cnt = length(d[v])-1
                active_j = get_active_partition_idx(d, m.bound_sol_pool[:sol][ws_idx][v], v)
                for j = 1:partition_cnt
                    j == active_j ? setvalue(α[v][j], 1.0) : setvalue(α[v][j], 0.0)
                end
            end
            m.bound_sol_pool[:stat][ws_idx] = :Warmstarter
            m.loglevel > 0 && println("!! WARM START bounding MIP using POOL SOL $(ws_idx) OBJ=$(m.bound_sol_pool[:obj][ws_idx])")
        end
    end

    return
end

function amp_post_ac_cuts(m::PODNonlinearModel, α::Dict)

    for c in keys(m.arc_consistency_cuts)
        D = length(c)
        @constraint(m.model_mip, sum(α[i[1]][i[2]] for i in c) <= D - 1)
    end

    return
end

"""
    Method for general multilinear terms with/without integer variables
"""
function amp_post_convhull_constrs(m::PODNonlinearModel, λ::Dict, α::Dict, indices::Any, dim::Tuple, ext_cnt::Int, d::Dict)

    # Adding λ constraints
    @constraint(m.model_mip, sum(λ[indices][:vars]) == 1)
    @constraint(m.model_mip, Variable(m.model_mip, λ[indices][:lifted_var_idx]) == dot(λ[indices][:vars], reshape(λ[indices][:vals], ext_cnt)))

    # Add links on each dimension
    for (cnt, i) in enumerate(indices)
        l_cnt = length(d[i])
        if m.var_type[i] == :Cont
            amp_post_inequalities_cont(m, d, λ, α, indices, dim, i, cnt)        # Add links between λ and α
        elseif m.var_type[i] == :Int
            amp_post_inequalities_int(m, d, λ, α, indices, dim, i, cnt)         # Add links between λ and α
        else
            error("EXCEPTION: unexpected variable type during integer related realxation")
        end
        sliced_indices = [collect_indices(λ[indices][:indices], cnt, [k], dim) for k in 1:l_cnt] # Add x = f(λ) for convex representation of x value
        @constraint(m.model_mip, Variable(m.model_mip, i) == sum(dot(repmat([d[i][k]],length(sliced_indices[k])), λ[indices][:vars][sliced_indices[k]]) for k in 1:l_cnt))
    end

    return
end

"""
    Method for sin/cos terms
"""
function amp_post_constrs_sincos(m::PODNonlinearModel, λ::Dict, α::Dict, indices::Any, dim::Tuple, ext_cnt::Int, d::Dict, opt::Symbol)

    # Some basic assertion as the function has limited scope
    @assert length(indices) == length(dim) == 1

    # Adding λ constraints
    @constraint(m.model_mip, sum(λ[indices][:vars]) == 1)

    # Special Treatment for λ constriants
    @constraint(m.model_mip, Variable(m.model_mip, λ[indices][:lifted_var_idx]) <= dot(λ[indices][:vars], ones(ext_cnt)))
    @constraint(m.model_mip, Variable(m.model_mip, λ[indices][:lifted_var_idx]) >= dot(λ[indices][:vars], -ones(ext_cnt)))

    # Upper & Lower Bound Regionss
    UBvec, LBvec = populate_bound_regions(λ, indices, d, opt)
    m.loglevel > 99 && println("Ave Regions Length = $(mean(UBvec-LBvec))")
    @constraint(m.model_mip, Variable(m.model_mip, λ[indices][:lifted_var_idx]) <= dot(λ[indices][:vars], UBvec))
    @constraint(m.model_mip, Variable(m.model_mip, λ[indices][:lifted_var_idx]) >= dot(λ[indices][:vars], LBvec))

    # Add links on each dimension
    for (cnt, i) in enumerate(indices)
        l_cnt = length(d[i])
        cnt <= 2 || error("EXCEPTION: two variable included in a sin/cos term, there is a bug.")
        m.var_type[i] == :Cont || error("EXCEPTION: unexpected variable type during integer related realxation")
        # Add links between λ and α
        amp_post_inequalities_cont(m, d, λ, α, indices, dim, i, cnt)
        # Add x = f(λ) for convex representation of x value
        sliced_indices = [collect_indices(λ[indices][:indices], cnt, [k], dim) for k in 1:l_cnt]
        @constraint(m.model_mip, Variable(m.model_mip, i) == sum(dot(repmat([d[i][k]],length(sliced_indices[k])), λ[indices][:vars][sliced_indices[k]]) for k in 1:l_cnt))
    end

    return
end

"""
    Method for power-2 term
"""
function amp_post_convhull_constrs(m::PODNonlinearModel, λ::Dict, α::Dict, monomial_idx::Int, dim::Tuple, discretization::Dict)

    partition_cnt = length(discretization[monomial_idx])-1
    lambda_cnt = length(discretization[monomial_idx])

    # Adding λ constraints
    @constraint(m.model_mip, sum(λ[monomial_idx][:vars]) == 1)
    @constraint(m.model_mip, Variable(m.model_mip, λ[monomial_idx][:lifted_var_idx]) <= dot(λ[monomial_idx][:vars], λ[monomial_idx][:vals]))
    @constraint(m.model_mip, Variable(m.model_mip, λ[monomial_idx][:lifted_var_idx]) >= Variable(m.model_mip, monomial_idx)^2)

    # Add SOS-2 Constraints with basic encoding
    if m.convhull_ebd && partition_cnt > 2
        ebd_map = embedding_map(lambda_cnt, m.convhull_ebd_encode, m.convhull_ebd_ibs)
        YCnt = Int(ebd_map[:L])
        @assert YCnt == length(α[monomial_idx])
        for i in 1:YCnt
            @constraint(m.model_mip, sum(λ[monomial_idx][:vars][collect(ebd_map[i])]) <= α[monomial_idx][i])
            @constraint(m.model_mip, sum(λ[monomial_idx][:vars][collect(ebd_map[i+YCnt])]) <= 1-α[monomial_idx][i])
        end
    else
        for i in 1:lambda_cnt
            if i == 1
                @constraint(m.model_mip, λ[monomial_idx][:vars][i] <= α[monomial_idx][i])
            elseif i == lambda_cnt
                @constraint(m.model_mip, λ[monomial_idx][:vars][i] <= α[monomial_idx][i-1])
            else
                @constraint(m.model_mip, λ[monomial_idx][:vars][i] <= α[monomial_idx][i-1] + α[monomial_idx][i])
            end
        end
        # Add x = f(α) for regulating the domains
        @constraint(m.model_mip, Variable(m.model_mip, monomial_idx) >= sum(α[monomial_idx][j]*discretization[monomial_idx][j] for j in 1:lambda_cnt-1))
        @constraint(m.model_mip, Variable(m.model_mip, monomial_idx) <= sum(α[monomial_idx][j-1]*discretization[monomial_idx][j] for j in 2:lambda_cnt))
    end

    # Add x = f(λ) for convex representation
    @constraint(m.model_mip, Variable(m.model_mip, monomial_idx) == dot(λ[monomial_idx][:vars], discretization[monomial_idx]))

    return
end

"""
    Method for regular multilinear terms
"""
function amp_post_inequalities_cont(m::PODNonlinearModel, discretization::Dict, λ::Dict, α::Dict, ml_indices::Any, dim::Tuple, var_ind::Int, cnt::Int)

    lambda_cnt = length(discretization[var_ind])
    partition_cnt = lambda_cnt - 1

    # Embedding formulation
    if m.convhull_formulation == "sos2" && m.convhull_ebd && partition_cnt > 2
        ebd_map = embedding_map(lambda_cnt, m.convhull_ebd_encode, m.convhull_ebd_ibs)
        YCnt = Int(ebd_map[:L])
        @assert YCnt == length(α[var_ind])
        for i in 1:YCnt
            p_sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, collect(ebd_map[i]), dim)
            n_sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, collect(ebd_map[i+YCnt]), dim)
            @constraint(m.model_mip, sum(λ[ml_indices][:vars][p_sliced_indices]) <= α[var_ind][i])
            @constraint(m.model_mip, sum(λ[ml_indices][:vars][n_sliced_indices]) <= 1-α[var_ind][i])
        end
        m.convhull_ebd_link && ebd_link_xα(m, α[var_ind], lambda_cnt, discretization[var_ind], ebd_map[:H_orig], var_ind)
        return
    end

    # SOS-2 Formulation
    if m.convhull_formulation == "sos2"
        for j in 1:lambda_cnt
            sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, [j], dim)
            if (j == 1)
                @constraint(m.model_mip, sum(λ[ml_indices][:vars][sliced_indices]) <= α[var_ind][j])
            elseif (j == lambda_cnt)
                @constraint(m.model_mip, sum(λ[ml_indices][:vars][sliced_indices]) <= α[var_ind][partition_cnt])
            else
                @constraint(m.model_mip, sum(λ[ml_indices][:vars][sliced_indices]) <= sum(α[var_ind][(j-1):j]))
            end
        end
        return
    elseif m.convhull_formulation == "facet"
        for j in 1:(partition_cnt-1) # Constraint cluster of α >= f(λ)
            sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, [1:j;], dim)
            @constraint(m.model_mip, sum(α[var_ind][1:j]) >= sum(λ[ml_indices][:vars][sliced_indices]))
        end
        for j in 1:(partition_cnt-1) # Constriant cluster of α <= f(λ)
            sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, [1:(j+1);], dim)
            @constraint(m.model_mip, sum(α[var_ind][1:j]) <= sum(λ[ml_indices][:vars][sliced_indices]))
        end
        return
    elseif m.convhull_formulation == "mini"
        for j in 1:min(partition_cnt, 1) # Constraint cluster of α >= f(λ)
            sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, [1;], dim)
            @constraint(m.model_mip, sum(α[var_ind][1:j]) >= sum(λ[ml_indices][:vars][sliced_indices]))
            sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, [lambda_cnt;], dim)
            @constraint(m.model_mip, sum(α[var_ind][(dim[cnt]-j):(dim[cnt]-1)]) >= sum(λ[ml_indices][:vars][sliced_indices]))
        end
        for j in 1:partition_cnt         # Constriant cluster of α <= f(λ)
            for i in 1:max(1, min(partition_cnt-j+1, 1)) # At least one
                sliced_indices = collect_indices(λ[ml_indices][:indices], cnt, [j:(j+i);], dim)
                @constraint(m.model_mip, sum(α[var_ind][j:(j+i-1)]) <= sum(λ[ml_indices][:vars][sliced_indices]))
            end
        end
        return
    else
        error("Must indicate a choice of convex hull formulation. ?(mini, sos2, facet)")
    end

    return
end

"""
    Method for multilinear terms with discrete variables
"""
function amp_post_inequalities_int(m::PODNonlinearModel, d::Dict, λ::Dict, α::Dict, indices::Any, dim::Tuple, var_ind::Int, cnt::Int)

    l_cnt = length(d[var_ind])
    p_cnt = l_cnt - 1

    # Embedding formulation
    m.convhull_ebd && warn("Embedding is currently not supported for multilinear terms with discrete variables", once=true)

    # SOS-2 Formulation
    if m.convhull_formulation == "sos2"
        for j in 1:l_cnt
            sliced_indices = collect_indices(λ[indices][:indices], cnt, [j], dim)
            if (j == 1)
                rate_r = amp_pick_ratevec(d[var_ind], j)
                @constraint(m.model_mip, sum(λ[indices][:vars][sliced_indices]) <= α[var_ind][j] * rate_r)
            elseif (j == l_cnt)
                rate_l = amp_pick_ratevec(d[var_ind], j)
                @constraint(m.model_mip, sum(λ[indices][:vars][sliced_indices]) <= α[var_ind][p_cnt] * rate_l)
            else
                rate_l, rate_r = amp_pick_ratevec(d[var_ind], j)
                @constraint(m.model_mip, sum(λ[indices][:vars][sliced_indices]) <= α[var_ind][j-1] * rate_l + α[var_ind][j] * rate_r)
            end
        end
    else
        error("Only SOS-2 formulation support convexification of terms with integer variables.")
    end

    return
end

function amp_pick_ratevec(partvec::Vector, i::Int)

    λCnt = length(partvec)

    if i == 1
        mod(partvec[i], 0.5) == 0.0 && partvec[i+1] - partvec[i] == 1.0 && return 0.5
    elseif i == λCnt
        mod(partvec[i], 0.5) == 0.0 && partvec[i] - partvec[i-1] == 1.0 && return 0.5
    else
        if mod(partvec[i], 0.5) == 0.0 && partvec[i] - partvec[i-1] == 1.0
            l = 0.5
        else
            l = 1.0
        end
        if mod(partvec[i], 0.5) == 0.0 && partvec[i+1] - partvec[i] == 1.0
            r = 0.5
        else
            r = 1.0
        end
        return l, r
    end

    return 1.0
end

"""
    Method for INTPROD convexification
"""
function amp_post_special_λ_ub(m::PODNonlinearModel, intprod_indices::Any, dim::Tuple, λ::Dict, d::Dict)

    D = length(intprod_indices)
    tight_regions = [amp_collect_tight_regions(d[i]) for i in intprod_indices]
    checker = [length(d[var])-1 == length(tight_regions[cnt]) for (cnt, var) in enumerate(intprod_indices)]
    prod(checker) ? amp_post_λ_upperbound(m, λ, intprod_indices, (1/2)^D) : amp_post_λ_upperbound(m, λ, intprod_indices, dim, d, tight_regions)

    return
end

"""
    Method for INTLIN convexification
"""
function amp_post_special_λ_lb(m::PODNonlinearModel, intlin_indices::Any, dim::Tuple, λ::Dict, α::Dict, d::Dict)

    D = length(intlin_indices)
    @assert D == 2

    # Functional Scope Check
    for i in intlin_indices
        if m.var_type[i] == :Int && i > m.num_var_orig
            error("Limited support for INTLIN term where INT variable is a lifted variable.")
        end
    end

    tight_regions = Dict(cnt=>amp_collect_tight_regions(d[i]) for (cnt,i) in enumerate(intlin_indices) if m.var_type[i] == :Int)
    amp_post_λ_lowerbound(m, λ, α, intlin_indices, dim, tight_regions)

    return
end

function amp_collect_tight_regions(partvec::Vector)

    PCnt = length(partvec) - 1
    PCnt == 1 && return []

    tight_regions = []
    for i in 1:PCnt
        if i == 1
            (partvec[i]+1.0==partvec[i+1]) && (partvec[i+1]+1.0==partvec[i+2]) && push!(tight_regions,i)
        elseif i == PCnt
            (partvec[i]+1.0==partvec[i+1]) && (partvec[i]-1.0==partvec[i-1]) && push!(tight_regions, i)
        else
            (partvec[i]+1.0==partvec[i+1]) && (partvec[i+1]+1.0==partvec[i+2]) && (partvec[i]-1.0==partvec[i-1]) && push!(tight_regions, i)
        end
    end

    return tight_regions
end

"""
    Limited on two dimensions
    Experimental Code
"""
function amp_post_λ_lowerbound(m::PODNonlinearModel, λ::Dict, α::Dict, indices::Any, dim::Tuple, tregions::Dict)

    @assert length(keys(tregions)) == 1
    @assert length(dim) == 2

    intvars = [i for i in indices if m.var_type[i] == :Int]
    @assert length(intvars) == 1
    intvar = intvars[1]

    for cnt in keys(tregions)
        for i in tregions[cnt]
            sliced_idxs_left = sort(collect_indices(λ[indices][:indices], cnt, [i;], dim))
            sliced_idxs_right = sort(collect_indices(λ[indices][:indices], cnt, [i+1;], dim))
            @assert length(sliced_idxs_left) == length(sliced_idxs_right)
            for k in 1:length(sliced_idxs_left)
                @constraint(m.model_mip, λ[indices][:vars][sliced_idxs_left[k]] <= λ[indices][:vars][sliced_idxs_right[k]] + (1-α[intvar][i]))
                @constraint(m.model_mip, λ[indices][:vars][sliced_idxs_left[k]] >= λ[indices][:vars][sliced_idxs_right[k]] - (1-α[intvar][i]))
            end
        end
    end

    return
end

function amp_post_λ_upperbound(m::PODNonlinearModel, λ::Dict, indices::Any, dim::Tuple, d::Dict, tregions::Vector, reg=[], level=0)

    if level == length(indices)
        isempty(tregions[level]) && return
        sliced_indices = Set(collect_indices(λ[indices][:indices], 1, [reg[1],reg[1]+1;], dim))
        for i in 2:length(reg)
            sliced_indices = intersect(sliced_indices, Set(collect_indices(λ[indices][:indices], i, [reg[i],reg[i]+1], dim)))
        end
        for i in sliced_indices
            setupperbound(λ[indices][:vars][i], (1/2)^level)
        end
        return
    end

    for i in 1:length(tregions[level+1])
        push!(reg, tregions[level+1][i])
        amp_post_λ_upperbound(m, λ, indices, dim, d, tregions, reg, level+1)
        length(reg) < level && error("Something is wrong")
        length(reg) > level && pop!(reg)
    end

    return
end

function amp_post_λ_upperbound(m::PODNonlinearModel, λ::Dict, indices::Any, ub::Float64)

    for i in λ[indices][:vars] setupperbound(i, ub) end

    return
end

function collect_indices(l::Array, locator::Tuple, dim::Tuple)

    k = 0
    indices = Vector{Int}(2^length(dim))
    for i in 1:prod(dim)
        ind = ind2sub(l, i)
        diff = [((ind[i] - locator[i] == 0) || (ind[i] - locator[i] == 1)) for i in 1:length(dim)]
        if prod(diff)
            k +=1
            indices[k] = i
        end
    end

    return indices
end

function collect_indices(l::Array, fixed_dim::Int, fixed_partition::Array, dim::Tuple)

	k = 0
	indices = Vector{Int}(Int(prod(dim)/dim[fixed_dim]*length(fixed_partition)))
	for i in 1:prod(dim)
		ind = ind2sub(l, i)
		if ind[fixed_dim] in fixed_partition
			k += 1
			indices[k] = i
		end
	end

	return indices
end
