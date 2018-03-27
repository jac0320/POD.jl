function initialize_slackness_link(m::PODNonlinearModel)

    m.slack_links = Dict()
    var_idx = m.num_var_orig+m.num_var_linear_mip+m.num_var_nonlinear_mip + 1
    for i in 1:m.num_constr_orig
        if m.constr_type_orig[i] in [:(>=), :(<=)]
            m.slack_links[i] = (var_idx)
            var_idx += 1
        else
            m.slack_links[i] = (var_idx, var_idx+1)
            var_idx += 2
        end
    end

    # Stores the maximum slackness variable index on 0
    m.slack_links[0] = (var_idx)

    # Measure the total number of slack variables
    m.num_slack_vars = var_idx - (m.num_var_orig+m.num_var_linear_mip+m.num_var_nonlinear_mip)
    return
end


function amp_post_affine_slack_constraint(m::PODNonlinearModel, constr_id::Int)

    affine = m.bounding_constr_mip[constr_id]
    if affine[:sense] == :(>=)
        @constraint(m.model_mip, sum(affine[:coefs][j]*Variable(m.model_mip, affine[:vars][j].args[2]) for j in 1:affine[:cnt]) +
                                 Variable(m.model_mip, m.slack_links[constr_id][1]) >= affine[:rhs])
        @constraint(m.model_mip, Variable(m.model_mip, m.slack_links[constr_id][1]) <= Variable(m.model_mip, m.slack_links[0]))
    elseif affine[:sense] == :(<=)
        @constraint(m.model_mip, sum(affine[:coefs][j]*Variable(m.model_mip, affine[:vars][j].args[2]) for j in 1:affine[:cnt]) -
                                 Variable(m.model_mip, m.slack_links[constr_id][1]) <= affine[:rhs])
        @constraint(m.model_mip, Variable(m.model_mip, m.slack_links[constr_id][1]) <= Variable(m.model_mip, m.slack_links[0]))
    elseif affine[:sense] == :(==)
        @constraint(m.model_mip, sum(affine[:coefs][j]*Variable(m.model_mip, affine[:vars][j].args[2]) for j in 1:affine[:cnt])
                                + Variable(m.model_mip, m.slack_links[constr_id][1])
                                - Variable(m.model_mip, m.slack_links[constr_id][2]) == affine[:rhs])
        @constraint(m.model_mip, Variable(m.model_mip, m.slack_links[constr_id][1]) <= Variable(m.model_mip, m.slack_links[0]))
        @constraint(m.model_mip, Variable(m.model_mip, m.slack_links[constr_id][2]) <= Variable(m.model_mip, m.slack_links[0]))
    end

    return
end

function amp_post_slackness_objective(m::PODNonlinearModel)

    @objective(m.model_mip, Min, Variable(m.model_mip, m.slack_links[0]))

    return
end

function amp_post_zero_objective(m::PODNonlinearModel)

    @objective(m.model_mip, Min, 0.0)

    return
end
