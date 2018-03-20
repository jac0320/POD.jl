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
    m.num_slack_vars = var_idx - (m.num_var_orig+m.num_var_linear_mip+m.num_var_nonlinear_mip) - 1
    return
end


function amp_post_affine_slack_constraint(m::PODNonlinearModel, constr_id::Int)

    affine = m.bounding_constr_mip[constr_id]
    if affine[:sense] == :(>=)
        @constraint(m.model_mip, sum(affine[:coefs][j]*Variable(m.model_mip, affine[:vars][j].args[2]) for j in 1:affine[:cnt]) +
                                 Variable(m.model_mip, m.slack_links[constr_id][1]) >= affine[:rhs])
    elseif affine[:sense] == :(<=)
        @constraint(m.model_mip, sum(affine[:coefs][j]*Variable(m.model_mip, affine[:vars][j].args[2]) for j in 1:affine[:cnt]) -
                                 Variable(m.model_mip, m.slack_links[constr_id][1]) <= affine[:rhs])
    elseif affine[:sense] == :(==)
        @constraint(m.model_mip, sum(affine[:coefs][j]*Variable(m.model_mip, affine[:vars][j].args[2]) for j in 1:affine[:cnt])
                                + Variable(m.model_mip, m.slack_links[constr_id][1])
                                - Variable(m.model_mip, m.slack_links[constr_id][2]) == affine[:rhs])
    end

    return
end

function amp_post_slackness_objective(m::PODNonlinearModel)

    slack_start = m.num_var_orig + m.num_var_linear_mip + m.num_var_nonlinear_mip + 1
    slack_end = m.num_var_orig + m.num_slack_vars
    @objective(m.model_mip, Min, sum(Variable(m.model_mip, i) for i in slack_start:slack_end))

    return
end
