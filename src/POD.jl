__precompile__()

module POD

PODDEBUG = true

using JuMP
using MathProgBase
using Compat

using Requires: @require

@require Gurobi begin
    import Gurobi
end

include("const.jl")

# Engine for High-level Algorithmic Control and User-interface
include("solver.jl")
include("mpb2moi.jl") # Transition file

# Engine for expression handling
include("nlexpr.jl")
include("operators.jl")

# Main Algorithmic Process
include("algorithm.jl")
include("vars.jl")
include("presolve.jl")
include("amp.jl")
include("embedding.jl")
include("heuristics.jl")
include("feasibility.jl")

# Convexification method
include("multi.jl")
include("tmc.jl")

# Model Manipulation and utilities
include("bounds.jl")
include("utility.jl")
include("integer.jl")
include("trigonometry.jl")

# Othes
include("log.jl")

end # module
