# Include global environment so RCall is accessible
push!(LOAD_PATH, "@v#.#")
include("test_core.jl")
