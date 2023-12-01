module LRP

using Clustering
using CSV
using DataFrames
using Distributions
using ElasticArrays
using OffsetArrays
using Plots
using ProgressMeter
using Random
using StatsBase

ElasticArrays.ElasticMatrix(A::OffsetMatrix) = OffsetMatrix(ElasticArray(A), A.offsets)
Base.append!(A::OffsetMatrix, items) = (append!(A.parent, items); A)

global φᵀ = false::Bool

include("sample.jl")
include("datastructure.jl")
include("functions.jl")
include("initialize.jl")
include("operations.jl")
include("remove.jl")
include("insert.jl")
include("localsearch.jl")
include("parameters.jl")
include("ALNS.jl")
include("visualize.jl")

export  initialize, vectorize, f, isfeasible,
        ALNSparameters, ALNS, visualize, animate, pltcnv

end

# TODO: Improve efficiency of local search methods to avoid enumeration of customers within a route. (move!(rng, k̅, s); opt!(rng, k̅, s))
# TODO: Test randomizing depot insertion position in a route. (swapdepot!(rng, k̅, s))
# TODO: Test updating constraint violation penalty evaluation. f(s)
# TODO: Test minimalizing addroute through improved constraint violation penalty evaluation. addroute(r, s)