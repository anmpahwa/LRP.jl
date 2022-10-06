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

const M = typemax(Int64)

include("sample.jl")
include("datastructure.jl")
include("instance.jl")
include("objfunction.jl")
include("feasible.jl")
include("initialize.jl")
include("operations.jl")
include("relatedness.jl")
include("remove.jl")
include("insert.jl")
include("localsearch.jl")
include("parameters.jl")
include("ALNS.jl")
include("visualize.jl")

export  f, isfeasible, isopt,
        ALNSParameters, initialsolution, ALNS, 
        visualize, vectorize, animate, pltcnv

end