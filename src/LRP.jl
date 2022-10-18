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

export  build, initialsolution, f, isfeasible,
        ALNSParameters, ALNS, 
        vectorize, visualize, animate, pltcnv

end