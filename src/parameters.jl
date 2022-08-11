@doc """
    ObjectiveFunctionParameters

Parameters for objective function evaluation.

- d     :   Depot constraint parameter
- v     :   Vehicle constraint parameter
- r     :   Route constraint parameter
- c     :   Customer constraint parameter
"""
Base.@kwdef mutable struct ObjectiveFunctionParameters
    d::Float64
    v::Float64
    r::Float64
    c::Float64
end

@doc """
    ALNSParameters

Optimization parameters for Adaptive Large Neighborhood Search (ALNS).

- k̲     :   ALNS segment size
- k̅     :   ALNS iterations
- k̲ₛ    :   Local Search segment size
- k̅ₛ    :   Local Search iterations 
- Ψᵣ    :   Vector of removal operators
- Ψᵢ    :   Vector of insertion operators
- Ψₛ    :   Vector of local search operators
- σ₁    :   Score for a new best solution
- σ₂    :   Score for a new better solution
- σ₃    :   Score for a new worse but accepted solution
- ω     :   Start tempertature control threshold 
- τ     :   Start tempertature control probability
- 𝜃     :   Cooling rate
- C̲     :   Minimum customer nodes removal
- C̅     :   Maximum customer nodes removal
- μ̲     :   Minimum removal fraction
- μ̅     :   Maximum removal fraction
- ρ     :   Reaction factor
- χₒ    :   Objective function parameters
"""
Base.@kwdef struct ALNSParameters
    k̲::Int64
    k̅::Int64
    k̲ₛ::Int64
    k̅ₛ::Int64
    Ψᵣ::Vector{Symbol}
    Ψᵢ::Vector{Symbol}
    Ψₛ::Vector{Symbol}
    σ₁::Float64
    σ₂::Float64
    σ₃::Float64
    ω::Float64
    τ::Float64
    𝜃::Float64
    C̲::Int64
    C̅::Int64
    μ̲::Float64
    μ̅::Float64
    ρ::Float64
    χₒ::ObjectiveFunctionParameters
end
