"""
    ALNSparameters

Optimization parameters for Adaptive Large Neighborhood Search (ALNS).

- j     :   Number of segments in the ALNS
- k     :   Number of segments to reset ALNS
- n     :   Number of iterations in an ALNS segment
- m     :   Number of iterations in a local search
- Ψᵣ    :   Vector of removal operators
- Ψᵢ    :   Vector of insertion operators
- Ψₗ    :   Vector of local search operators
- σ₁    :   Score for a new best solution
- σ₂    :   Score for a new better solution
- σ₃    :   Score for a new worse but accepted solution
- μ̲     :   Minimum removal fraction
- C̲     :   Minimum customer nodes removed
- μ̅     :   Maximum removal fraction
- C̅     :   Maximum customer nodes removed
- ω̅     :   Initial temperature deviation parameter
- τ̅     :   Initial temperatureprobability parameter
- ω̲     :   Final temperature deviation parameter
- τ̲     :   Final temperature probability parameter
- φ     :   Local search trigger
- θ     :   Cooling rate
- ρ     :   Reaction factor
"""
Base.@kwdef struct ALNSparameters
    j::Int64
    k::Int64
    n::Int64
    m::Int64
    Ψᵣ::Vector{Symbol}
    Ψᵢ::Vector{Symbol}
    Ψₗ::Vector{Symbol}
    σ₁::Float64
    σ₂::Float64
    σ₃::Float64
    μ̲::Float64
    C̲::Int64
    μ̅::Float64
    C̅::Int64
    ω̅::Float64
    τ̅::Float64
    ω̲::Float64
    τ̲::Float64
    φ::Float64
    θ::Float64
    ρ::Float64
end