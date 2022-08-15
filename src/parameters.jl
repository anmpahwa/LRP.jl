@doc """
    ALNSParameters

Optimization parameters for Adaptive Large Neighborhood Search (ALNS).

- k̲     :   Number of ALNS iterations triggering operator probability update (segment size)
- l̲     :   Number of ALNS iterations triggering local search
- l̅     :   Number of local search iterations
- k̅     :   Number of ALNS iterations
- Ψᵣ    :   Vector of removal operators
- Ψᵢ    :   Vector of insertion operators
- Ψₗ    :   Vector of local search operators
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
"""
Base.@kwdef struct ALNSParameters
    k̲::Int64
    l̲::Int64
    l̅::Int64
    k̅::Int64
    Ψᵣ::Vector{Symbol}
    Ψᵢ::Vector{Symbol}
    Ψₗ::Vector{Symbol}
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
end