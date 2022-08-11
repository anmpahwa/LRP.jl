@doc """
    ALNSParameters

- k̲   :   ALNS segment size
- k̅   :   ALNS iterations
- k̲ₛ  :   Local Search segment size
- k̅ₛ  :   Local Search iterations 
- Ψᵣ  :   Vector of removal operators
- Ψᵢ  :   Vector of insertion operators
- Ψₛ  :   Vector of local search operators
- σ₁  :   Score for a new best solution
- σ₂  :   Score for a new better solution
- σ₃  :   Score for a new worse but accepted solution
- ω   :   Start tempertature control threshold 
- τ   :   Start tempertature control probability
- 𝜃   :   Cooling rate
- C̲   :   Minimum customer nodes removal
- C̅   :   Maximum customer nodes removal
- μ̲   :   Minimum removal fraction
- μ̅   :   Maximum removal fraction
- ρ   :   reaction factor

"""
Base.@kwdef struct ALNSParameters
    k̲::Int64                                    # ALNS segment size
    k̅::Int64                                    # ALNS iterations
    k̲ₛ::Int64                                   # Local Search segment size
    k̅ₛ::Int64                                   # Local Search iterations
    Ψᵣ::Vector{Symbol}                          # Vector of removal operators
    Ψᵢ::Vector{Symbol}                          # Vector of insertion operators
    Ψₛ::Vector{Symbol}                          # Vector of local search operators
    σ₁::Float64                                 # Score for a new best solution
    σ₂::Float64                                 # Score for a new better solution
    σ₃::Float64                                 # Score for a new worse solution
    ω::Float64                                  # Start temperature control threshold
    τ::Float64                                  # Start temperature control probability
    𝜃::Float64                                  # Cooling rate
    C̲::Int64                                    # Minimum customer nodes removal
    C̅::Int64                                    # Maximum customer nodes removal
    μ̲::Float64                                  # Minimum removal fraction
    μ̅::Float64                                  # Maximum removal fraction
    ρ::Float64                                  # Reaction factor
end