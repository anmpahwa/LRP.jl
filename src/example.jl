using LRP
using Revise
using Random

let
# Developing an optimal TSP route 
    # Define instance
    instance = "cmt10"
    # Define a random number generator
    rng = MersenneTwister(1234)
    # Build instance as graph
    G = build(instance)
    D, C, A = G
    # Define ALNS parameters
    χₒ  = ObjectiveFunctionParameters(
        d = 0.                          ,
        v = 100000.                     ,
        r = 0.                          ,
        c = 0.                          ,
    )
    χ   = ALNSParameters(
        k̲   =   6                       ,
        k̅   =   1500                    ,
        k̲ₛ  =   240                     ,
        k̅ₛ  =   750                     ,   
        Ψᵣ  =   [
                    :randomnode!    , 
                    :relatedpair!   ,
                    :relatednode!   , 
                    :worstnode!     ,  
                    :randomroute!   ,
                    :relatedroute!  ,
                    :worstroute!    ,
                    :randomvehicle! 
                ]                       , 
        Ψᵢ  =   [
                    :best!          ,
                    :greedy!        ,
                    :regret₂insert! ,
                    :regret₃insert!
                ]                       ,
        Ψₛ  =   [
                    :move!          ,
                    :split!         ,
                    :swap!
                ]                       ,
        σ₁  =   33                      ,
        σ₂  =   9                       ,
        σ₃  =   13                      ,
        ω   =   0.05                    ,
        τ   =   0.5                     ,
        𝜃   =   0.99975                 ,
        C̲   =   30                      ,
        C̅   =   60                      ,
        μ̲   =   0.1                     ,
        μ̅   =   0.4                     ,
        ρ   =   0.1                     ,
        χₒ  =   χₒ  
    )
    # Define inital solution method and build the initial solution
    method = :regret₂init
    sₒ = initialsolution(rng, G, χₒ, method)
    # Run ALNS and fetch best solution
    S = ALNS(rng, sₒ, χ)
    s⃰ = S[end]
            
# Fetch objective function values
    println("Initial: $(f(sₒ, χₒ))")
    println("Optimal: $(f(s⃰,  χₒ))")

# Visualizations
    # Visualize initial solution
    display(visualize(sₒ)) 
    # Visualize best solution
    display(visualize(s⃰))
    # Animate ALNS solution search process from inital to best solution
    display(animate(S))
    # Show convergence plots
    display(convergence(S, χₒ))
    
    return
end