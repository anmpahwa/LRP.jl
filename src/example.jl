using LRP
using Revise
using Random

let
# Developing an optimal TSP route 
    # Define instance
    instance = "att48"
    # Define a random number generator
    rng = MersenneTwister(1234)
    # Build instance as graph
    G = build(instance)
    D, C, A = G
    # Define ALNS parameters
    χₒ  = ObjectiveFunctionParameters(
        d = 0.                          ,
        v = 0.                          ,
        r = 0.                          ,
        c = 0.                          ,
    )
    χ   = ALNSParameters(
        k̲   =   2                       ,
        k̅   =   500                     ,
        k̲ₛ  =   80                      ,
        k̅ₛ  =   250                     ,   
        Ψᵣ  =   [
                    :node_remove! , 
                    :worst_remove!  , 
                    :shaw_remove!
                ]                       , 
        Ψᵢ  =   [
                    :best_insert!   ,
                    :greedy_insert! ,
                    :regret₂insert! ,
                    :regret₃insert!
                ]                       ,
        Ψₛ  =   [
                    :move!              ,
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
    method = :cw_init
    sₒ = initialsolution(rng, G, χₒ, method)
    # Run ALNS and fetch best solution
    S = ALNS(rng, sₒ, χ)
    s⃰ = S[end]

# Visualizations
    # Visualize initial solution
    display(visualize(sₒ)) 
    # Visualize best solution   
    display(visualize(s⃰))
    # Animate ALNS solution search process from inital to best solution
    display(animate(S))
    # Show convergence plots
    display(convergence(S, χₒ))

# Fetch objective function values
    println("Initial: $(f(sₒ, χₒ))")
    println("Optimal: $(f(s⃰,  χₒ))")
    
    return
end