using LRP
using Revise
using Random

let 
# Developing an optimal Traveling Salesman Problem solution
    # Define instance
    instance = "att48"
    # Define a random number generator
    rng = MersenneTwister(1234)
    # Build instance as graph
    G = build(instance)
    D, C, A = G
    # Define ALNS parameters
    χ   = ALNSParameters(
        k̲   =   6                       ,
        l̲   =   240                     ,
        l̅   =   750                     ,
        k̅   =   1500                    ,
        Ψᵣ  =   [
                    :randomnode!    , 
                    :relatedpair!   ,
                    :relatednode!   , 
                    :worstnode!     ,  
                ]                       , 
        Ψᵢ  =   [
                    :best!          ,
                    :greedy!        ,
                    :regret₂insert! ,
                    :regret₃insert!
                ]                       ,
        Ψₗ  =   [
                    :move!          ,
                    :intraopt!      ,
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
        ρ   =   0.1
    )
    # Define inital solution method and build the initial solution
    method = :cw
    sₒ = initialsolution(rng, G, method)
    # Run ALNS and fetch best solution
    S = ALNS(rng, sₒ, χ)
    s⃰ = S[end]
            
# Fetch objective function values
    println("Initial: $(f(sₒ))")
    println("Optimal: $(f(s⃰))")

# Visualizations
    # Visualize initial solution
    display(visualize(sₒ)) 
    # Visualize best solution
    display(visualize(s⃰))
    # Animate ALNS solution search process from inital to best solution
    display(animate(S))
    # Show convergence plots
    display(convergence(S))
    return
end

let
# Developing an optimal Single-Depot Vehicle Routing Problem solution 
    # Define instance
    instance = "cmt10"
    # Define a random number generator
    rng = MersenneTwister(1234)
    # Build instance as graph
    G = build(instance)
    D, C, A = G
    # Define ALNS parameters
    χ   = ALNSParameters(
        k̲   =   6                       ,
        l̲   =   240                     ,
        l̅   =   750                     ,
        k̅   =   1500                    ,
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
        Ψₗ  =   [
                    :move!          ,
                    :intraopt!      ,
                    :interopt!      ,
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
        ρ   =   0.1
    )
    # Define inital solution method and build the initial solution
    method = :random
    sₒ = initialsolution(rng, G, method)
    # Run ALNS and fetch best solution
    S = ALNS(rng, sₒ, χ)
    s⃰ = S[end]
            
# Fetch objective function values
    println("Initial: $(f(sₒ))")
    println("Optimal: $(f(s⃰))")

# Visualizations
    # Visualize initial solution
    display(visualize(sₒ)) 
    # Visualize best solution
    display(visualize(s⃰))
    # Animate ALNS solution search process from inital to best solution
    display(animate(S))
    # Show convergence plots
    display(convergence(S))
    return
end