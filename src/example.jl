using LRP
using Revise
using Random

let
# Developing an optimal solution 
    # Define instance
    instance = "x-n303-k21"
    # Define a random number generator
    rng = MersenneTwister(1234)
    # Build instance as graph
    G = build(instance)
    D, C, A = G
    # Visualize instance
    display(visualize(Solution(G...)))
    # Define inital solution method and build the initial solution
    method = :cw
    sₒ = initialsolution(rng, G, method)
    # Define ALNS parameters
    x = length(C) + length(D)
    n = ceil(x, digits=-(length(digits(x))-1))
    χ   = ALNSParameters(
        k̲   =   n ÷ 25                  ,
        l̲   =   2n                      ,
        l̅   =   5n                      ,
        k̅   =   10n                     ,
        Ψᵣ  =   [
                    :randomnode!    , 
                    :randomroute!   ,
                    :randomvehicle! ,
                    :randomdepot!   ,
                    :relatednode!   , 
                    :relatedroute!  ,  
                    :relatedvehicle!,
                    :relateddepot!  ,
                    :worstnode!     ,
                    :worstroute!    ,
                    :worstvehicle!  ,
                    :worstdepot!
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
    # Run ALNS and fetch best solution
    S = ALNS(rng, χ, sₒ)
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
