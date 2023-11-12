using LRP
using Revise
using Random

let
# Developing an optimal solution 
    # Define instance
    instance = "prins200-10-3";
    # Visualize instance
    display(visualize(instance))
    # Define a random number generator
    seeds = [1010, 1104, 1905, 2104, 2412, 2703, 2704, 2710, 2806, 3009]
    for seed ∈ seeds
        println("\nseed: $seed")
        rng = MersenneTwister(seed);
        # Define inital solution method and build the initial solution
        G  = build(instance)
        sₒ = initialsolution(rng, G, :cluster);
        # Define ALNS parameters
        x = length(sₒ.D)+length(sₒ.C);
        n = max(1000, ceil(x, digits=-(length(digits(x))-1)));
        χ = ALNSparameters(
            j   =   250                     ,
            k   =   125                     ,
            n   =   n ÷ 25                  ,
            m   =   2n                      ,
            Ψᵣ  =   [
                        :randomcustomer!    ,
                        :randomroute!       ,
                        :randomvehicle!     ,
                        :randomdepot!       ,
                        :relatedcustomer!   ,
                        :relatedroute!      ,
                        :relatedvehicle!    ,
                        :relateddepot!      ,
                        :worstcustomer!     ,
                        :worstroute!        ,
                        :worstvehicle!      ,
                        :worstdepot!
                    ]                       ,
            Ψᵢ  =   [
                        :best!              ,
                        :precise!           ,
                        :perturb!           ,
                        :regret2!           ,
                        :regret3!
                    ]                       ,
            Ψₗ  =   [
                        :intraopt!          ,
                        :interopt!          ,
                        :move!              ,
                        :split!             ,
                        :swap!              
                    ]                       ,
            σ₁  =   15                      ,
            σ₂  =   10                      ,
            σ₃  =   3                       ,
            ω̅   =   0.1                     ,
            τ̅   =   0.5                     ,
            ω̲   =   0.01                    ,
            τ̲   =   0.01                    ,
            𝜃   =   0.9975                  ,
            μ̲   =   0.1                     ,
            C̲   =   4                       ,
            μ̅   =   0.4                     ,
            C̅   =   60                      ,
            ρ   =   0.1
        );
        # Run ALNS and fetch best solution
        S = ALNS(rng, χ, sₒ);
        s⃰ = S[end];
    # Fetch objective function values
        println("Objective function value:")
        println("   Initial: $(f(sₒ; penalty=false))")
        println("   Optimal: $(f(s⃰ ; penalty=false))")
    #= Fetch fixed costs
        println("Fixed costs:")
        println("   Initial: $(f(sₒ; operational=false, penalty=false))")
        println("   Optimal: $(f(s⃰ ; operational=false, penalty=false))")
    # Fetch operational costs
        println("Operational costs:")
        println("   Initial: $(f(sₒ; fixed=false, penalty=false))")
        println("   Optimal: $(f(s⃰ ; fixed=false, penalty=false))")
    =#
    # Check if the solutions are feasible
        println("Solution feasibility:")
        println("   Initial: $(isfeasible(sₒ))")
        println("   Optimal: $(isfeasible(s⃰))")
    #= Optimal solution characteristics
        println("Optimal solution characteristics:")
        println("   Number of depots: $(sum([LRP.isopt(d) for d ∈ s⃰.D]))")
        println("   Number of vehicles: $(sum([LRP.isopt(v) for d ∈ s⃰.D for v ∈ d.V]))")
        println("   Number of routes: $(sum([LRP.isopt(r) for d ∈ s⃰.D for v ∈ d.V for r ∈ v.R]))")
    =#
    # Visualizations
        # Visualize initial solution
        display(visualize(sₒ))
        # Visualize best solution
        display(visualize(s⃰))
        # Animate ALNS solution search process from inital to best solution
        #display(animate(S))
        # Show convergence plot
        display(pltcnv(S; penalty=true))
        display(pltcnv(S; penalty=false))
    end
    return
end
