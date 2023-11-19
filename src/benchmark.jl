using LRP
using CSV
using Revise
using Random
using CPUTime
using DataFrames

let
# Developing an optimal solution 
    # Define instances
    instances = ["c101", "min134-8", "daskin150-10"];
    # Define a random number generators
    seeds = [1010, 1104, 1905, 2104, 2412, 2703, 2704, 2710, 2806, 3009]
    # Dataframes
    df₁ = DataFrame([instances, [zeros(length(instances)) for _ ∈ seeds]...], [iszero(i) ? "instance" : "$(seeds[i])" for i ∈ 0:length(seeds)])
    df₂ = DataFrame([instances, [zeros(length(instances)) for _ ∈ seeds]...], [iszero(i) ? "instance" : "$(seeds[i])" for i ∈ 0:length(seeds)])
    for i ∈ eachindex(instances)
        instance = instances[i]
        # Visualize instance
        display(visualize(instance))
        for j ∈ eachindex(seeds)
            seed = seeds[j]
            println("\n instance: $instance | seed: $seed")
            rng = MersenneTwister(seed);
            # Define inital solution method and build the initial solution
            D, C, A = build(instance)
            G = (D, C, A)
            sₒ= initialsolution(rng, G, :cluster);
            # Define ALNS parameters
            x = max(100, lastindex(C))
            χ = ALNSparameters(
                j   =   50                      ,
                k   =   15                      ,
                n   =   x                       ,
                m   =   25x                     ,
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
                μ̲   =   0.1                     ,
                C̲   =   4                       ,
                μ̅   =   0.4                     ,
                C̅   =   60                      ,
                ω̅   =   0.2                     ,
                τ̅   =   0.5                     ,
                ω̲   =   0.01                    ,
                τ̲   =   0.01                    ,
                φ   =   0.25                    ,
                θ   =   0.998                   ,
                ρ   =   0.1
            );
            # Run ALNS and fetch best solution
            t = @CPUelapsed S = ALNS(rng, χ, sₒ);
            s⃰ = S[argmin(f.(S))];
        # Fetch objective function values
            println("Objective function value:")
            println("   Initial: $(f(sₒ; penalty=false))")
            println("   Optimal: $(f(s⃰ ; penalty=false))")
        # Fetch fixed costs
            println("Fixed costs:")
            println("   Initial: $(f(sₒ; operational=false, penalty=false))")
            println("   Optimal: $(f(s⃰ ; operational=false, penalty=false))")
        # Fetch operational costs
            println("Operational costs:")
            println("   Initial: $(f(sₒ; fixed=false, penalty=false))")
            println("   Optimal: $(f(s⃰ ; fixed=false, penalty=false))")
        # Check if the solutions are feasible
            println("Solution feasibility:")
            println("   Initial: $(isfeasible(sₒ))")
            println("   Optimal: $(isfeasible(s⃰))")
        # Optimal solution characteristics
            println("Optimal solution characteristics:")
            println("   Number of depots: $(sum([LRP.isopt(d) for d ∈ s⃰.D]))")
            println("   Number of vehicles: $(sum([LRP.isopt(v) for d ∈ s⃰.D for v ∈ d.V]))")
            println("   Number of routes: $(sum([LRP.isopt(r) for d ∈ s⃰.D for v ∈ d.V for r ∈ v.R]))")
        # Visualizations
            # Visualize initial solution
            display(visualize(sₒ))
            # Visualize best solution
            display(visualize(s⃰))
            # Animate ALNS solution search process from inital to best solution
            display(animate(S))
            # Show convergence plot
            display(pltcnv(S))
        # Store Results
            df₁[i,j+1] = f(s⃰)
            df₂[i,j+1] = t
        end
    end
    display(df₁)
    display(df₂)
    return
end
