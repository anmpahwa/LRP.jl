using LRP
using Revise
using Random
using CPUTime
using DataFrames

let
    # Set A
    A = ["prins20-5-1", "gaskell36-5", "prins50-5-1b", "daskin88-8"];
    # Set B
    B = ["prins100-5-2", "prins100-10-2b", "christofides100-10"]
    # Set C
    C = ["min134-8", "daskin150-10", "prins200-10-3"] 
    # Define instances
    instances = [A..., B..., C...]
    # Define random number generators
    seeds = [1010, 1106, 1509, 1604, 1905, 2104, 2412, 2703, 2710, 2807]
    # Dataframes to store solution quality and run time
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
            sₒ = initialize(rng, instance; method=:global);
            # Visualize initial solution
            display(visualize(sₒ))
            # Define ALNS parameters
            x = max(100, lastindex(sₒ.C))
            χ = ALNSparameters(
                j   =   50                      ,
                k   =   5                       ,
                n   =   x                       ,
                m   =   100x                    ,
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
                            :intramove!         ,
                            :intraswap!         ,
                            :intraopt!          ,
                            :intermove!         ,
                            :interswap!         ,
                            :interopt!          ,
                            :swapdepot!
                        ]                       ,
                σ₁  =   15                      ,
                σ₂  =   10                      ,
                σ₃  =   3                       ,
                μ̲   =   0.1                     ,
                c̲   =   4                       ,
                μ̅   =   0.4                     ,
                c̅   =   60                      ,
                ω̅   =   0.05                    ,
                τ̅   =   0.5                     ,
                ω̲   =   0.01                    ,
                τ̲   =   0.01                    ,
                θ   =   0.9985                  ,
                ρ   =   0.1
            );
            # Run ALNS and fetch best solution
            t = @CPUelapsed s⃰ = ALNS(rng, χ, sₒ);
            # Visualize best solution
            display(visualize(s⃰))
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
            println("   Number of depots: $(sum([!iszero(d.n) for d ∈ s⃰.D]))")
            println("   Number of vehicles: $(sum([!iszero(v.n) for d ∈ s⃰.D for v ∈ d.V]))")
            println("   Number of routes: $(sum([!iszero(r.n) for d ∈ s⃰.D for v ∈ d.V for r ∈ v.R]))")
            # Store Results
            df₁[i,j+1] = f(s⃰)
            df₂[i,j+1] = t
            println(df₁)
            println(df₂)
        end
    end
    return
end