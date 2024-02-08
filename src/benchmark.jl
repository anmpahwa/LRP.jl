using LRP
using Revise
using Random
using CPUTime
using DataFrames

let
    # Set A
    A = ["c101", "c201", "r101", "r201", "rc101", "rc201"]
    # Set B
    B = ["prins20-5-1", "gaskell36-5", "prins50-5-1b", "daskin88-8"];
    # Set C
    C = ["prins100-5-2", "prins100-10-2b", "christofides100-10"]
    # Set D
    D = ["min134-8", "daskin150-10", "prins200-10-3"] 
    # Define instances
    instances = [A..., B..., C..., D...]
    # Define random number generators
    seeds = [1010, 1104, 1509, 1604, 1905, 2104, 2412, 2703, 2710, 2807]
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
            s₁ = initialize(rng, instance; method=:global);
            # Visualize initial solution
            display(visualize(s₁))
            # Define ALNS parameters
            x = max(100, lastindex(s₁.C))
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
            t = @CPUelapsed s₂ = ALNS(rng, χ, s₁);
            # Visualize best solution
            display(visualize(s₂))
            # Fetch objective function values
            println("Objective function value:")
            println("   Initial: $(round(s₁.πᶠ + s₁.πᵒ, digits=3))")
            println("   Optimal: $(round(s₂.πᶠ + s₂.πᵒ, digits=3))")
            # Check if the solutions are feasible
            println("Solution feasibility:")
            println("   Initial: $(isfeasible(s₁)) | $(round(s₁.πᵖ, digits=3))")
            println("   Optimal: $(isfeasible(s₂)) | $(round(s₂.πᵖ, digits=3))")
            # Optimal solution characteristics
            println("Optimal solution characteristics:")
            nᵈ, nᵛ, nʳ = 0, 0, 0
            for d ∈ s₂.D nᵈ += LRP.isopt(d) end
            for d ∈ s₂.D for v ∈ d.V nᵛ += LRP.isopt(v) end end
            for d ∈ s₂.D for v ∈ d.V for r ∈ v.R nʳ += LRP.isopt(r) end end end
            println("   Number of depots: $nᵈ")
            println("   Number of vehicles: $nᵛ")
            println("   Number of routes: $nʳ")
            # Store Results
            df₁[i,j+1] = f(s₂)
            df₂[i,j+1] = t
            println(df₁)
            println(df₂)
        end
    end
    return
end
