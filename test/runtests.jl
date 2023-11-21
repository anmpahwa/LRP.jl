using LRP
using Revise
using Test
using Random

let
    # ALNS parameters
    χ = ALNSparameters(
        j   =   50                      ,
        k   =   20                      ,
        n   =   20                      ,
        m   =   500                     ,
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
                    :split!             ,
                    :intraopt!          ,
                    :move!              ,
                    :swap!              ,
                    :interopt!          
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
        θ   =   0.9985                  ,
        ρ   =   0.1
    );

    # Vehicle Routing Problem with time-windows
    @testset "VRPTW" begin
        instances = ["r101", "c101"]
        methods = [:cluster, :random]
        for instance ∈ instances
            visualize(instance)
            for method ∈ methods
                println("\n Solving $instance using ALNS metaheuristic initialized with $method method")
                G = build(instance)
                sₒ = initialsolution(G, method)     
                s⃰ = ALNS(χ, sₒ)
                visualize(s⃰)
                @test isfeasible(s⃰)
                @test f(s⃰) ≤ f(sₒ)
            end
        end
    end


    # Location Routing Problem
    @testset "LRP" begin
        instances = ["prins20-5-1", "prins50-5-1b"]
        methods = [:cluster, :random]
        for instance ∈ instances
            visualize(instance)
            for method ∈ methods
                println("\n Solving $instance using ALNS metaheuristic initialized with $method method")
                G = build(instance)
                sₒ = initialsolution(G, method)     
                s⃰ = ALNS(χ, sₒ)
                visualize(s⃰)
                @test isfeasible(s⃰)
                @test f(s⃰) ≤ f(sₒ)
            end
        end
    end
    
    return
end
        
