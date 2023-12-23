using LRP
using Revise
using Test
using Random

let
    # ALNS parameters
    χ = ALNSparameters(
        j   =   50                      ,
        k   =   5                       ,
        n   =   10                      ,
        m   =   1000                    ,
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
     # Vehicle Routing Problem
     @testset "VRP" begin
        instances = ["m-n101-k10", "tai150a", "cmt10"]
        for instance ∈ instances
            visualize(instance)
            println(instance)
            s₁ = initialize(instance)
            s₂ = ALNS(χ, s₁)
            visualize(s₂)
            @test isfeasible(s₂)
            @test f(s₂) ≤ f(s₁)
        end
    end
    # Vehicle Routing Problem with time-windows
    @testset "VRPTW" begin
        instances = ["r101", "c101", "rc101"]
        for instance ∈ instances
            visualize(instance)
            println(instance)
            s₁ = initialize(instance)
            s₂ = ALNS(χ, s₁)
            visualize(s₂)
            @test isfeasible(s₂)
            @test f(s₂) ≤ f(s₁)
        end
    end
    # Location Routing Problem
    @testset "LRP" begin
        instances = ["prins20-5-1", "gaskell36-5", "prins50-5-1b"]
        for instance ∈ instances
            visualize(instance)
            println("\n $instance")
            s₁ = initialize(instance)
            s₂ = ALNS(χ, s₁)
            visualize(s₂)
            @test isfeasible(s₂)
            @test f(s₂) ≤ f(s₁)
        end
    end
    return
end
        
