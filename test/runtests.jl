using LRP
using Revise
using Test
using Random

let
    # Vehicle Routing Problem with time-windows
    @testset "VRPTW" begin
        χ   = ALNSparameters(
            n   =   4                       ,
            k   =   250                     ,
            m   =   200                     ,
            j   =   125                     ,
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
                        :movecustomer!      ,
                        :movedepot!         ,
                        :swapcustomers!     ,
                        :swapdepots!
                    ]                       ,
            σ₁  =   15                      ,
            σ₂  =   10                      ,
            σ₃  =   3                       ,
            ω   =   0.05                    ,
            τ   =   0.5                     ,
            𝜃   =   0.9975                  ,
            C̲   =   4                       ,
            C̅   =   60                      ,
            μ̲   =   0.1                     ,
            μ̅   =   0.4                     ,
            ρ   =   0.1
        )
        instances = ["r101", "c101"]
        methods = [:cluster, :random]
        for k ∈ 1:2
            instance = instances[k]
            method = methods[k]
            println("\nSolving $instance")
            visualize(instance)
            G  = build(instance)
            sₒ = initialsolution(G, method)         
            S  = ALNS(χ, sₒ)
            s⃰  = S[end]
            visualize(s⃰)
            pltcnv(S)
            @test isfeasible(s⃰)
            @test f(s⃰) ≤ f(sₒ)
        end
    end

    # Location Routing Problem
    @testset "LRP" begin
        χ   = ALNSparameters(
            n   =   4                       ,
            k   =   250                     ,
            m   =   200                     ,
            j   =   125                     ,
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
                        :movecustomer!      ,
                        :movedepot!         ,
                        :swapcustomers!     ,
                        :swapdepots!
                    ]                       ,
            σ₁  =   15                      ,
            σ₂  =   10                      ,
            σ₃  =   3                       ,
            ω   =   0.05                    ,
            τ   =   0.5                     ,
            𝜃   =   0.9975                  ,
            C̲   =   4                       ,
            C̅   =   60                      ,
            μ̲   =   0.1                     ,
            μ̅   =   0.4                     ,
            ρ   =   0.1
        )
        instances = ["prins20-5-1", "prins50-5-1b"]
        methods = [:cluster, :random]
        for k ∈ 1:2
            instance = instances[k]
            method = methods[k]
            println("\nSolving $instance")
            visualize(instance)
            G  = build(instance)
            sₒ = initialsolution(G, method)     
            S  = ALNS(χ, sₒ)
            s⃰  = S[end]
            visualize(s⃰)
            pltcnv(S)
            @test isfeasible(s⃰)
            @test f(s⃰) ≤ f(sₒ)
        end
    end
    return
end
        
