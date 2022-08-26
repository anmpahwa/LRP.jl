using LRP
using Revise
using Test
using Random

let
    χ   = ALNSParameters(
        k̲   =   50                      ,
        l̲   =   50                      ,
        l̅   =   250                     ,
        k̅   =   500                     ,
        Ψᵣ  =   [
                    :randomnode!    , 
                    :relatednode!   , 
                    :worstnode!     ,
                    :randomroute!   ,
                    :relatedroute!  ,  
                    :worstroute!    ,
                    :randomvehicle! ,
                    :relatedvehicle!,
                    :worstvehicle!  ,
                    :randomdepot!   ,
                    :relateddepot!  ,
                    :worstdepot!
                ]                       , 
        Ψᵢ  =   [
                    :best!          ,
                    :greedy!        ,
                    :regret2!       ,
                    :regret3!
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
    # Traveling Salesman Problem instances
    @testset "TSP" begin
        instances = ["att48", "eil101", "ch150", "d198", "a280"]
        for instance ∈ instances
            println("\nSolving $instance")
            sₒ = initialsolution(instance, :random)     
            S  = ALNS(χ, sₒ)
            s⃰  = S[end]
            @test isfeasible(sₒ)
            @test isfeasible(s⃰)
            @test f(s⃰) ≤ f(sₒ)
        end
    end
    # Vehicle Routing Problem
    @testset "VRP" begin
        instances = ["m-n101-k10", "tai150a", "cmt10", "x-n251-k28", "x-n303-k21"]
        for instance ∈ instances
            println("\nSolving $instance")
            sₒ = initialsolution(instance, :random)     
            S  = ALNS(χ, sₒ)
            s⃰  = S[end]
            @test isfeasible(s⃰)
            @test f(s⃰) ≤ f(sₒ)
        end
    end
    return
end