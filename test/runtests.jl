using LRP
using Revise
using Test
using Random

# Traveling Salesman Problem instances
@testset "TSP" begin
    K = 5
    instances = ["att48", "eil101", "ch150", "d198", "a280"]
    methods   = [:cw, :nn, :random, :regret₂init, :regret₃init]
    χₒ  = ObjectiveFunctionParameters(
        d = 0.                          ,
        v = 0.                          ,
        r = 0.                          ,
        c = 0.                          ,
    )
    χ   = ALNSParameters(
        k̲   =   2                       ,
        k̅   =   500                     ,
        k̲ₛ  =   80                      ,
        k̅ₛ  =   250                     ,   
        Ψᵣ  =   [
                    :randomnode!    , 
                    :relatedpair!   ,
                    :relatednode!   , 
                    :worstnode!     
                ]                       , 
        Ψᵢ  =   [
                    :best!          ,
                    :greedy!        ,
                    :regret₂insert! ,
                    :regret₃insert!
                ]                       ,
        Ψₛ  =   [
                    :move!              ,
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
        ρ   =   0.1                     ,
        χₒ  =   χₒ  
    )
    for k ∈ 1:K
        instance = instances[k]
        method = methods[k]
        println("\n Solving $instance")
        G = build(instance)
        sₒ= initialsolution(G, χₒ, method)     
        @test isfeasible(sₒ)
        S = ALNS(sₒ, χ)
        s⃰ = S[end]
        @test isfeasible(s⃰)
        @test f(s⃰, χₒ) ≤ f(sₒ, χₒ)
    end
    return
end

# Single Depot Vehicle Routing Problem
@testset "SDVRP" begin
    K = 3
    instances = ["m-n101-k10", "tai150a", "cmt10"]
    methods   = [:cw, :regret₂init, :regret₃init]
    χₒ  = ObjectiveFunctionParameters(
        d = 0.                          ,
        v = 10000.                      ,
        r = 0.                          ,
        c = 0.                          ,
    )
    χ   = ALNSParameters(
        k̲   =   2                       ,
        k̅   =   500                     ,
        k̲ₛ  =   80                      ,
        k̅ₛ  =   250                     ,   
        Ψᵣ  =   [
                    :randomnode!    , 
                    :relatedpair!   ,
                    :relatednode!   , 
                    :worstnode!     ,  
                    :randomroute!   ,
                    :relatedroute!  ,
                    :worstroute!    ,
                    :randomvehicle! ,  
                ]                       , 
        Ψᵢ  =   [
                    :best!          ,
                    :greedy!        ,
                    :regret₂insert! ,
                    :regret₃insert!
                ]                       ,
        Ψₛ  =   [
                    :move!              ,
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
        ρ   =   0.1                     ,
        χₒ  =   χₒ  
    )
    for k ∈ 1:K
        instance = instances[k]
        method = methods[k]
        println("\n Solving $instance")
        G = build(instance)
        sₒ= initialsolution(G, χₒ, method)     
        @test isfeasible(sₒ)
        S = ALNS(sₒ, χ)
        s⃰ = S[end]
        @test isfeasible(s⃰)
        @test f(s⃰, χₒ) ≤ f(sₒ, χₒ)
    end
    return
end