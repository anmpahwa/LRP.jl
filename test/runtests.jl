using LRP
using Revise
using Test
using Random

@testset "TSP.jl" begin
    K = 5
    instances = ["att48", "eil101", "ch150", "d198", "a280"]
    methods = [:cw_init, :nn_init, :random_init, :regret₂init, :regret₃init]
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
                    :node_remove! , 
                    :worst_remove!  , 
                    :shaw_remove!
                ]                       , 
        Ψᵢ  =   [
                    :best_insert!   ,
                    :greedy_insert! ,
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
