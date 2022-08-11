"""
    remove!([rng], q, s::Solution, χₒ::ObjectiveFunctionParameters, method::Symbol)

Return solution removing q nodes from solution s using the given `method`.
`χₒ` includes the objective function parameters for objective function evaluation.

Available methods include,
- Random Node Removal   : `node_remove!`
- Related Node Removal  : `shaw_remove!`
- Worst Node Removal    : `worst_remove!`

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
remove!(rng::AbstractRNG, q, s::Solution, χₒ::ObjectiveFunctionParameters, method::Symbol)::Solution = getfield(LRP, method)(rng, q, s, χₒ)
remove!(q, s::Solution, χₒ::ObjectiveFunctionParameters, method::Symbol) = remove!(Random.GLOBAL_RNG, q, s, χₒ, method)

# Random Node Removal
# Randomly select q customer nodes to remove
function node_remove!(rng::AbstractRNG, q, s::Solution)
    D = s.D
    C = s.C
    # Step 1: Randomly select customer nodes to remove until q customer nodes have been removed
    n = 0
    while n < q
        c  = rand(rng, C)
        if isopen(c) continue end
        r  = c.r
        nₜ = isequal(r.s, c.i) ? D[c.t] : C[c.t]
        nₕ = isequal(r.e, c.i) ? D[c.h] : C[c.h]
        removenode!(c, nₜ, nₕ, r, s)
        n += 1
    end
    # Step 2: Return solution
    return s
end
node_remove!(rng::AbstractRNG, q, s::Solution, χₒ::ObjectiveFunctionParameters) = node_remove!(rng, q, s)

# Related Node Remove
# For a randomly selected customer node, remove q most related customer nodes
function shaw_remove!(rng::AbstractRNG, q, s::Solution)
    D = s.D
    C = s.C
    A = s.A
    V = s.V
    # Step 1: Randomly select a pivot customer node
    j  = rand(rng, eachindex(C))
    cₒ = C[j]
    # Step 2: For each customer node, evaluate relatedness to this pivot customer node
    x = fill(-Inf, eachindex(C))    # x[k]: relatedness of customer node C[i] with customer node C[j]  
    for (i,c) ∈ pairs(C)
        r = c.r
        a = A[(i,j)]
        v = V[r.o] 
        x[i] = 1/(0.482 * a.l + 35.0 * a.t + 3.826 * a.f - 1.0 * (isequal(c.q, cₒ.q)) - 1.0 * (isequal(c.r, cₒ.r)))
    end
    # Step 3: Remove q most related customer nodes
    for _ ∈ 1:q
        k  = argmax(x)
        c  = C[k]
        r  = c.r
        nₜ = isequal(r.s, c.i) ? D[c.t] : C[c.t]
        nₕ = isequal(r.e, c.i) ? D[c.h] : C[c.h]
        removenode!(c, nₜ, nₕ, r, s)
        x[k] = -Inf
    end
    # Step 4: Return solution
    return s
end
shaw_remove!(rng::AbstractRNG, q, s::Solution, χₒ::ObjectiveFunctionParameters) = shaw_remove!(rng, q, s)

# Worst Node Remove
# Remove q customer nodes with highest removal cost
function worst_remove!(q, s::Solution, χₒ::ObjectiveFunctionParameters)
    D = s.D
    C = s.C
    V = s.V
    x = fill(Inf, eachindex(C))     # x[k]: removal cost of customer node C[k]
    ϕ = ones(Int64, length(V))      # ϕ[k]: selection weight for vehicle V[k]
    # Step 1: Iterate until q customer nodes have been removed
    for _ ∈ 1:q
        # Step 1.1: For every closed customer node evaluate removal cost
        z = f(s, χₒ)
        for (k,c) ∈ pairs(C)
            if isopen(c) continue end
            if iszero(ϕ[c.r.o]) continue end
            # Step 1.1.1: Remove closed customer node c between tail node nₜ and head node nₕ in route r
            r  = c.r
            nₜ = isequal(r.s, c.i) ? D[c.t] : C[c.t]
            nₕ = isequal(r.e, c.i) ? D[c.h] : C[c.h]
            removenode!(c, nₜ, nₕ, r, s)
            # Step 1.1.2: Evaluate the removal cost
            z⁻ = f(s, χₒ)
            x[k] = z - z⁻
            # Step 1.1.3: Re-insert customer node c between tail node nₜ and head node nₕ in route r
            insertnode!(c, nₜ, nₕ, r, s)
        end
        # Step 1.2: Remove the customer node with highest removal cost
        k  = argmax(x)
        c  = C[k]
        r  = c.r
        nₜ = isequal(r.s, c.i) ? D[c.t] : C[c.t]
        nₕ = isequal(r.e, c.i) ? D[c.h] : C[c.h]
        removenode!(c, nₜ, nₕ, r, s)
        # Step 1.3: Update cost and selection weight vectors
        x[k] = -Inf
        for (k,v) ∈ pairs(V) ϕ[k] = isequal(r.o, v.i) ? 1 : 0 end
    end
    # Step 2: Return solution
    return s
end
worst_remove!(rng::AbstractRNG, q, s::Solution, χₒ::ObjectiveFunctionParameters) = worst_remove!(q, s, χₒ)