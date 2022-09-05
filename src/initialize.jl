# Initial solution
"""
    initialsolution([rng], instance, method)

Returns initial LRP solution for the given `instance` using the given `method`.

Available methods include,
- Random Initialization                 : `:random`

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
initialsolution(rng::AbstractRNG, instance, method::Symbol)::Solution = getfield(LRP, method)(rng, instance)
initialsolution(instance, method::Symbol) = initialsolution(Random.GLOBAL_RNG, instance, method)

# Random Initialization
# Create initial solution with randomly selcted node-route combination until all customer nodes have been added to the solution
function random(rng::AbstractRNG, instance)
    G = build(instance)
    s = Solution(G...)
    D = s.D
    C = s.C
    for d ∈ D for v ∈ d.V push!(v.R, Route(v, d)) end end
    # Step 1: Initialize
    w = ones(Int64, eachindex(C))                      # w[i]: selection weight for customer node C[i]
    # Step 2: Iteratively append randomly selected customer node in randomly selected route
    while any(isopen, C)
        c = sample(rng, C, OffsetWeights(w))
        d = sample(rng, D)
        v = sample(rng, d.V)
        r = sample(rng, v.R)
        nᵗ = d
        nʰ = isopt(r) ? C[r.iˢ] : D[r.iˢ]
        insertnode!(c, nᵗ, nʰ, r, s)
        iⁿ = c.iⁿ
        w[iⁿ] = 0
    end
    # Step 3: Return initial solution
    for d ∈ D 
        for v ∈ d.V
            deleteat!(v.R, deleteroute.(v.R))
            for (iʳ,r) ∈ pairs(v.R) r.iʳ = iʳ end
        end
    end
    return s
end
