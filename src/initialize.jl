# Initial solution
"""
    initialsolution([rng], instance, method)

Returns initial LRP solution for the given `instance` using the given `method`.

Available methods include,
- Random Initialization             : `:random`
- K-means Clustering Initialization : `:cluster`

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
    # Step 1: Initialize
    preinitialize!(s)
    d = sample(rng, D)
    v = sample(rng, d.V)
    r = sample(rng, v.R)
    w = ones(Int64, eachindex(C))                      # w[i]: selection weight for customer node C[i]
    # Step 2: Iteratively append randomly selected customer node in randomly selected route
    while any(isopen, C)
        c = sample(rng, C, OffsetWeights(w))
        nᵗ = d
        nʰ = isopt(r) ? C[r.iˢ] : D[r.iˢ]
        insertnode!(c, nᵗ, nʰ, r, s)
        iⁿ = c.iⁿ
        w[iⁿ] = 0
    end
    postinitialize!(s)
    # Step 4: Return initial solution
    return s
end

# k-means clustering Initialization
# Create initial solution using k-means clustering algorithm
function cluster(rng::AbstractRNG, instance)
    G = build(instance)
    s = Solution(G...)
    D = s.D
    C = s.C
    V = [v for d ∈ D for v ∈ d.V]
    # Step 1: Initialize
    preinitialize!(s)
    X = zeros(4, eachindex(C))
    for (iⁿ,c) ∈ pairs(C) X[:,iⁿ] = [c.x, c.y, c.tᵉ, c.tˡ] end
    # Step 2: Clustering
    k = ceil(Int64, sum(getproperty.(C, :q))/mean(getproperty.(V, :q)))
    Y = kmeans(X.parent, k)
    # Step 3: Assignment
    A  = OffsetVector(Y.assignments, eachindex(C))
    Cᵒ = Y.centers
    Iᵒ = 1:(size(Cᵒ)[2])
    w  = ones(Int64, length(D))
    for iᵒ ∈ Iᵒ
        Z = fill(Inf, length(D))
        for d ∈ D
            iⁿ = d.iⁿ
            if iszero(w[iⁿ]) continue end
            xᵒ = Cᵒ[1,iᵒ]
            yᵒ = Cᵒ[2,iᵒ]
            xᵈ = d.x
            yᵈ = d.y
            Z[iⁿ] = sqrt((xᵒ-xᵈ)^2 + (yᵒ-yᵈ)^2)
        end
        iⁿ = argmin(Z)
        d  = D[iⁿ]
        v  = sample(rng, d.V, Weights((!isopt).(d.V)))
        r  = sample(rng, v.R)
        nᵗ = d
        Cᶜ = filter(c -> isequal(A[c.iⁿ], iᵒ), C)
        for c ∈ Cᶜ
            insertnode!(c, nᵗ, d, r, s)
            nᵗ = c
        end
        if all(isopt, d.V) w[iⁿ] = 0 end
    end
    postinitialize!(s)
    # Step 4: Return initial solution
    return s
end