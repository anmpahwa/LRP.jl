"""
    insert!([rng], s::Solution, method::Symbol)

Return solution inserting open customer nodes to the solution `s` using the given `method`.

Available methods include,
- Precise Best Insertion    : `:bestprecise!`
- Perturb Best Insertion    : `:bestperturb!`
- Precise Greedy Insertion  : `:greedyprecise!`
- Perturb Greedy Insertion  : `:greedyperturb!`
- Regret-two Insertion      : `:regret2!`
- Regret-three Insertion    : `:regret3!`

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
insert!(rng::AbstractRNG, s::Solution, method::Symbol)::Solution = getfield(LRP, method)(rng, s)
insert!(s::Solution, method::Symbol) = insert!(Random.GLOBAL_RNG, s, method)

# Best insertion
# Iteratively insert randomly selected customer node at its best position until all open customer nodes have been added to the solution
function bestinsert!(rng::AbstractRNG, s::Solution; noise=false)
    if all(isclose, s.C) return s end
    D  = s.D
    C  = s.C
    ϕⁿ = noise
    # Step 1: Initialize
    preinsertion!(s)
    R = [r for d ∈ D for v ∈ d.V for r ∈ v.R]
    L = [c for c ∈ C if isopen(c)]
    I = eachindex(L)
    J = eachindex(R)
    X = ElasticMatrix(fill(Inf, (I,J)))     # X[i,j]: insertion cost of customer node L[i] at best position in route R[j]
    P = ElasticMatrix(fill((0, 0), (I,J)))  # P[i,j]: best insertion postion of customer node L[i] in route R[j]
    W = ones(Int64, I)                      # W[j]  : selection weight for customer node L[i]
    ϕ = ones(Int64, J)                      # ϕ[j]  : selection weight for route R[j]
    # Step 2: Iterate until all open customer nodes have been inserted into the route
    for _ ∈ I
        # Step 2.1: Iterate through all open customer nodes and every possible insertion position in each route
        zᵒ = f(s)
        for (i,c) ∈ pairs(L)
            if !isopen(c) continue end
            for (j,r) ∈ pairs(R)
                if iszero(ϕ[j]) continue end
                d  = s.D[r.iᵈ]
                nˢ = isopt(r) ? C[r.iˢ] : D[r.iˢ]
                nᵉ = isopt(r) ? C[r.iᵉ] : D[r.iᵉ]
                nᵗ = d
                nʰ = nˢ
                while true
                    # Step 2.1.1: Insert customer node c between tail node nᵗ and head node nʰ in route r
                    insertnode!(c, nᵗ, nʰ, r, s)
                    # Step 2.1.2: Compute the insertion cost
                    z⁺ = f(s) * (1 + ϕⁿ * rand(rng, Uniform(-0.2, 0.2)))
                    Δ  = z⁺ - zᵒ
                    # Step 2.1.3: Revise least insertion cost in route r and the corresponding best insertion position in route r
                    if Δ < X[i,j] X[i,j], P[i,j] = Δ, (nᵗ.iⁿ, nʰ.iⁿ) end
                    # Step 2.1.4: Remove customer node c from its position between tail node nᵗ and head node nʰ
                    removenode!(c, nᵗ, nʰ, r, s)
                    if isequal(nᵗ, nᵉ) break end
                    nᵗ = nʰ
                    nʰ = isequal(r.iᵉ, nᵗ.iⁿ) ? D[nᵗ.iʰ] : C[nᵗ.iʰ]
                end
            end
        end
        # Step 2.2: Randomly select a customer node to insert at its best position
        i = sample(rng, I, Weights(W))
        j = argmin(X[i,:])
        c = L[i]
        r = R[j]
        d = s.D[r.iᵈ]
        v = d.V[r.iᵛ]
        iᵗ = P[i,j][1]
        iʰ = P[i,j][2]
        nᵗ = iᵗ ≤ length(D) ? D[iᵗ] : C[iᵗ]
        nʰ = iʰ ≤ length(D) ? D[iʰ] : C[iʰ]
        insertnode!(c, nᵗ, nʰ, r, s)
        # Step 2.3: Revise vectors appropriately
        X[i,:] .= Inf
        P[i,:] .= ((0, 0), )
        W[i] = 0
        ϕ .= 0
        for (j,r) ∈ pairs(R) 
            if !isequal(r.iᵛ, v.iᵛ) continue end
            X[:,j] .= Inf
            P[:,j] .= ((0, 0), )
            ϕ[j] = 1  
        end
        if addroute(r, s)
            r = Route(v, d)
            push!(v.R, r) 
            push!(R, r)
            append!(X, fill(Inf, (I,1)))
            append!(P, fill((0, 0), (I,1)))
            push!(ϕ, 1)
        end
        if addvehicle(v, s)
            v = Vehicle(v, d)
            r = Route(v, d)
            push!(d.V, v)
            push!(v.R, r) 
            push!(R, r)
            append!(X, fill(Inf, (I,1)))
            append!(P, fill((0, 0), (I,1)))
            push!(ϕ, 1)
        end
    end
    postinsertion!(s)
    # Step 3: Return solution
    return s
end
bestprecise!(rng::AbstractRNG, s::Solution) = bestinsert!(rng, s; noise=false)
bestperturb!(rng::AbstractRNG, s::Solution) = bestinsert!(rng, s; noise=true)

# Greedy insertion
# Iteratively insert customer nodes with least insertion cost until all open customer nodes have been added to the solution
function greedyinsert!(rng::AbstractRNG, s::Solution; noise=false)
    if all(isclose, s.C) return s end
    D  = s.D
    C  = s.C
    ϕⁿ = noise
    # Step 1: Initialize
    preinsertion!(s)
    R = [r for d ∈ D for v ∈ d.V for r ∈ v.R]
    L = [c for c ∈ C if isopen(c)]
    I = eachindex(L)
    J = eachindex(R)
    X = ElasticMatrix(fill(Inf, (I,J)))     # X[i,j]: insertion cost of customer node L[i] at best position in route R[j]
    P = ElasticMatrix(fill((0, 0), (I,J)))  # P[i,j]: best insertion postion of customer node L[i] in route R[j]
    ϕ = ones(Int64, J)                      # ϕ[j]  : selection weight for route R[j]
    # Step 2: Iterate until all open customer nodes have been inserted into the route
    for _ ∈ I
        # Step 2.1: Iterate through all open customer nodes and every possible insertion position in each route
        zᵒ = f(s)
        for (i,c) ∈ pairs(L)
            if !isopen(c) continue end
            for (j,r) ∈ pairs(R)
                if iszero(ϕ[j]) continue end
                d  = s.D[r.iᵈ]
                nˢ = isopt(r) ? C[r.iˢ] : D[r.iˢ]
                nᵉ = isopt(r) ? C[r.iᵉ] : D[r.iᵉ]
                nᵗ = d
                nʰ = nˢ
                while true
                    # Step 2.1.1: Insert customer node c between tail node nᵗ and head node nʰ in route r
                    insertnode!(c, nᵗ, nʰ, r, s)
                    # Step 2.1.2: Compute the insertion cost
                    z⁺ = f(s) * (1 + ϕⁿ * rand(rng, Uniform(-0.2, 0.2)))
                    Δ  = z⁺ - zᵒ
                    # Step 2.1.3: Revise least insertion cost in route r and the corresponding best insertion position in route r
                    if Δ < X[i,j] X[i,j], P[i,j] = Δ, (nᵗ.iⁿ, nʰ.iⁿ) end
                    # Step 2.1.4: Remove customer node c from its position between tail node nᵗ and head node nʰ
                    removenode!(c, nᵗ, nʰ, r, s)
                    if isequal(nᵗ, nᵉ) break end
                    nᵗ = nʰ
                    nʰ = isequal(r.iᵉ, nᵗ.iⁿ) ? D[nᵗ.iʰ] : C[nᵗ.iʰ]
                end
            end
        end
        # Step 2.2: Randomly select a customer node to insert at its best position        
        i,j = Tuple(argmin(X))
        c = L[i]
        r = R[j]
        d = s.D[r.iᵈ]
        v = d.V[r.iᵛ]
        iᵗ = P[i,j][1]
        iʰ = P[i,j][2]
        nᵗ = iᵗ ≤ length(D) ? D[iᵗ] : C[iᵗ]
        nʰ = iʰ ≤ length(D) ? D[iʰ] : C[iʰ]
        insertnode!(c, nᵗ, nʰ, r, s)
        # Step 2.3: Revise vectors appropriately
        X[i,:] .= Inf
        P[i,:] .= ((0, 0), )
        ϕ .= 0
        for (j,r) ∈ pairs(R) 
            if !isequal(r.iᵛ, v.iᵛ) continue end
            X[:,j] .= Inf
            P[:,j] .= ((0, 0), )
            ϕ[j] = 1  
        end
        if addroute(r, s)
            r = Route(v, d)
            push!(v.R, r) 
            push!(R, r)
            append!(X, fill(Inf, (I,1)))
            append!(P, fill((0, 0), (I,1)))
            push!(ϕ, 1)
        end
        if addvehicle(v, s)
            v = Vehicle(v, d)
            r = Route(v, d)
            push!(d.V, v)
            push!(v.R, r) 
            push!(R, r)
            append!(X, fill(Inf, (I,1)))
            append!(P, fill((0, 0), (I,1)))
            push!(ϕ, 1)
        end
    end
    postinsertion!(s)
    # Step 3: Return solution
    return s
end
greedyprecise!(rng::AbstractRNG, s::Solution) = greedyinsert!(rng, s; noise=false)
greedyperturb!(rng::AbstractRNG, s::Solution) = greedyinsert!(rng, s; noise=true)

# Regret-N insertion
# Iteratively add customer nodes with highest regret cost at its best position until all open customer nodes have been added to the solution
function regretNinsert!(rng::AbstractRNG, N::Int64, s::Solution)
    if all(isclose, s.C) return s end
    D = s.D
    C = s.C
    # Step 1: Initialize
    preinsertion!(s)
    R = [r for d ∈ D for v ∈ d.V for r ∈ v.R]
    L = [c for c ∈ C if isopen(c)]
    I = eachindex(L)
    J = eachindex(R)
    X = ElasticMatrix(fill(Inf, (I,J)))     # X[i,j]: insertion cost of customer node L[i] at best position in route R[j]
    P = ElasticMatrix(fill((0, 0), (I,J)))  # P[i,j]: best insertion postion of customer node L[i] in route R[j]
    Y = fill(Inf, (I,N))                    # Y[i,n]: insertion cost of customer node L[i] at nᵗʰ best position
    W = zeros(Int64, (I,N))                 # W[i,n]: route index of customer node L[j] at nᵗʰ best position
    Z = fill(-Inf, I)                       # Z[i]  : regret-N cost of customer node L[i]
    ϕ = ones(Int64, J)                      # ϕ[j]  : selection weight for route R[j]
    # Step 2: Iterate until all open customer nodes have been inserted into the route
    for _ ∈ I
        # Step 2.1: Iterate through all open customer nodes and every route
        zᵒ = f(s)
        for (i,c) ∈ pairs(L)
            if !isopen(c) continue end
            for (j,r) ∈ pairs(R)
                # Step 2.1.1: Iterate through all possible insertion position in route r
                if iszero(ϕ[j]) continue end
                d  = s.D[r.iᵈ]
                nˢ = isopt(r) ? C[r.iˢ] : D[r.iˢ]
                nᵉ = isopt(r) ? C[r.iᵉ] : D[r.iᵉ]
                nᵗ = d
                nʰ = nˢ
                while true
                    # Step 2.1.1.1: Insert customer node c between tail node nᵗ and head node nʰ in route r
                    insertnode!(c, nᵗ, nʰ, r, s)
                    # Step 2.1.1.2: Compute the insertion cost
                    z⁺ = f(s)
                    Δ  = z⁺ - zᵒ
                    # Step 2.1.1.3: Revise least insertion cost in route r and the corresponding best insertion position in route r
                    if Δ < X[i,j] X[i,j], P[i,j] = Δ, (nᵗ.iⁿ, nʰ.iⁿ) end
                    # Step 2.1.1.4: Revise N least insertion costs
                    n̲ = 1
                    for n ∈ 1:N 
                        n̲ = n
                        if Δ < Y[i,n] break end
                    end
                    for n ∈ N:-1:n̲ 
                        Y[i,n] = isequal(n, n̲) ? Δ : Y[i,n-1]::Float64
                        W[i,n] = isequal(n, n̲) ? r.iʳ : W[i,n-1]::Int64
                    end
                    # Step 2.1.1.5: Remove customer node c from its position between tail node nᵗ and head node nʰ in route r
                    removenode!(c, nᵗ, nʰ, r, s)
                    if isequal(nᵗ, nᵉ) break end
                    nᵗ = nʰ
                    nʰ = isequal(r.iᵉ, nᵗ.iⁿ) ? D[nᵗ.iʰ] : C[nᵗ.iʰ]
                end
            end
            # Step 2.1.2: Compute regret cost for customer node c
            Z[i] = 0.
            for n ∈ 1:N Z[i] += Y[i,n] - Y[i,1] end
        end
        # Step 2.2: Insert customer node with highest regret cost in its best position (break ties by inserting the node with the lowest insertion cost)
        I̲ = findall(i -> i == maximum(Z), Z)
        i,j = Tuple(argmin(X[I̲,:]))
        i = I̲[i]
        c = L[i]
        r = R[j]
        d = s.D[r.iᵈ]
        v = d.V[r.iᵛ]
        iᵗ = P[i,j][1]
        iʰ = P[i,j][2]
        nᵗ = iᵗ ≤ length(D) ? D[iᵗ] : C[iᵗ]
        nʰ = iʰ ≤ length(D) ? D[iʰ] : C[iʰ]
        insertnode!(c, nᵗ, nʰ, r, s)
        # Step 2.3: Revise vectors appropriately
        X[i,:] .= Inf
        P[i,:] .= ((0, 0), )
        Y[i,:] .= Inf
        W[i,:] .= 0
        Z .= -Inf 
        for (i,c) ∈ pairs(L)
            for n ∈ 1:N
                if iszero(w[i,n]) break end
                k = findfirst(x -> x.iʳ == W[i,n], R)
                r = R[k]
                if isequal(r.iᵛ, v.iᵛ) Y[i,n], W[i,n] = Inf, 0 end
            end
            ix = sortperm(y[i,:])
            Y[i,:] .= y[i,ix]
            W[i,:] .= w[i,ix]
        end
        ϕ .= 0
        for (j,r) ∈ pairs(R) 
            if !isequal(r.iᵛ, v.iᵛ) continue end
            X[:,j] .= Inf
            P[:,j] .= ((0, 0), )
            ϕ[j] = 1  
        end
        if addroute(r, s)
            r = Route(v, d)
            push!(v.R, r) 
            push!(R, r)
            append!(X, fill(Inf, (I,1)))
            append!(P, fill((0, 0), (I,1)))
            push!(ϕ, 1)
        end
        if addvehicle(v, s)
            v = Vehicle(v, d)
            r = Route(v, d)
            push!(d.V, v)
            push!(v.R, r) 
            push!(R, r)
            append!(X, fill(Inf, (I,1)))
            append!(P, fill((0, 0), (I,1)))
            push!(ϕ, 1)
        end
    end
    postinsertion!(s)
    # Step 3: Return solution
    return s
end
regret2!(rng::AbstractRNG, s::Solution) = regretNinsert!(rng, Int64(2), s)
regret3!(rng::AbstractRNG, s::Solution) = regretNinsert!(rng, Int64(3), s)