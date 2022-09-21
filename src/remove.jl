"""
    remove!([rng], q::Int64, s::Solution, method::Symbol)

Return solution removing q customer nodes from solution s using the given `method`.

Available methods include,
- Random Node Removal       : `:randomnode!`
- Random Route Removal      : `:randomroute!`
- Random Vehicle Removal    : `:randomvehicle!`
- Random Depot Removal      : `:randomdepot!` 
- Related Node Removal      : `:relatednode!`
- Related Route removal     : `:relatedroute!`
- Related Vehicle Removal   : `:relatedvehicle!`
- Related Depot Removal     : `:relateddepot!`
- Worst Node Removal        : `:worstnode!`
- Worst Route Removal       : `:worstroute!`
- Worst Vehicle Removal     : `:worstvehicle!`
- Worst Depot Removal       : `:worstdepot!`

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
remove!(rng::AbstractRNG, q::Int64, s::Solution, method::Symbol)::Solution = getfield(LRP, method)(rng, q, s)
remove!(q::Int64, s::Solution, method::Symbol) = remove!(Random.GLOBAL_RNG, q, s, method)

# -------------------------------------------------- NODE REMOVAL --------------------------------------------------
# Random Node Removal
# Randomly select q customer nodes to remove
function randomnode!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    preremoval(s)
    # Step 1: Randomly select customer nodes to remove until q customer nodes have been removed
    n = 0
    w = isclose.(C)
    while n < q
        iⁿ = sample(rng, eachindex(C), OffsetWeights(w))
        c  = C[iⁿ]
        if isopen(c) continue end
        r  = c.r
        nᵗ = isequal(r.iˢ, c.iⁿ) ? D[c.iᵗ] : C[c.iᵗ]
        nʰ = isequal(r.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ]
        removenode!(c, nᵗ, nʰ, r, s)
        n += 1
        w[iⁿ] = 0
    end
    postremoval(s)
    # Step 2: Return solution
    return s
end

# Related Node Removal (related to pivot)
# For a randomly selected customer node, remove q most related customer nodes
function relatednode!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    A = s.A
    preremoval(s)
    # Step 1: Randomly select a pivot customer node
    iᵒ = rand(rng, eachindex(C))
    # Step 2: For each customer node, evaluate relatedness to this pivot customer node
    x  = fill(-Inf, eachindex(C))   # x[iⁿ]: relatedness of customer node C[iⁿ] with customer node C[iᵒ]  
    for iⁿ ∈ eachindex(C) x[iⁿ] = relatedness(C[iⁿ], C[iᵒ], A[(iⁿ,iᵒ)]) end
    # Step 3: Remove q most related customer nodes
    n = 0
    while n < q
        iⁿ = argmax(x)
        c  = C[iⁿ]
        r  = c.r
        nᵗ = isequal(r.iˢ, c.iⁿ) ? D[c.iᵗ] : C[c.iᵗ]
        nʰ = isequal(r.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ]
        removenode!(c, nᵗ, nʰ, r, s)
        n += 1
        x[iⁿ] = -Inf
    end
    # Step 4: Remove redundant vehicles and routes
    postremoval(s)
    # Step 5: Return solution
    return s
end

# Worst Node Removal
# Remove q customer nodes with highest removal cost (savings)
function worstnode!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    preremoval(s)
    V = [v for d ∈ D for v ∈ d.V]
    x = fill(-Inf, eachindex(C))    # x[i]: removal cost of customer node C[i]
    ϕ = ones(Int64, eachindex(V))   # ϕ[j]: selection weight for vehicle V[j]
    # Step 1: Iterate until q customer nodes have been removed
    n = 0
    while n < q
        # Step 1.1: For every closed customer node evaluate removal cost
        zᵒ = f(s)
        for (i,c) ∈ pairs(C)
            if isopen(c) continue end
            r = c.r
            d = s.D[r.iᵈ]
            v = d.V[r.iᵛ]
            j = findfirst(isequal(v), V)
            if iszero(ϕ[j]) continue end
            # Step 1.1.1: Remove closed customer node c between tail node nᵗ and head node nʰ in route r
            nᵗ = isequal(r.iˢ, c.iⁿ) ? D[c.iᵗ] : C[c.iᵗ]
            nʰ = isequal(r.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ]
            removenode!(c, nᵗ, nʰ, r, s)
            # Step 1.1.2: Evaluate the removal cost
            z⁻ = f(s)
            Δ  = z⁻ - zᵒ
            x[i] = -Δ
            # Step 1.1.3: Re-insert customer node c between tail node nᵗ and head node nʰ in route r
            insertnode!(c, nᵗ, nʰ, r, s)
        end
        # Step 1.2: Remove the customer node with highest removal cost (savings)
        iⁿ = argmax(x)
        c  = C[iⁿ]
        r  = c.r
        d  = s.D[r.iᵈ]
        v  = d.V[r.iᵛ]
        nᵗ = isequal(r.iˢ, c.iⁿ) ? D[c.iᵗ] : C[c.iᵗ]
        nʰ = isequal(r.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ]
        removenode!(c, nᵗ, nʰ, r, s)
        n += 1
        # Step 1.3: Update cost and selection weight vectors
        x[iⁿ] = -Inf
        for (j,v) ∈ pairs(V) ϕ[j] = isequal(r.iᵛ, v.iᵛ) ? 1 : 0 end 
    end
    postremoval(s)
    # Step 2: Return solution
    return s
end

# -------------------------------------------------- ROUTE REMOVAL --------------------------------------------------
# Random Route Removal
# Iteratively select a random route and remove customer nodes from it until at least q customer nodes are removed
function randomroute!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    preremoval(s)
    R = [r for d ∈ D for v ∈ d.V for r ∈ v.R]
    # Step 1: Iteratively select a random route and remove customer nodes from it until at least q customer nodes are removed
    n = 0
    w = isopt.(R)
    while n < q
        if isone(sum(w)) break end
        iʳ = sample(rng, eachindex(R), Weights(w))
        r  = R[iʳ]
        d  = D[r.iᵈ]
        while true
            if n ≥ q break end
            nᵗ = d
            c  = C[r.iˢ]
            nʰ = isequal(r.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ]
            removenode!(c, nᵗ, nʰ, r, s)
            n += 1
            if isequal(nʰ, d) break end
        end
        w[iʳ] = 0
    end
    postremoval(s)
    # Step 2: Return solution
    return s
end

# Related Route Removal
# For a randomly selected route, remove customer nodes from most related route until q customer nodes are removed
function relatedroute!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    preremoval(s)
    R = [r for d ∈ D for v ∈ d.V for r ∈ v.R]
    # Step 1: Randomly select a pivot route
    iᵒ = sample(rng, eachindex(R), Weights(isopt.(R)))  
    # Step 2: For each route, evaluate relatedness to this pivot route
    x  = fill(-Inf, eachindex(R))
    for iʳ ∈ eachindex(R) x[iʳ] = relatedness(R[iʳ], R[iᵒ]) end
    # Step 3: Remove at least q customers from most related route to this pivot route
    n = 0
    w = isopt.(R)
    while n < q
        if isone(sum(w)) break end
        iʳ = argmax(x)
        r  = R[iʳ]
        d  = D[r.iᵈ]
        while true
            if n ≥ q break end
            nᵗ = d
            c  = C[r.iˢ]
            nʰ = isequal(r.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ]
            removenode!(c, nᵗ, nʰ, r, s)
            n += 1
            if isequal(nʰ, d) break end
        end 
        x[iʳ] = -Inf
        w[iʳ] = 0
    end
    postremoval(s)
    # Step 4: Return solution
    return s
end

# Worst Route Removal
# Iteratively select low-utilization route and remove customer nodes from it until at least q customer nodes are removed
function worstroute!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    preremoval(s)
    R = [r for d ∈ D for v ∈ d.V for r ∈ v.R]
    # Step 1: Evaluate utilization of each route
    x = fill(Inf, eachindex(R))
    for (iʳ,r) ∈ pairs(R)
        if !isopt(r) continue end
        d = s.D[r.iᵈ]
        v = d.V[r.iᵛ]
        x[iʳ] = r.q/v.q
    end
    # Step 2: Iteratively select low-utilization route and remove customer nodes from it until at least q customer nodes are removed
    n = 0
    w = isopt.(R)
    while n < q
        if isone(sum(w)) break end
        iʳ = argmin(x)
        r  = R[iʳ]
        d  = D[r.iᵈ]
        while true
            if n ≥ q break end
            nᵗ = d
            c  = C[r.iˢ]
            nʰ = isequal(r.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ]
            removenode!(c, nᵗ, nʰ, r, s)
            n += 1
            if isequal(nʰ, d) break end
        end
        x[iʳ] = Inf
        w[iʳ] = 0
    end
    postremoval(s)
    # Step 3: Return solution
    return s
end
    
# -------------------------------------------------- VEHICLE REMOVAL --------------------------------------------------
# Random Vehicle Removal
# Iteratively select a random vehicle and remove customer nodes from it until at least q customer nodes are removed
function randomvehicle!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    preremoval(s)
    V = [v for d ∈ D for v ∈ d.V]
    # Step 1: Iteratively select a random vehicle and remove customer nodes from it until at least q customer nodes are removed
    n = 0
    w = ones(Int64, eachindex(V))
    while n < q
        if isone(sum(w)) break end
        iᵛ = sample(rng, eachindex(V), Weights(w))
        v  = V[iᵛ]
        d  = D[v.iᵈ]
        for r ∈ v.R
            if n ≥ q break end
            if !isopt(r) continue end
            while true
                nᵗ = d
                c  = C[r.iˢ]
                nʰ = isequal(r.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ] 
                removenode!(c, nᵗ, nʰ, r, s)
                n += 1
                if isequal(nʰ, d) break end
            end
        end
        w[iᵛ] = 0
    end
    postremoval(s)
    # Step 2: Return solution
    return s
end

# Related Vehicle Removal
# For a randomly selected vehicle, remove customer nodes from most related vehicles until q customer nodes are removed
function relatedvehicle!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    preremoval(s)
    V = [v for d ∈ D for v ∈ d.V]
    # Step 1: Select a random closed depot node
    iᵒ = sample(rng, eachindex(V), Weights(isopt.(V)))
    # Step 2: For each vehicle, evaluate relatedness to this pivot vehicle
    x  = fill(-Inf, eachindex(V))
    for iᵛ ∈ eachindex(V) x[iᵛ] = relatedness(V[iᵛ], V[iᵒ]) end
    # Step 3: Remove at least q customers from the most related vehicles to this pivot vehicle
    n = 0
    w = ones(Int64, eachindex(V))
    while n < q
        if isone(sum(w)) break end
        iᵛ = argmax(x)
        v  = V[iᵛ]
        d  = D[v.iᵈ] 
        for r ∈ v.R
            if n ≥ q break end
            if !isopt(r) continue end
            while true
                nᵗ = d
                c  = C[r.iˢ]
                nʰ = isequal(r.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ]
                removenode!(c, nᵗ, nʰ, r, s)
                n += 1
                if isequal(nʰ, d) break end
            end
        end
        x[iᵛ] = -Inf
        w[iᵛ] = 0
    end
    postremoval(s)
    # Step 4: Return solution
    return s
end

# Worst Vehicle Removal
# Iteratively select low-utilization vehicle and remove customer nodes from it until at least q customer nodes are removed
function worstvehicle!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    preremoval(s)
    V = [v for d ∈ D for v ∈ d.V]
    # Step 1: Evaluate utilization for each vehicle
    x = fill(Inf, eachindex(V))
    for (iᵛ,v) ∈ pairs(V)
        a = 0
        b = 0
        for r ∈ v.R
            a += r.q
            b += v.q
        end
        x[iᵛ] = a/b
    end
    # Step 2: Iteratively select low-utilization route and remove customer nodes from it until at least q customer nodes are removed
    n = 0
    w = ones(Int64, eachindex(V))
    while n < q
        if isone(sum(w)) break end
        iᵛ = argmin(x)
        v  = V[iᵛ]
        d  = D[v.iᵈ]
        for r ∈ v.R
            if n ≥ q break end
            if !isopt(r) continue end
            while true
                nᵗ = d
                c  = C[r.iˢ]
                nʰ = isequal(r.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ]
                removenode!(c, nᵗ, nʰ, r, s)
                n += 1
                if isequal(nʰ, d) break end
            end
        end
        x[iᵛ] = Inf
        w[iᵛ] = 0
    end
    postremoval(s)
    # Step 3: Return solution
    return s
end

# -------------------------------------------------- DEPOT REMOVAL --------------------------------------------------
# Random Depot Removal
# Iteratively select a random depot and remove customer nodes from it until at least q customer nodes are removed
function randomdepot!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    preremoval(s)
    # Step 1: Iteratively select a random depot and remove customer nodes from it until at least q customer nodes are removed
    n = 0
    w = ones(Int64, eachindex(D))
    while n < q
        if isone(sum(w)) break end
        iᵈ = sample(rng, eachindex(D), Weights(w))
        d  = D[iᵈ]
        for v ∈ d.V
            if n ≥ q break end
            for r ∈ v.R
                if !isopt(r) continue end
                while true
                    nᵗ = d
                    c  = C[r.iˢ]
                    nʰ = isequal(r.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ]
                    removenode!(c, nᵗ, nʰ, r, s)
                    n += 1
                    if isequal(nʰ, d) break end
                end
            end
        end
        w[iᵈ] = 0
    end
    postremoval(s)
    # Step 2: Return solution
    return s
end

# Related Depot Removal
# Select a random closed depot node to open and remove q customer nodes most related to this depot node
function relateddepot!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    A = s.A
    preremoval(s)
    # Step 1: Select a random closed depot node
    iᵒ = sample(rng, eachindex(D), Weights(isclose.(D)))
    # Step 2: Evaluate relatedness of this depot node to every customer node
    x  = fill(-Inf, eachindex(C))
    for iⁿ ∈ eachindex(C) x[iⁿ] = relatedness(C[iⁿ], D[iᵒ], A[(iⁿ,iᵒ)]) end
    # Step 3: Remove at least q customer nodes most related to this pivot depot node
    n = 0
    while n < q 
        iⁿ = argmax(x)
        c  = C[iⁿ]
        r  = c.r
        nᵗ = isequal(r.iˢ, c.iⁿ) ? D[c.iᵗ] : C[c.iᵗ]
        nʰ = isequal(r.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ] 
        removenode!(c, nᵗ, nʰ, r, s)
        n += 1
        x[iⁿ] = -Inf
    end
    postremoval(s)
    # Step 4: Return solution
    return s
end

# Worst Depot Removal
# Iteratively select low-utilization depot and remove customer nodes from it until at least q customer nodes are removed
function worstdepot!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    preremoval(s)
    # Step 1: Evaluate utilization for each depot
    x = fill(Inf, eachindex(D))
    for (iᵈ,d) ∈ pairs(D)
        u = 0.
        for v ∈ d.V for r ∈ v.R u += r.q/d.q end end
        x[iᵈ] = u
    end
    # Step 2: Iteratively select low-utilization route and remove customer nodes from it until at least q customer nodes are removed
    n = 0
    w = ones(Int64, eachindex(D))
    while n < q
        if isone(sum(w)) break end
        iᵈ = argmin(x)
        d  = D[iᵈ]
        for v ∈ d.V
            if n ≥ q break end
            for r ∈ v.R
                if !isopt(r) continue end
                while true
                    nᵗ = d
                    c  = C[r.iˢ]
                    nʰ = isequal(r.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ]
                    removenode!(c, nᵗ, nʰ, r, s)
                    n += 1
                    if isequal(nʰ, d) break end
                end
            end
        end
        x[iᵈ] = Inf
        w[iᵈ] = 0
    end
    postremoval(s)
    # Step 3: Return solution
    return s
end