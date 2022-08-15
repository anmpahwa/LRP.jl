"""
    remove!([rng], q::Int64, s::Solution, method::Symbol)

Return solution removing q customer nodes from solution s using the given `method`.

Available methods include,
- Random Node Removal   : `randomnode!`
- Related Pair Removal  : `relatedpair!`
- Related Node Removal  : `relatednode!`
- Worst Node Removal    : `worstnode!`
- Random Route Removal  : `randomroute!`
- Worst Route Removal   : `worstroute!`
- Random Vehicle Removal: `randomvehicle!`
- Random Depot Removal  : `randomdepot!` 
- Related Depot Removal : `relateddepot!`

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
    # Step 1: Randomly select customer nodes to remove until q customer nodes have been removed
    n = 0
    while n < q
        c = rand(rng, C)
        if isopen(c) continue end
        r = c.r
        nₜ = isequal(r.s, c.i) ? D[c.t] : C[c.t]
        nₕ = isequal(r.e, c.i) ? D[c.h] : C[c.h]
        removenode!(c, nₜ, nₕ, r, s)
        n += 1
    end
    # Step 2: Return solution
    return s
end

# Related Node Pair Removal (related in pairs)
# Randomly remove q/2 customer nodes and the corresponding most related customer nodes
function relatedpair!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    A = s.A
    # Step 1: Randomly select initial q/2 customer nodes to remove
    n = ceil(Int64, q/2)
    randomnode!(rng, n, s)
    # Step 2: Compute relatedness for every closed customer node to every open customer node
    x = fill(-Inf, (eachindex(C), eachindex(C)))
    for (i,c₁) ∈ pairs(C)
        if isclose(c₁) continue end
        for (j,c₂) ∈ pairs(C)
            if isopen(c₂) continue end
            a = A[(i,j)]
            x[i,j] = relatedness(c₁, c₂, a)
        end
    end 
    # Step 3: For each q/2 customer node initially removed, remove the most related closed customer node  
    while n < q
        i,j = Tuple(argmax(x))
        c = C[j]
        r = c.r
        nₜ = isequal(r.s, c.i) ? D[c.t] : C[c.t]
        nₕ = isequal(r.e, c.i) ? D[c.h] : C[c.h]
        removenode!(c, nₜ, nₕ, r, s)
        n += 1
        x[i,:] .= -Inf
        x[:,j] .= -Inf
    end
    # Step 4: Return solution
    return s 
end

# Related Node Removal (related to pivot)
# For a randomly selected customer node, remove q most related customer nodes
function relatednode!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    A = s.A
    # Step 1: Randomly select a pivot customer node
    i = rand(rng, eachindex(C))
    # Step 2: For each customer node, evaluate relatedness to this pivot customer node
    x = fill(-Inf, eachindex(C))    # x[j]: relatedness of customer node C[i] with customer node C[j]  
    for j ∈ eachindex(C) x[j] = relatedness(C[i], C[j], A[(i,j)]) end
    # Step 3: Remove q most related customer nodes
    for _ ∈ 1:q
        k = argmax(x)
        c = C[k]
        r = c.r
        nₜ = isequal(r.s, c.i) ? D[c.t] : C[c.t]
        nₕ = isequal(r.e, c.i) ? D[c.h] : C[c.h]
        removenode!(c, nₜ, nₕ, r, s)
        x[k] = -Inf
    end
    # Step 4: Return solution
    return s
end

# Worst Node Removal
# Remove q customer nodes with highest removal cost
function worstnode!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    V = s.V
    x = fill(Inf, eachindex(C))     # x[k]: removal cost of customer node C[k]
    ϕ = ones(Int64, length(V))      # ϕ[k]: selection weight for vehicle V[k]
    # Step 1: Iterate until q customer nodes have been removed
    for _ ∈ 1:q
        # Step 1.1: For every closed customer node evaluate removal cost
        z = f(s)
        for (k,c) ∈ pairs(C)
            r = c.r
            if isopen(c) || iszero(ϕ[r.o]) continue end
            # Step 1.1.1: Remove closed customer node c between tail node nₜ and head node nₕ in route r
            nₜ = isequal(r.s, c.i) ? D[c.t] : C[c.t]
            nₕ = isequal(r.e, c.i) ? D[c.h] : C[c.h]
            removenode!(c, nₜ, nₕ, r, s)
            # Step 1.1.2: Evaluate the removal cost
            z⁻ = f(s)
            x[k] = z - z⁻
            # Step 1.1.3: Re-insert customer node c between tail node nₜ and head node nₕ in route r
            insertnode!(c, nₜ, nₕ, r, s)
        end
        # Step 1.2: Remove the customer node with highest removal cost
        k = argmax(x)
        c = C[k]
        r = c.r
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

# -------------------------------------------------- ROUTE REMOVAL --------------------------------------------------
# Random Route Removal
# Iteratively select a random route and remove customer nodes from it until at least q customer nodes are removed
function randomroute!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    V = s.V
    # Step 1: Iteratively select a random route and remove customer nodes from it until at least q customer nodes are removed
    n = 0
    while n < q
        v = rand(rng, V)
        r = rand(rng, v.R)
        if isclose(r) continue end
        v = V[r.o]
        d = D[v.o]
        while true
            nₜ = d
            c  = C[r.s]
            nₕ = isequal(r.e, c.i) ? D[c.h] : C[c.h]
            removenode!(c, nₜ, nₕ, r, s)
            n += 1
            if isequal(nₕ, d) break end
        end
    end
    # Step 2: Return solution
    return s
end

# Related Route Removal
# For a randomly selected route, remove customer nodes from most related route until q customer nodes are removed
function relatedroute!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    V = s.V
    R = [r for v ∈ V for r ∈ v.R]
    # Step 1: Randomly select a pivot route
    w = isopen.(R)
    i = sample(rng, eachindex(R), Weights(w))  
    # Step 2: For each route, evaluate relatedness to this pivot route
    x = fill(-Inf, length(R))
    for j ∈ eachindex(R) x[j] = iszero(w[j]) ? -Inf : relatedness(R[i], R[j]) end
    n = 0
    while n < q
        k = argmax(x)
        r = R[k]
        v = V[r.o]
        d = D[v.o]
        while true
            nₜ = d
            c  = C[r.s]
            nₕ = isequal(r.e, c.i) ? D[c.h] : C[c.h]
            removenode!(c, nₜ, nₕ, r, s)
            n += 1
            if isequal(nₕ, d) break end
        end 
        x[k] = -Inf
    end
    # Step 3: Return solution
    return s
end

# Worst Route Removal
# Iteratively select low-utilization route and remove customer nodes from it until at least q customer nodes are removed
function worstroute!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    V = s.V
    R = [r for v ∈ V for r ∈ v.R]
    # Step 1: Iteratively select low-utilization route and remove customer nodes from it until at least q customer nodes are removed
    n = 0
    x = fill(Inf, length(R))
    for (k,r) ∈ pairs(R) x[k] = isclose(r) ? Inf : r.l end
    while n < q
        k = argmin(x)
        r = R[k]
        v = V[r.o]
        d = D[v.o]
        while true
            nₜ = d
            c  = C[r.s]
            nₕ = isequal(r.e, c.i) ? D[c.h] : C[c.h]
            removenode!(c, nₜ, nₕ, r, s)
            n += 1
            if isequal(nₕ, d) break end
        end
        x[k] = Inf
    end
    # Step 2: Return solution
    return s
end
    
# -------------------------------------------------- VEHICLE REMOVAL --------------------------------------------------
# Random Vehicle Removal
# Iteratively select a random vehicle and remove customer nodes from it until at least q customer nodes are removed
function randomvehicle!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    V = s.V
    # Step 1: Iteratively select a random vehicle and remove customer nodes from it until at least q customer nodes are removed
    n = 0
    while n < q
        v = rand(rng, V)
        if isclose(v) continue end
        for r ∈ v.R
            if isclose(r) continue end
            d = D[v.o]
            while true
                nₜ = d
                c  = C[r.s]
                nₕ = isequal(r.e, c.i) ? D[c.h] : C[c.h]
                removenode!(c, nₜ, nₕ, r, s)
                n += 1
                if isequal(nₕ, d) break end
            end
        end
    end
    # Step 2: Return solution
    return s
end

# -------------------------------------------------- DEPOT REMOVAL --------------------------------------------------
# Random Depot Removal
# Iteratively select a random depot and remove customer nodes from it until at least q customer nodes are removed
function randomdepot!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    # Step 1: Iteratively select a random depot and remove customer nodes from it until at least q customer nodes are removed
    n = 0
    while n < q
        d = rand(rng, D)
        if isclose(d) continue end
        for v ∈ d.V
            for r ∈ v.R
                if isclose(r) continue end
                while true
                    nₜ = d
                    c  = C[r.s]
                    nₕ = isequal(r.e, c.i) ? D[c.h] : C[c.h]
                    removenode!(c, nₜ, nₕ, r, s)
                    n += 1
                    if isequal(nₕ, d) break end
                end
            end
        end
    end
    # Step 2: Return solution
    return s
end

# Related Depot Removal
# Select a random closed depot node to open and remove q customer nodes most related to this depot node
function relateddepot!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    A = s.A
    if all(isopen, D) return s end
    # Step 1: Select a random closed depot node
    w = isclose.(D)
    i = sample(rng, eachindex(D), Weights(w))
    # Step 2: Evaluate relatedness of this depot node to every customer node
    x = fill(-Inf, eachindex(C))
    for j ∈ eahchindex(C) x[j] = relatedness(D[i], C[j], A[(i,j)]) end
    # Step 3: Remove q customer nodes most related to this depot node
    for _ ∈ 1:q   
        k = argmax(x)
        c = C[k]
        r = c.r
        nₜ = isequal(r.s, c.i) ? D[c.t] : C[c.t]
        nₕ = isequal(r.e, c.i) ? D[c.h] : C[c.h]
        removenode!(c, nₜ, nₕ, r, s)
        x[k] = -Inf
    end
    # Step 4: Return solution
    return s
end

# -------------------------------------------------------------------------------------------------------------------
function relatedness(d₁::DepotNode, d₂::DepotNode, a::Arc) 
    l = a.l
    t = a.t
    f = a.f
    q = abs(d₁.q - d₂.q)
    ϕ = false
    z = 1/(0.482 * l + 35.0 * t + 3.826 * f - 1.0 * q - 1.0 * ϕ)
    return z
end
function relatedness(c₁::CustomerNode, c₂::CustomerNode, a::Arc)
    l = a.l
    t = a.t
    f = a.f
    q = abs(c₁.q - c₂.q)
    ϕ = isequal(c₁.r, c₂.r)
    r = 1/(0.482 * l + 35.0 * t + 3.826 * f - 1.0 * q - 1.0 * ϕ)
    return r
end
function relatedness(c::CustomerNode, d::DepotNode, a::Arc)
    l = a.l
    t = a.t
    f = a.f
    q = abs(c.q - d.q)
    ϕ = false
    r = c.r
    for v ∈ d.V 
        ϕ = isequal(r.o, v.i) 
        ϕ ? break : continue
    end
    z = 1/(0.482 * l + 35.0 * t + 3.826 * f - 1.0 * q - 1.0 * ϕ)
    return z
end
function relatedness(d::DepotNode, c::CustomerNode, a::Arc)
    l = a.l
    t = a.t
    f = a.f
    q = abs(c.q - d.q)
    ϕ = false
    r = c.r
    for v ∈ d.V 
        ϕ = isequal(r.o, v.i) 
        ϕ ? break : continue
    end
    z = 1/(0.482 * l + 35.0 * t + 3.826 * f - 1.0 * q - 1.0 * ϕ)
    return z
end
function relatedness(r₁::Route, r₂::Route)
    l = abs(r₁.l - r₂.l)
    t = abs(r₁.t - r₂.t)
    f = abs(r₁.f - r₂.f)
    q = abs(r₁.q - r₂.q)
    ϕ = isequal(r₁.o, r₂.o)
    z = 1/(0.482 * l + 35.0 * t + 3.826 * f - 1.0 * q - 1.0 * ϕ)
    return z
end