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
    V = s.V
    R = s.R
    # Step 1: Randomly select customer nodes to remove until q customer nodes have been removed
    n = 0
    w = isclose.(C)
    while n < q
        k = sample(rng, eachindex(C), OffsetWeights(w))
        c = C[k]
        if isopen(c) continue end
        r = c.r
        nₜ = isequal(r.s, c.i) ? D[c.t] : C[c.t]
        nₕ = isequal(r.e, c.i) ? D[c.h] : C[c.h]
        removenode!(c, nₜ, nₕ, r, s)
        n += 1
        w[k] = 0
    end
    # Step 2: Return solution
    deleteat!(R, findall(!isopt, R))
    for (k,v) ∈ pairs(V) deleteat!(v.R, findall(!isopt, v.R)) end
    return s
end

# Related Node Removal (related to pivot)
# For a randomly selected customer node, remove q most related customer nodes
function relatednode!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    A = s.A
    V = s.V
    R = s.R
    # Step 1: Randomly select a pivot customer node
    i = rand(rng, eachindex(C))
    # Step 2: For each customer node, evaluate relatedness to this pivot customer node
    x = fill(-Inf, eachindex(C))    # x[j]: relatedness of customer node C[i] with customer node C[j]  
    for j ∈ eachindex(C) x[j] = relatedness(C[i], C[j], A[(i,j)]) end
    # Step 3: Remove q most related customer nodes
    n = 0
    while n < q
        k = argmax(x)
        c = C[k]
        r = c.r
        nₜ = isequal(r.s, c.i) ? D[c.t] : C[c.t]
        nₕ = isequal(r.e, c.i) ? D[c.h] : C[c.h]
        removenode!(c, nₜ, nₕ, r, s)
        n += 1
        x[k] = -Inf
    end
    # Step 4: Return solution
    deleteat!(R, findall(!isopt, R))
    for (k,v) ∈ pairs(V) deleteat!(v.R, findall(!isopt, v.R)) end
    return s
end

# Worst Node Removal
# Remove q customer nodes with highest removal cost (savings)
function worstnode!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    V = s.V
    R = s.R
    x = fill(-Inf, eachindex(C))    # x[i]: removal cost of customer node C[i]
    ϕ = ones(Int64, eachindex(V))   # ϕ[j]: selection weight for route V[j]
    # Step 1: Iterate until q customer nodes have been removed
    n = 0
    while n < q
        # Step 1.1: For every closed customer node evaluate removal cost
        zᵒ = f(s)
        for (i,c) ∈ pairs(C)
            if isopen(c) continue end
            r = c.r
            j = r.o
            if iszero(ϕ[j]) continue end
            # Step 1.1.1: Remove closed customer node c between tail node nₜ and head node nₕ in route r
            nₜ = isequal(r.s, c.i) ? D[c.t] : C[c.t]
            nₕ = isequal(r.e, c.i) ? D[c.h] : C[c.h]
            removenode!(c, nₜ, nₕ, r, s)
            # Step 1.1.2: Evaluate the removal cost
            z⁻ = f(s)
            Δ  = z⁻ - zᵒ
            x[i] = -Δ
            # Step 1.1.3: Re-insert customer node c between tail node nₜ and head node nₕ in route r
            insertnode!(c, nₜ, nₕ, r, s)
        end
        # Step 1.2: Remove the customer node with highest removal cost (savings)
        i = argmax(x)
        c = C[i]
        r = c.r
        v = V[r.o]
        nₜ = isequal(r.s, c.i) ? D[c.t] : C[c.t]
        nₕ = isequal(r.e, c.i) ? D[c.h] : C[c.h]
        removenode!(c, nₜ, nₕ, r, s)
        n += 1
        # Step 1.3: Update cost and selection weight vectors
        x[i] = -Inf
        for (j,v) ∈ pairs(V) ϕ[j] = isequal(r.o, v.i) ? 1 : 0 end
        if isopt(r) continue end
    end
    # Step 2: Return solution
    deleteat!(R, findall(!isopt, R))
    for (k,v) ∈ pairs(V) deleteat!(v.R, findall(!isopt, v.R)) end
    return s
end

# -------------------------------------------------- ROUTE REMOVAL --------------------------------------------------
# Random Route Removal
# Iteratively select a random route and remove customer nodes from it until at least q customer nodes are removed
function randomroute!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    V = s.V
    R = s.R
    # Step 1: Iteratively select a random route and remove customer nodes from it until at least q customer nodes are removed
    n = 0
    w = isopt.(R)
    while n < q
        if isone(sum(w)) break end
        k = sample(rng, eachindex(R), Weights(w))
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
        w[k] = 0
    end
    # Step 2: Return solution
    deleteat!(R, findall(!isopt, R))
    for (k,v) ∈ pairs(V) deleteat!(v.R, findall(!isopt, v.R)) end
    return s
end

# Related Route Removal
# For a randomly selected route, remove customer nodes from most related route until q customer nodes are removed
function relatedroute!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    V = s.V
    R = s.R
    # Step 1: Randomly select a pivot route
    i = sample(rng, eachindex(R), Weights(isopt.(R)))  
    # Step 2: For each route, evaluate relatedness to this pivot route
    x = fill(-Inf, eachindex(R))
    for (j,r) ∈ pairs(R) x[j] = !isopt(r) ? -Inf : relatedness(R[i], R[j]) end
    # Step 3: Remove at least q customers from most related route to this pivot route
    n = 0
    w = isopt.(R)
    while n < q
        if isone(sum(w)) break end
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
        w[k] = 0
    end
    # Step 4: Return solution
    deleteat!(R, findall(!isopt, R))
    for (k,v) ∈ pairs(V) deleteat!(v.R, findall(!isopt, v.R)) end
    return s
end

# Worst Route Removal
# Iteratively select low-utilization route and remove customer nodes from it until at least q customer nodes are removed
function worstroute!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    V = s.V
    R = s.R
    # Step 1: Evaluate utilization of each route
    x = fill(Inf, eachindex(R))
    for (k,r) ∈ pairs(R)
        if !isopt(r) continue end
        v = V[r.o]
        u = r.l/v.q
        x[k] = u
    end
    # Step 2: Iteratively select low-utilization route and remove customer nodes from it until at least q customer nodes are removed
    n = 0
    w = isopt.(R)
    while n < q
        if isone(sum(w)) break end
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
        w[k] = 0
    end
    # Step 3: Return solution
    deleteat!(R, findall(!isopt, R))
    for (k,v) ∈ pairs(V) deleteat!(v.R, findall(!isopt, v.R)) end
    return s
end
    
# -------------------------------------------------- VEHICLE REMOVAL --------------------------------------------------
# Random Vehicle Removal
# Iteratively select a random vehicle and remove customer nodes from it until at least q customer nodes are removed
function randomvehicle!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    V = s.V
    R = s.R
    # Step 1: Iteratively select a random vehicle and remove customer nodes from it until at least q customer nodes are removed
    n = 0
    w = isopt.(V)
    while n < q
        if isone(sum(w)) break end
        k = sample(rng, eachindex(V), Weights(w))
        v = V[k]
        d = D[v.o]
        for r ∈ v.R
            if !isopt(r) continue end
            while true
                nₜ = d
                c  = C[r.s]
                nₕ = isequal(r.e, c.i) ? D[c.h] : C[c.h]
                removenode!(c, nₜ, nₕ, r, s)
                n += 1
                if isequal(nₕ, d) break end
            end
        end
        w[k] = 0
    end
    # Step 2: Return solution
    deleteat!(R, findall(!isopt, R))
    for (k,v) ∈ pairs(V) deleteat!(v.R, findall(!isopt, v.R)) end
    return s
end

function relatedvehicle!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    V = s.V
    R = s.R
    # Step 1: Select a random closed depot node
    i = sample(rng, eachindex(V), Weights(isopt.(V)))
    # Step 2: For each vehicle, evaluate relatedness to this pivot vehicle
    x = fill(-Inf, eachindex(V))
    for (j,v) ∈ pairs(V) x[j] = !isopt(v) ? -Inf : relatedness(V[i], V[j]) end
    # Step 3: Remove at least q customers from the most related vehicles to this pivot vehicle
    n = 0
    w = isopt.(V)
    while n < q
        if isone(sum(w)) break end
        k = argmax(x)
        v = V[k]
        d = D[v.o]
        for r ∈ v.R
            if !isopt(r) continue end
            while true
                nₜ = d
                c  = C[r.s]
                nₕ = isequal(r.e, c.i) ? D[c.h] : C[c.h]
                removenode!(c, nₜ, nₕ, r, s)
                n += 1
                if isequal(nₕ, d) break end
            end
        end
        x[k] = -Inf
        w[k] = 0
    end
    # Step 4: Return solution
    deleteat!(R, findall(!isopt, R))
    for (k,v) ∈ pairs(V) deleteat!(v.R, findall(!isopt, v.R)) end
    return s
end

# Worst Vehicle Removal
# Iteratively select low-utilization vehicle and remove customer nodes from it until at least q customer nodes are removed
function worstvehicle!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    V = s.V
    R = s.R
    # Step 1: Evaluate utilization for each vehicle
    x = fill(Inf, eachindex(V))
    for v ∈ V
        if !isopt(v) continue end
        a = 0
        b = 0
        k = v.i
        for r ∈ v.R
            if !isopt(r) continue end
            a += r.l
            b += v.q
        end
        u = a/b
        x[k] = u
    end
    # Step 2: Iteratively select low-utilization route and remove customer nodes from it until at least q customer nodes are removed
    n = 0
    w = isopt.(V)
    while n < q
        if isone(sum(w)) break end
        k = argmin(x)
        v = V[k]
        d = D[v.o]
        for r ∈ v.R
            if !isopt(r) continue end
            while true
                nₜ = d
                c  = C[r.s]
                nₕ = isequal(r.e, c.i) ? D[c.h] : C[c.h]
                removenode!(c, nₜ, nₕ, r, s)
                n += 1
                if isequal(nₕ, d) break end
            end
        end
        x[k] = Inf
        w[k] = 0
    end
    # Step 3: Return solution
    deleteat!(R, findall(!isopt, R))
    for (k,v) ∈ pairs(V) deleteat!(v.R, findall(!isopt, v.R)) end
    return s
end

# -------------------------------------------------- DEPOT REMOVAL --------------------------------------------------
# Random Depot Removal
# Iteratively select a random depot and remove customer nodes from it until at least q customer nodes are removed
function randomdepot!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    V = s.V
    R = s.R
    # Step 1: Iteratively select a random depot and remove customer nodes from it until at least q customer nodes are removed
    n = 0
    w = isopt.(D)
    while n < q
        if isone(sum(w)) break end
        k = sample(rng, eachindex(D), Weights(w))
        d = D[k]
        for v ∈ d.V
            if !isopt(v) continue end
            for r ∈ v.R
                if !isopt(r) continue end
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
        w[k] = 0
    end
    # Step 2: Return solution
    deleteat!(R, findall(!isopt, R))
    for (k,v) ∈ pairs(V) deleteat!(v.R, findall(!isopt, v.R)) end
    return s
end

# Related Depot Removal
# Select a random closed depot node to open and remove q customer nodes most related to this depot node
function relateddepot!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    A = s.A
    V = s.V
    R = s.R
    # Step 1: Select a random closed depot node
    i = sample(rng, eachindex(D), Weights(isclose.(D)))
    # Step 2: Evaluate relatedness of this depot node to every customer node
    x = fill(-Inf, eachindex(C))
    for j ∈ eachindex(C) x[j] = relatedness(D[i], C[j], A[(i,j)]) end
    # Step 3: Remove at least q customer nodes most related to this pivot depot node
    n = 0
    while n < q 
        k = argmax(x)
        c = C[k]
        r = c.r
        v = V[r.o]
        d = D[v.o]
        nₜ = isequal(r.s, c.i) ? D[c.t] : C[c.t]
        nₕ = isequal(r.e, c.i) ? D[c.h] : C[c.h]
        removenode!(c, nₜ, nₕ, r, s)
        n += 1
        x[k] = -Inf
    end
    # Step 4: Return solution
    deleteat!(R, findall(!isopt, R))
    for (k,v) ∈ pairs(V) deleteat!(v.R, findall(!isopt, v.R)) end
    return s
end

# Worst Depot Removal
# Iteratively select low-utilization depot and remove customer nodes from it until at least q customer nodes are removed
function worstdepot!(rng::AbstractRNG, q::Int64, s::Solution)
    D = s.D
    C = s.C
    V = s.V
    R = s.R
    # Step 1: Evaluate utilization for each depot
    x = fill(Inf, eachindex(D))
    for d ∈ D
        if !isopt(d) continue end
        u = 0.
        k = d.i
        for v ∈ d.V
            if !isopt(v) continue end
            for r ∈ v.R
                if !isopt(r) continue end
                u += r.l/d.q
            end
        end
        x[k] = u
    end
    # Step 2: Iteratively select low-utilization route and remove customer nodes from it until at least q customer nodes are removed
    n = 0
    w = isopt.(D)
    while n < q
        if isone(sum(w)) break end
        k = argmin(x)
        d = D[k]
        for v ∈ d.V
            for r ∈ v.R
                if !isopt(r) continue end
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
        x[k] = Inf
        w[k] = 0
    end
    # Step 3: Return solution
    deleteat!(R, findall(!isopt, R))
    for (k,v) ∈ pairs(V) deleteat!(v.R, findall(!isopt, v.R)) end
    return s
end