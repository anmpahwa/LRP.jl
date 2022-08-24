# Initial solution
"""
    initialsolution([rng], instance, method)

Returns initial LRP solution for the given `instance` using the given `method`.

Available methods include,
- Clarke and Wright Savings Algorithm   : `:cw`
- Nearest Neighborhood Algorithm        : `:nn`
- Random Initialization                 : `:random`
- Regret N Insertion                    : `:regret₂init`, `:regret₃init`

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
initialsolution(rng::AbstractRNG, instance, method::Symbol)::Solution = getfield(LRP, method)(rng, instance)
initialsolution(instance, method::Symbol) = initialsolution(Random.GLOBAL_RNG, instance, method)

# Clarke and Wright Savings Algorithm
# Create initial solution merging routes that render the most savings until no merger can render further savings
function cw(rng::AbstractRNG, instance)
    G = build(instance)
    s = Solution(G...)
    D = s.D
    C = s.C
    V = s.V
    R = s.R
    # Step 1: Initialize by assigning customer node to the vehicle-depot pair that results in least assignment cost
    I = eachindex(V)
    J = eachindex(C)
    x = fill(Inf, (I,J))                # x[i,j]: Cost of assigning customer node C[j] to vehicle V[i]
    # Step 2: Iterate through each customer node
    for (j,c) ∈ pairs(C)
        # Step 2.1: For customer node c, iterate through every vehicle-depot pair evaluating assignment cost
        zᵒ = f(s)
        for (i,v) ∈ pairs(V)
            # Step 2.1.1: Assign customer node c to vehicle v of depot node d
            d = D[v.o] 
            r = Route(j, v, d)
            push!(v.R, r)
            insertnode!(c, d, d, r, s)
            # Step 2.1.2: Compute assignment cost
            z⁺ = f(s)
            Δ  = (z⁺ + (d.πᶠ/d.q + v.πᶠ/v.q)) - (zᵒ - (d.πᶠ + v.πᶠ))
            x[i,j] = Δ
            # Step 2.1.3: Unassign customer node c from vehicle v of depot node d
            removenode!(c, d, d, r, s)
            pop!(v.R)
        end
        # Step 2.2: Assign customer node c to vehicle v of depot node v that results in least assignment cost
        i = argmin(x[:,j])
        v = V[i]
        d = D[v.o]
        r = Route(rand(rng, 1:M), v, d)
        push!(v.R, r)
        insertnode!(c, d, d, r, s)
        # Step 2.3: Revise vectors appropriately
        x[:,j] .= Inf
    end
    for v ∈ V append!(R, v.R) end
    # Step 3: Merge routes iteratively until no merger can render further savings
    K = eachindex(R)
    y = fill(-Inf, (K,K))               # y[i,j]: Savings from merging route R[i] into R[j] 
    ϕ = ones(Int64, K)                  # ϕ[k]  : selection weight for route R[k]  
    while true
        # Step 3.1: Iterate through every route-pair combination
        zᵒ = f(s)
        for (i,r₁) ∈ pairs(R)
            if !isopt(r₁) continue end
            for (j,r₂) ∈ pairs(R)
                #Δ = 0
                # Step 3.1.1: Merge route r₁ into route r₂
                if !isopt(r₂) continue end
                if isequal(r₁, r₂) continue end
                if iszero(ϕ[i]) & iszero(ϕ[j]) continue end
                v₁, v₂ = V[r₁.o], V[r₂.o]
                d₁, d₂ = D[v₁.o], D[v₂.o]
                cₛ, cₑ = C[r₁.s], C[r₁.e]
                while true
                    c  = C[r₁.s]
                    nₜ = d₁
                    nₕ = isequal(r₁.e, c.i) ? D[c.h] : C[c.h]
                    removenode!(c, nₜ, nₕ, r₁, s)
                    nₜ = C[r₂.e]
                    nₕ = d₂
                    insertnode!(c, nₜ, nₕ, r₂, s)
                    if isequal(c, cₑ) break end
                end
                # Step 3.1.2: Compute savings from merging route r₁ into route r₂
                z⁻ = f(s)
                Δ  = z⁻ - zᵒ
                y[i,j] = -Δ
                # Step 3.1.3: Unmerge routes r₁ and r₂
                while true
                    c  = C[r₂.e] 
                    nₜ = C[c.t]
                    nₕ = d₂
                    removenode!(c, nₜ, nₕ, r₂, s)
                    nₜ = d₁
                    nₕ = isopt(r₁) ? C[r₁.s] : D[r₁.s]
                    insertnode!(c, nₜ, nₕ, r₁, s)
                    if isequal(c, cₛ) break end
                end
            end
        end
        # Step 3.2: Merge routes that render highest savings. If no route render savings, go to step 4. 
        if maximum(y) < 0 break end
        i,j = Tuple(argmax(y))
        r₁, r₂ = R[i], R[j]
        v₁, v₂ = V[r₁.o], V[r₂.o]
        d₁, d₂ = D[v₁.o], D[v₂.o]
        cₛ, cₑ = C[r₁.s], C[r₁.e]
        while true
            c  = C[r₁.s]
            nₜ = d₁
            nₕ = isequal(r₁.e, c.i) ? D[c.h] : C[c.h]
            removenode!(c, nₜ, nₕ, r₁, s)
            nₜ = C[r₂.e]
            nₕ = d₂
            insertnode!(c, nₜ, nₕ, r₂, s)
            if isequal(c, cₑ) break end
        end
        # Step 3.3: Revise savings and selection vectors appropriately
        y[i,:] .= -Inf
        y[:,i] .= -Inf
        for (j,r) ∈ pairs(R) 
            if isequal(r.o, v₂.i) y[j,:] .= -Inf end 
            if isequal(r.o, v₂.i) y[:,j] .= -Inf end
        end
        for (k,r) ∈ pairs(R) ϕ[k] = isequal(r.o, r₂.o) ? 1 : 0 end
    end
    # Step 5: Return initial solution
    deleteat!(R, findall(!isopt, R))
    for (k,v) ∈ pairs(V) deleteat!(v.R, findall(!isopt, v.R)) end
    return s
end

# Nearest Neighborhood Algorithm
# Create initial solution with node-route combination that results in least increase in cost until all customer nodes have been added to the solution
function nn(rng::AbstractRNG, instance)
    G = build(instance)
    s = Solution(G...)
    D = s.D
    C = s.C
    V = s.V
    R = s.R
    # Step 1: Initialize an empty route for every vehicle
    for d ∈ D for v ∈ d.V push!(v.R, Route(rand(rng, 1:M), v, d)) end end
    for v ∈ V append!(R, v.R) end
    I = eachindex(C)
    J = eachindex(R)
    x = ElasticMatrix(fill(Inf, (I,J)))     # x[i,j]: insertion cost of customer node C[i] in route R[j]
    ϕ = ones(Int64, J)                      # ϕ[j]  : selection weight for route R[j]
    # Step 2: Iterate until all customer nodes have been added to the routes
    for _ ∈ I
        # Step 2.1: Iteratively compute cost of appending each open customer node in each route
        zᵒ = f(s)
        for (i,c) ∈ pairs(C)
            if !isopen(c) continue end
            for (j,r) ∈ pairs(R)   
                # Step 2.2.1: Append customer node c in the route r
                if iszero(ϕ[j]) continue end
                v  = V[r.o]
                d  = D[v.o]
                nₜ = isopt(r) ? C[r.e] : D[r.e]
                nₕ = d
                insertnode!(c, nₜ, nₕ, r, s)
                # Step 2.2.2: Compute increase in cost
                z⁺ = f(s)
                Δ  = (z⁺ + (d.πᶠ/d.q + v.πᶠ/v.q)) - (zᵒ - (d.πᶠ + v.πᶠ))
                x[i,j] = Δ
                # Step 2.2.3: Pop customer node c from the route r
                removenode!(c, nₜ, nₕ, r, s)
            end
        end
        # Step 2.2: Choose node-route combination that results in least increase in cost and append this customer node into this route
        i,j = Tuple(argmin(x))
        c = C[i]
        r = R[j]
        v = V[r.o]
        d = D[v.o]
        nₜ = isopt(r) ? C[r.e] : D[r.e]
        nₕ = d
        insertnode!(c, nₜ, nₕ, r, s)
        # Step 2.3: Revise cost and selection vectors appropriately
        ϕ .= 0
        x[i,:] .= Inf
        for (j,r) ∈ pairs(R)
            if !isequal(r.o, v.i) continue end
            x[:,j] .= Inf
            ϕ[j] = 1
        end 
        for (k,v) ∈ pairs(V)
            if any(!isopt, v.R) continue end
            r = Route(rand(rng, 1:M), V[k], D[v.o])
            push!(v.R, r) 
            push!(R, r)
            append!(x, fill(Inf, (I,1)))
            push!(ϕ, 1)
        end
    end
    # Step 3: Return initial solution
    deleteat!(R, findall(!isopt, R))
    for (k,v) ∈ pairs(V) deleteat!(v.R, findall(!isopt, v.R)) end
    return s
end

# Random Initialization
# Create initial solution with randomly selcted node-route combination until all customer nodes have been added to the solution
function random(rng::AbstractRNG, instance)
    G = build(instance)
    s = Solution(G...)
    D = s.D
    C = s.C
    V = s.V
    R = s.R
    # Step 1: Initialize an empty route for every vehicle
    for d ∈ D for v ∈ d.V push!(v.R, Route(rand(rng, 1:M), v, d)) end end
    for v ∈ V append!(R, v.R) end
    I = eachindex(C)
    J = eachindex(R)
    w = ones(Int64, I)                      # w[i]: selection weight for customer node C[i]
    # Step 2: Iteratively append randomly selected customer node in randomly selected route
    for _ ∈ I
        i = sample(rng, I, OffsetWeights(w))
        j = sample(rng, J)
        c = C[i]
        r = R[j]
        v = V[r.o]
        d = D[v.o]
        nₜ = d
        nₕ = isopt(r) ? C[r.s] : D[r.s]
        insertnode!(c, nₜ, nₕ, r, s)
        w[i] = 0
    end
    # Step 3: Return initial solution
    deleteat!(R, findall(!isopt, R))
    for (k,v) ∈ pairs(V) deleteat!(v.R, findall(!isopt, v.R)) end
    return s
end

# Regret-N Insertion
# Create initial solution by iteratively adding customer nodes with highest regret cost at its best position until all customer nodes have been added to the solution
function regretₙinit(rng::AbstractRNG, N::Int64, instance)
    G = build(instance)
    s = Solution(G...)
    D = s.D
    C = s.C
    V = s.V
    R = s.R
    # Step 1: Initialize an empty route for every vehicle
    for d ∈ D for v ∈ d.V push!(v.R, Route(rand(rng, 1:M), v, d)) end end
    for v ∈ V append!(R, v.R) end
    I = eachindex(C)
    J = eachindex(R)
    x = ElasticMatrix(fill(Inf, (I,J)))     # x[i,j]: insertion cost of customer node C[i] at best position in route R[j]
    p = ElasticMatrix(fill((0, 0), (I,J)))  # p[i,j]: best insertion postion of customer node C[i] in route R[j]
    y = fill(Inf, (I,N))                    # y[i,n]: insertion cost of customer node C[i] at nᵗʰ best position
    w = zeros(Int64, (I,N))                 # w[i,n]: route index of customer node C[j] at nᵗʰ best position
    z = fill(-Inf, I)                       # z[i]  : regret-N cost of customer node C[i]
    ϕ = ones(Int64, J)                      # ϕ[j]  : selection weight for route R[j]
    # Step 2: Iterate until all customer nodes have been inserted into the route
    for _ ∈ I
        # Step 2.1: Iterate through all open customer nodes and every route
        zᵒ = f(s)
        for (i,c) ∈ pairs(C)
            if !isopen(c) continue end
            for (j,r) ∈ pairs(R)
                # Step 2.1.1: Iterate through all possible insertion position in route r
                if iszero(ϕ[j]) continue end
                v  = V[r.o]
                d  = D[v.o]
                nₛ = isopt(r) ? C[r.s] : D[r.s]
                nₑ = isopt(r) ? C[r.e] : D[r.e]
                nₜ = d
                nₕ = nₛ
                while true
                    # Step 2.1.1.1: Insert customer node c between tail node nₜ and head node nₕ in route r, and compute the insertion cost
                    insertnode!(c, nₜ, nₕ, r, s)
                    # Step 2.1.1.2: Compute the insertion cost
                    z⁺ = f(s)
                    Δ  = (z⁺ + (d.πᶠ/d.q + v.πᶠ/v.q)) - (zᵒ - (d.πᶠ + v.πᶠ))
                    # Step 2.1.1.3: Revise least insertion cost in route r and the corresponding best insertion position in route r
                    if Δ < x[i,j] x[i,j], p[i,j] = Δ, (nₜ.i, nₕ.i) end
                    # Step 2.1.1.4: Revise N least insertion costs
                    n̲ = 1
                    for n ∈ 1:N 
                        n̲ = n
                        if Δ < y[i,n] break end
                    end
                    for n ∈ N:-1:n̲ 
                        y[i,n] = isequal(n, n̲) ? Δ : y[i,n-1]::Float64
                        w[i,n] = isequal(n, n̲) ? r.i : w[i,n-1]::Int64
                    end
                    # Step 2.1.1.5: Remove customer node c from its position between tail node nₜ and head node nₕ in route r
                    removenode!(c, nₜ, nₕ, r, s)
                    if isequal(nₜ, nₑ) break end
                    nₜ = nₕ
                    nₕ = isequal(r.e, nₜ.i) ? D[nₜ.h] : C[nₜ.h]
                end
            end
            # Step 2.1.2: Compute regret cost for customer node c
            z[i] = 0.
            for n ∈ 1:N z[i] += y[i,n] - y[i,1] end
        end
        # Step 2.2: Insert customer node with highest regret cost in its best position (break ties by inserting the node with the lowest insertion cost)
        I̲ = findall(i -> i == maximum(z), z)
        i,j = Tuple(argmin(x[I̲,:]))
        i = I̲[i]
        c = C[i]
        r = R[j]
        v = V[r.o]
        t, h = p[i,j]
        nₜ = t ≤ length(D) ? D[t] : C[t]
        nₕ = h ≤ length(D) ? D[h] : C[h]
        insertnode!(c, nₜ, nₕ, r, s)
        # Step 2.3: Revise vectors appropriately
        x[i,:] .= Inf
        p[i,:] .= ((0, 0), )
        y[i,:] .= Inf
        w[i,:] .= 0
        z .= -Inf 
        for (i,c) ∈ pairs(C)
            for n ∈ 1:N
                if iszero(w[i,n]) break end
                k = findfirst(x -> x.i == w[i,n], R)
                r = R[k]
                if isequal(r.o, v.i) y[i,n], w[i,n] = Inf, 0 end
            end
            ix = sortperm(y[i,:])
            y[i,:] .= y[i,ix]
            w[i,:] .= w[i,ix]
        end
        ϕ .= 0
        for (j,r) ∈ pairs(R)
            if !isequal(r.o, v.i) continue end
            x[:,j] .= Inf
            ϕ[j] = 1
        end
        for (k,v) ∈ pairs(V)
            if any(!isopt, v.R) continue end
            r = Route(rand(rng, 1:M), V[k], D[v.o])
            push!(v.R, r) 
            push!(R, r)
            append!(x, fill(Inf, (I,1)))
            append!(p, fill((0, 0), (I,1)))
            push!(ϕ, 1)
        end
    end
    # Step 3: Return initial solution
    deleteat!(R, findall(!isopt, R))
    for (k,v) ∈ pairs(V) deleteat!(v.R, findall(!isopt, v.R)) end
    return s
end
regret₂init(rng::AbstractRNG, instance) = regretₙinit(rng, Int64(2), instance)
regret₃init(rng::AbstractRNG, instance) = regretₙinit(rng, Int64(3), instance)
