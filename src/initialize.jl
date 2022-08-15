# Initial solution
"""
    initialsolution([rng], (D::Vector{DepotNode}, C:Vector{CustomerNode}, A::Dict{Tuple{Int64, Int64}, Arc}), method)

Returns initial LRP solution using given `method` on graph with depot nodes `D`, customer nodes `C`, and arcs `A`.

Available methods include,
- Clarke and Wright Savings Algorithm   : `:cw`
- Nearest Neighborhood Algorithm        : `:nn`
- Random Initialization                 : `:random`
- Regret N Insertion                    : `:regret₂init`, `:regret₃init`

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
initialsolution(rng::AbstractRNG, G, method::Symbol)::Solution = getfield(LRP, method)(rng, G)
initialsolution(G, method::Symbol) = initialsolution(Random.GLOBAL_RNG, G, method)

# Clarke and Wright Savings Algorithm
# Create initial solution merging routes that render the most savings until no merger can render further savings
function cw(rng::AbstractRNG, G)
    s = Solution(G...)
    D = s.D
    C = s.C
    V = s.V
    # Step 1: Initialize by assigning customer node to the vehicle-depot pair that results in least assignment cost
    M = typemax(Int64)
    I = length(V)
    J = eachindex(C)
    x = fill(Inf, (I,J))                # x[i,j]: Cost of assigning customer node C[j] to vehicle V[i] of depot node d
    # Step 2: Iterate through each customer node
    for (j,c) ∈ pairs(C)
        # Step 2.1: For customer node c, iterate through every vehicle-depot pair evaluating assignment cost
        for (i,v) ∈ pairs(V)
            # Step 2.1.1: Assign customer node c to vehicle v of depot node d
            d = D[v.o]
            r = Route(j, v, d)
            push!(v.R, r)
            insertnode!(c, d, d, r, s)
            # Step 2.1.2: Compute assignment cost
            x[i,j] = f(s)
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
        # Step 2.3: Close customer node c
        x[:,j] .= Inf
    end
    R = [r for v ∈ V for r ∈ v.R]
    # Step 3: Merge routes iteratively until no merger can render further savings
    K = length(R)
    y = fill(-Inf, (K,K))               # y[i,j]: Savings from merging route R[i] into R[j] 
    ϕ = ones(Int64, K)                  # ϕ[k]  : selection weight for route R[k]  
    while true
        # Step 3.1: Iterate through every route-pair combination
        z = f(s)
        for (i,r₁) ∈ pairs(R)
            if isclose(r₁) continue end
            for (j,r₂) ∈ pairs(R)
                # Step 3.1.1: Merge route r₁ into route r₂
                if isclose(r₂) continue end
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
                Δ  = z - z⁻
                y[i,j] = Δ
                # Step 3.1.3: Unmerge routes r₁ and r₂
                while true
                    c  = C[r₂.e] 
                    nₜ = C[c.t]
                    nₕ = d₂
                    removenode!(c, nₜ, nₕ, r₂, s)
                    nₜ = d₁
                    nₕ = isclose(r₁) ? d₁ : C[r₁.s]
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
        for (j,_) ∈ pairs(v₂.R) y[j,:] .= -Inf end
        for (j,_) ∈ pairs(v₂.R) y[:,j] .= -Inf end
        for (k,r) ∈ pairs(R) ϕ[k] = isequal(r.o, r₂.o) ? 1 : 0 end
    end
    # Step 4: Return initial solution
    for d ∈ D for v ∈ d.V filter!(isopen, v.R) end end
    return s
end

# Nearest Neighborhood Algorithm
# Create initial solution with node-route combination that results in least increase in cost until all customer nodes have been added to the solution
function nn(rng::AbstractRNG, G)
    s = Solution(G...)
    D = s.D
    C = s.C
    V = s.V
    # Step 1: Initialize an empty route for every vehicle
    M = typemax(Int64)
    for d ∈ D for v ∈ d.V push!(v.R, Route(rand(rng, 1:M), v, d)) end end
    R = [r for v ∈ V for r ∈ v.R]
    I = length(R)
    J = eachindex(C)
    K = length(V)
    x = fill(Inf, (I,J))                # x[i,j]: insertion cost of customer node C[j] in route R[i]
    ϕ = ones(Int64, K)                  # ϕ[k]  : selection weight for vehicle V[k]
    # Step 2: Iterate until all customer nodes have been added to the routes
    for _ ∈ eachindex(C)
        # Step 2.1: Iteratively compute cost of appending each open customer node in each route
        z = f(s)
        for (j,c) ∈ pairs(C)
            if isclose(c) continue end
            for (i,r) ∈ pairs(R)
                # Step 2.2.1: Append customer node c in the route r
                if iszero(ϕ[r.o]) continue end
                v = V[r.o]
                d = D[v.o]
                nₜ = isclose(r) ? D[r.e] : C[r.e]
                nₕ = d
                insertnode!(c, nₜ, nₕ, r, s)
                # Step 2.2.2: Compute increase in cost
                z⁺ = f(s)
                Δ  = z⁺ - z
                x[i,j] = Δ
                # Step 2.2.3: Pop customer node c from the route r
                removenode!(c, nₜ, nₕ, r, s)
            end
        end
        # Step 2.2: Choose node-route combination that results in least increase in cost and append this customer node into this route
        i,j = Tuple(argmin(x))
        r = R[i]
        c = C[j]
        v = V[r.o]
        d = D[v.o]
        nₜ = isclose(r) ? D[r.e] : C[r.e]
        nₕ = d
        insertnode!(c, nₜ, nₕ, r, s)
        # Step 2.3: Revise cost and selection vectors appropriately
        for (i,_) ∈ pairs(v.R) x[i,:] .= Inf end
        x[:,j] .= Inf
        for (k,v) ∈ pairs(V) ϕ[k] = isequal(r.o, v.i) ? 1 : 0 end
    end
    # Step 3: Return initial solution
    return s
end

# Random Initialization
# Create initial solution with randomly selcted node-route combination until all customer nodes have been added to the solution
function random(rng::AbstractRNG, G)
    s = Solution(G...)
    D = s.D
    C = s.C
    V = s.V
    # Step 1: Initialize an empty route for every vehicle
    M = typemax(Int64)
    for d ∈ D for v ∈ d.V push!(v.R, Route(rand(rng, 1:M), v, d)) end end
    R = [r for v ∈ V for r ∈ v.R]
    w = ones(Int64, eachindex(C))       # w[j]: selection weight for customer node C[j]
    # Step 2: Iteratively append randomly selected ncustomer ode in randomly selected route
    for _ ∈ eachindex(C)
        i = sample(rng, eachindex(R))
        j = sample(rng, eachindex(C), OffsetWeights(w))
        r = R[i]
        c = C[j]
        v = V[r.o]
        d = D[v.o]
        nₜ = d
        nₕ = isclose(r) ? D[r.s] : C[r.s]
        insertnode!(c, nₜ, nₕ, r, s)
        w[j] = 0
    end
    # Step 3: Return initial solution
    return s
end

# Regret-N Insertion
# Create initial solution by iteratively adding customer nodes with highest regret cost at its best position until all customer nodes have been added to the solution
function regretₙinit(rng::AbstractRNG, N::Int64, G)
    s = Solution(G...)
    D = s.D
    C = s.C
    V = s.V
    # Step 1: Initialize an empty route for every vehicle
    M = typemax(Int64)
    for d ∈ D for v ∈ d.V push!(v.R, Route(rand(rng, 1:M), v, d)) end end
    R = [r for v ∈ V for r ∈ v.R]
    I = length(R)
    J = eachindex(C)
    K = length(V)
    xᵢ = fill(Inf, (I,J))               # xᵢ[i,j]: insertion cost of customer node C[j] at best position in route R[i]
    pᵢ = fill(Int64.((0, 0)), (I,J))    # pᵢ[i,j]: best insertion postion of customer node C[j] in route R[i]
    xₙ = fill(Inf, (N,J))               # xₙ[i,n]: insertion cost of customer node C[j] at nᵗʰ best position
    rₙ = zeros(Int64, (N,J))            # rₙ[i,j]: N best insertion route index of customer node C[j]
    xᵣ = fill(-Inf, J)                  # x[j]   : regret-N cost of customer node C[j]
    ϕ  = ones(Int64, K)                 # ϕ[k]   : selection weight for vehicle V[k]
    # Step 2: Iterate until all customer nodes have been inserted into the route
    for _ ∈ eachindex(C)
        # Step 2.1: Iterate through all open customer nodes and every route
        z = f(s)
        for (j,c) ∈ pairs(C)
            if isclose(c) continue end
            for (i,r) ∈ pairs(R)
                # Step 2.1.1: Iterate through all possible insertion position in route r
                if iszero(ϕ[r.o]) continue end
                v  = V[r.o]
                d  = D[v.o]
                nₛ = isclose(r) ? D[r.s] : C[r.s]
                nₑ = isclose(r) ? D[r.e] : C[r.e]
                nₜ = d
                nₕ = nₛ
                while true
                    # Step 2.1.1.1: Insert customer node c between tail node nₜ and head node nₕ in route r, and compute the insertion cost
                    insertnode!(c, nₜ, nₕ, r, s)
                    # Step 2.1.1.2: Compute the insertion cost
                    z⁺ = f(s)
                    Δ  = z⁺ - z
                    # Step 2.1.1.3: Revise least insertion cost in route r and the corresponding best insertion position in route r
                    if Δ < xᵢ[i,j] xᵢ[i,j], pᵢ[i,j] = Δ, (nₜ.i, nₕ.i) end
                    # Step 2.1.1.4: Revise N least insertion costs
                    n̲ = 1
                    for n ∈ 1:N
                        n̲ = n
                        if Δ < xₙ[n,j] break end
                    end
                    for n ∈ N:-1:n̲ 
                        xₙ[n,j] = isequal(n, n̲) ? Δ : xₙ[n-1,j]::Float64
                        rₙ[n,j] = isequal(n, n̲) ? r.i : rₙ[n-1,j]::Int64
                    end
                    # Step 2.1.1.5: Remove customer node c from its position between tail node nₜ and head node nₕ in route r
                    removenode!(c, nₜ, nₕ, r, s)
                    if isequal(nₜ, nₑ) break end
                    nₜ = nₕ
                    nₕ = isequal(r.e, nₜ.i) ? D[nₜ.h] : C[nₜ.h]
                end
            end
            # Step 2.1.2: Compute regret cost for customer node c
            xᵣ[j] = 0.
            for n ∈ 1:N xᵣ[j] += xₙ[n,j] - xₙ[1,j] end
        end
        # Step 2.2: Insert customer node with highest regret cost in its best position (break ties by inserting the node with the lowest insertion cost)
        J̲ = findall(j -> j == maximum(xᵣ), xᵣ)
        i,j = Tuple(argmin(xᵢ[:,J̲]))
        j = J̲[j]
        r = R[i]
        c = C[j]
        v = V[r.o]
        t, h = pᵢ[i,j]
        nₜ = t ≤ length(D) ? D[t] : C[t]
        nₕ = h ≤ length(D) ? D[h] : C[h]
        insertnode!(c, nₜ, nₕ, r, s)
        # Step 2.3: Revise vectors appropriately
        for (k,r) ∈ pairs(R) if isequal(r.o, v.i) xᵢ[k,:] .= Inf end end
        for (k,r) ∈ pairs(R) if isequal(r.o, v.i) pᵢ[k,:] .= ((0, 0), ) end end
        xᵢ[:,j] .= Inf
        pᵢ[:,j] .= ((0, 0), )
        xₙ[:,j] .= Inf
        rₙ[:,j] .= 0
        xᵣ .= -Inf 
        for (k,v) ∈ pairs(V) ϕ[k] = isequal(r.o, v.i) ? 1 : 0 end
        for j ∈ eachindex(C) 
            for n ∈ 1:N
                if iszero(rₙ[n,j]) break end
                k = findfirst(x -> x.i == rₙ[n,j], R)
                r = R[k]
                if isequal(r.o, v.i) xₙ[n,j], rₙ[n,j] = Inf, 0 end
            end
            ix = sortperm(xₙ[:,j])
            xₙ[:,j] .= xₙ[ix,j]
            rₙ[:,j] .= rₙ[ix,j]
        end
    end
    # Step 3: Return initial solution
    return s
end
regret₂init(rng::AbstractRNG, G) = regretₙinit(rng, Int64(2), G)
regret₃init(rng::AbstractRNG, G) = regretₙinit(rng, Int64(3), G)
