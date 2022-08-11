"""
    initialsolution([rng], (N::Vector{Node}, A::Dict{Tuple{Int64, Int64}, Arc}), method)

Returns initial TSP solution using given `method` on graph with nodes `N` and arcs `A`.
Available methods include,
- Clarke and Wright Savings Algorithm   : `:cw_init`
- Nearest Neighborhood Algorithm        : `:nn_init`
- Random Initialization                 : `:random_init`
- Regret K Insertion                    : `:regret₂init`, `:regret₃init`

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
initialsolution(rng::AbstractRNG, G, method::Symbol)::Solution = getfield(LRP, method)(rng, G)
initialsolution(G, method::Symbol) = initialsolution(Random.GLOBAL_RNG, G, method)

# Clarke and Wright Savings Algorithm
# Create initial solution by merging routes that render the most savings until single route traversing all nodes remains
function cw_init(rng::AbstractRNG, G)
    N, A = G
    s = Solution(N, A)
    d = N[1]
    # Step 1: Initialize solution with routes to every node from the depot node (first node)
    K = length(N)
    for k ∈ K:-1:1 insertnode!(N[k], d, d, s) end
    # Step 2: Merge routes iteratively until single route traversing all nodes remains
    x = fill(-Inf, (K,K))       # x[i,j]: Savings from merging route with tail node N[i] into route with tail node N[j]
    ϕ = ones(Int64, K)          # ϕ[k]  : if isone(ϕ[k]) implies route k is active else inactive
    for _ ∈ 3:K
        # Step 2.1: Iterate through every route-pair combination
        z = f(s)
        for i ∈ 2:K
            n₁ = N[i]
            if !isequal(n₁.h, d.i) continue end
            for j ∈ 2:K
                # Step 2.1.1: Merge routes with tail node N[i] into route with tail node N[j]
                n₂ = N[j]
                if isequal(n₁, n₂) continue end
                if !isequal(n₂.h, d.i) continue end
                if iszero(ϕ[i]) & iszero(ϕ[j]) continue end
                q  = d
                nₒ = n₁
                p  = N[nₒ.t]
                while true
                    removenode!(nₒ, p, d, s)
                    insertnode!(nₒ, n₂,q, s)
                    if isequal(p, d) break end
                    q  = nₒ
                    nₒ = p
                    p  = N[nₒ.t]
                end
                # Step 2.1.2: Compute savings from merging route with tail node N[i] into route with tail node N[j]
                z⁻ = f(s)
                Δ  = z - z⁻
                x[i,j] = Δ
                # Step 2.1.3: Unmerge routes with tail node N[i] into route with tail node N[j]
                p  = d
                nₒ = N[n₂.h]
                q  = N[nₒ.h]
                while true
                    removenode!(nₒ, n₂,q, s)
                    insertnode!(nₒ, p, d, s)
                    if isequal(q, d) break end
                    p  = nₒ
                    nₒ = q
                    q  = N[nₒ.h]
                end
            end
        end
        # Step 2.2: Merge routes that render highest savings
        i,j = Tuple(argmax(x))
        n₁ = N[i]
        n₂ = N[j]
        q  = d
        nₒ = n₁
        p  = N[nₒ.t]
        while true
            removenode!(nₒ, p, d, s)
            insertnode!(nₒ, n₂,q, s)
            if isequal(p, d) break end
            q  = nₒ
            nₒ = p
            p  = N[nₒ.t]
        end
        # Step 2.3: Revise savings and selection vectors appropriately
        x[i, :] .= -Inf
        x[:, i] .= -Inf
        x[j, :] .= -Inf
        x[:, j] .= -Inf
        for k ∈ 1:K ϕ[k] = isequal(k, i) ? 1 : 0 end
    end
    for n ∈ N d.t = isequal(n.h, d.i) ? n.i : d.t end
    for n ∈ N d.h = isequal(n.t, d.i) ? n.i : d.h end
    # Step 3: Return initial solution
    return s
end

# Nearest Neighborhood Algorithm
# Create initial solution appending nodes that result in least increase in cost until all nodes have been added to the solution
function nn_init(rng::AbstractRNG, G)
    N, A = G
    s = Solution(N, A)
    N = s.N
    d = N[1]
    # Step 1: Initialize
    I = length(N)
    x = fill(Inf, I)            # x[i]: cost of appending node N[i] in the route
    # Step 2: Iterate until all nodes have been added to the route
    for _ ∈ 2:I
        # Step 2.1: Iteratively compute cost of appending each open node in the route
        z = f(s)
        for i ∈ 2:I
            # Step 2.2.1: Append node N[i] in the route
            nₒ = N[i]
            if isclose(nₒ) continue end
            nₜ = N[d.t]
            nₕ = d
            insertnode!(nₒ, nₜ, nₕ, s)
            # Step 2.2.2: Compute increase in cost
            z⁺ = f(s)
            Δ  = z⁺ - z
            x[i] = Δ
            # Step 2.2.3: Pop node N[i] from the route
            removenode!(nₒ, nₜ, nₕ, s)
        end
        # Step 2.2: Append open node N[i] in the route with least increase in cost
        i = argmin(x)
        nₒ = N[i]
        nₜ = N[d.t]
        nₕ = d
        insertnode!(nₒ, nₜ, nₕ, s)
        # Step 2.3: Close customer N[i]
        x[i] = Inf
    end
    # Step 3: Return initial solution
    return s
end

# Random Initialization
# Create initial solution by iteratively appending randomly selected node until all nodes have been added to the solution
function random_init(rng::AbstractRNG, G)
    N, A = G
    s = Solution(N, A)
    d = N[1]
    # Step 1: Intiialize
    I = length(N)
    w = ones(Int64, I)          # w[i]: selection weight for node N[i]
    w[1] = 0
    # Step 2: Iteratively append randomly selected open node
    for _ ∈ 2:I
        i = sample(rng, 1:I, Weights(w))
        nₒ = N[i]
        nₜ = N[d.t]
        nₕ = d
        insertnode!(nₒ, nₜ, nₕ, s)
        w[i] = 0
    end
    # Step 3: Return initial solution
    return s
end

# Regret-K Insertion
# Create initial solution by iteratively adding nodes with highest regret cost as its best position until all nodes have been added to the solution
function regretₖinit(rng::AbstractRNG, K, G)
    N, A = G
    s = Solution(N, A)
    d = N[1]
    # Step 1: Initialize
    I = length(N)
    p = fill((0, 0), I)         # p[i]   : best insertion postion of node N[i]
    x = fill(Inf, I)            # x[i]   : insertion cost of node N[i] at best position
    y = fill(Inf, (K,I))        # y[k,i] : insertion cost of node N[i] at kᵗʰ best position
    r = fill(-Inf, I)           # r[i]   : regret-K cost of node N[i]
    # Step 2: Iterate until all nodes have been inserted into the route
    for _ ∈ 2:I
        # Step 2.1: Iterate through all open nodes
        z = f(s)
        for i ∈ 2:I
            # Step 2.1.1: For node N[i] compute insertion cost for every possible insertion position
            nₒ = N[i]
            if isclose(nₒ) continue end
            nₜ = d
            nₕ = N[nₜ.h]
            while true
                # Step 2.1.1.1: Insert node between tail node nₜ and head node nₕ
                insertnode!(nₒ, nₜ, nₕ, s)
                # Step 2.1.1.2: Compute the insertion cost
                z⁺ = f(s)
                Δ  = z⁺ - z
                # Step 2.1.1.3: Revise least insertion cost and the corresponding best insertion position
                if Δ < x[i] x[i], p[i] = Δ, (nₜ.i, nₕ.i) end
                # Step 2.1.1.4: Revise K least insertion costs
                k̲ = 1
                for k ∈ 1:K
                    k̲ = k
                    if Δ < y[k,i] break end
                end
                for k ∈ K:-1:k̲ y[k,i] = isequal(k, k̲) ? Δ : y[k-1,i]::Float64 end
                # Step 2.1.1.5: Remove node from its position between tail node nₜ and head node nₕ
                removenode!(nₒ, nₜ, nₕ, s)
                if isequal(nₕ, d) break end
                nₜ = nₕ
                nₕ = N[nₜ.h]
            end
            # Step 2.1.2: Compute regret cost for node N[i]
            r[i] = 0.
            for k ∈ 1:K r[i] += y[k,i] - y[1,i] end
        end
        # Step 2.2: Insert node with highest regret cost in the best insertion position (break ties by inserting the node with the lowest insertion cost)
        I̲  = findall(i -> i == maximum(r), r)
        i  = I̲[argmin(x[I̲])]
        nₒ = N[i]
        t  = p[i][1]
        h  = p[i][2]
        nₜ = N[t]
        nₕ = N[h]
        insertnode!(nₒ, nₜ, nₕ, s)
        # Step 2.3: Revise vectors appropriately
        p .= ((0, 0), )
        x .= Inf
        y[:,:] .= Inf
        r .= -Inf
    end
    # Step 3: Return initial solution
    return s
end
regret₂init(rng::AbstractRNG, G) = regretₖinit(rng, 2, G)
regret₃init(rng::AbstractRNG, G) = regretₖinit(rng, 3, G)
