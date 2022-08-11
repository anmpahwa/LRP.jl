"""
    insert!([rng], s::Solution, method::Symbol)

Return solution inserting open nodes to the solution `s` using the given `method`.

Available methods include,
- Best Insertion    : `best_insert!`
- Greedy Insertion  : `greedy_insert!`
- Regret Insertion  : `regret₂insert!`, `regret₃insert!`

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
insert!(rng::AbstractRNG, s::Solution, method::Symbol)::Solution = getfield(LRP, method)(rng, s)
insert!(s::Solution, method::Symbol) = insert!(Random.GLOBAL_RNG, s, method)

# Best insertion
# Iteratively insert randomly selected node at its best position until all open nodes have been added to the solution
function best_insert!(rng::AbstractRNG, s::Solution)
    N = s.N
    d = N[1]
    L = [n for n ∈ N if isopen(n)]
    # Step 1: Initialize
    I = length(L)
    p = fill((0, 0), I)         # p[i]: best insertion postion of node L[i]
    x = fill(Inf, I)            # x[i]: insertion cost of node L[i] at best position
    w = ones(Int64, I)          # w[i]: selection weight for node L[i]
    # Step 2: Iterate until all open nodes have been inserted into the route
    for _ ∈ 1:I
        # Step 2.1: Iterate through all open nodes and every possible insertion position
        z = f(s)
        for i ∈ 1:I
            nₒ = L[i]
            if isclose(nₒ) continue end
            nₜ = d
            nₕ = N[nₜ.h]
            while true
                # Step 2.1.1: Insert node between tail node nₜ and head node nₕ
                insertnode!(nₒ, nₜ, nₕ, s)
                # Step 2.1.2: Compute the insertion cost
                z⁺ = f(s)
                Δ  = z⁺ - z
                # Step 2.1.3: Revise least insertion cost and the corresponding best insertion position
                if Δ < x[i] x[i], p[i] = Δ, (nₜ.i, nₕ.i) end
                # Step 2.1.4: Remove node from its position between tail node nₜ and head node nₕ
                removenode!(nₒ, nₜ, nₕ, s)
                if isequal(nₕ, d) break end
                nₜ = nₕ
                nₕ = N[nₜ.h]
            end
        end
        # Step 2.2: Randomly select a node to insert at its best position
        i = sample(rng, 1:I, Weights(w))
        nₒ = L[i]
        t  = p[i][1]
        h  = p[i][2]
        nₜ = N[t]
        nₕ = N[h]
        insertnode!(nₒ, nₜ, nₕ, s)
        # Step 2.3: Revise vectors appropriately
        p .= ((0, 0), )
        x .= Inf
        w[i] = 0
    end
    return s
end

# Greedy insertion
# Iteratively insert nodes with least insertion cost at its best position until all open nodes have been added to the solution
function greedy_insert!(s::Solution)
    N = s.N
    d = N[1]
    L = [n for n ∈ N if isopen(n)]
    # Step 1: Initialize
    I = length(L)
    p = fill((0, 0), I)         # p[i]: best insertion postion of node L[i]
    x = fill(Inf, I)            # x[i]: insertion cost of node L[i] at best position
    # Step 2: Iterate until all open nodes have been inserted into the route
    for _ ∈ 1:I
        # Step 2.1: Iterate through all open nodes
        z = f(s)
        for i ∈ 1:I
            # Step 2.1.1: For open node L[i] compute insertion cost for every possible insertion position
            nₒ = L[i]
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
                # Step 2.1.1.4: Remove node from its position between tail node nₜ and head node nₕ
                removenode!(nₒ, nₜ, nₕ, s)
                if isequal(nₕ, d) break end
                nₜ = nₕ
                nₕ = N[nₜ.h]
            end
        end
        # Step 2.2: Insert node with least insertion cost at its best position
        i  = argmin(x)
        nₒ = L[i]
        t  = p[i][1]
        h  = p[i][2]
        nₜ = N[t]
        nₕ = N[h]
        insertnode!(nₒ, nₜ, nₕ, s)
        # Step 2.3: Revise vectors appropriately
        p .= ((0, 0), )
        x .= Inf
    end
    return s
end
greedy_insert!(rng::AbstractRNG, s::Solution) = greedy_insert!(s)

# Regret-K Insertion
# Iteratively add nodes with highest regret cost at its best position until all open nodes have been added to the solution
function regretₖinsert!(K, s::Solution)
    N = s.N
    d = N[1]
    L = [n for n ∈ N if isopen(n)]
    # Step 1: Initialize
    I = length(L)
    p = fill((0, 0), I)         # p[i]  : best insertion postion of node L[i]
    x = fill(Inf, I)            # x[i]  : insertion cost of node L[i] at best position
    y = fill(Inf, (K,I))        # y[k,i]: insertion cost of node L[i] at kᵗʰ best position
    r = fill(-Inf, I)           # r[i]  : regret-K cost of node L[i]
    # Step 2: Iterate until all open nodes have been inserted into the route
    for _ ∈ 1:I
        # Step 2.1: Iterate through all open nodes
        z = f(s)
        for i ∈ 1:I
            # Step 2.1.1: For node L[i] compute insertion cost for every possible insertion position
            nₒ = L[i]
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
            # Step 2.1.2: Compute regret cost for node L[i]
            r[i] = 0.
            for k ∈ 1:K r[i] += y[k,i] - y[1,i] end
        end
        # Step 2.2: Insert node with highest regret cost at its best position (break ties by inserting the node with the lowest insertion cost)
        I̲  = findall(i -> i == maximum(r), r)
        i  = I̲[argmin(x[I̲])]
        nₒ = L[i]
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
regret₂insert!(rng::AbstractRNG, s::Solution) = regretₖinsert!(2, s)
regret₃insert!(rng::AbstractRNG, s::Solution) = regretₖinsert!(3, s)
