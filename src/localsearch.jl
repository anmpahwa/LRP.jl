"""
    localsearch!(rng::AbstractRNG, k̅::Int64, s::Solution, method::Symbol)

Return solution `s` performing local seach on the solution using given `method` for `k̅` iterations 
until improvement.

Available methods include,
- Move      : `:move!`
- Inter-Opt : `:interopt!`
- Intra-Opt : `:intraopt!`
- Split     : `:split!`
- Swap      : `:swap!`

Optionally specify a random number generator `rng` as the first argument (defaults to `Random.GLOBAL_RNG`).
"""
localsearch!(rng::AbstractRNG, k̅::Int64, s::Solution, method::Symbol)::Solution = getfield(LRP, method)(rng, k̅, s)
localsearch!(k̅::Int64, s::Solution, method::Symbol) = localsearch!(Random.GLOBAL_RNG, k̅, s, method)

# Move
# Iteratively move a randomly seceted customer node in its best position if the move 
# results in reduction in objective function value for k̅ iterations until improvement
function move!(rng::AbstractRNG, k̅::Int64, s::Solution)
    z = f(s)
    D = s.D
    C = s.C
    V = s.V
    R = [r for v ∈ V for r ∈ v.R]
    # Step 1: Initialize
    I = eachindex(R)
    J = eachindex(C)
    x = fill(Inf, I)                # x[i]: insertion cost in route R[i]
    p = fill((0, 0), I)             # p[i]: best insertion postion in route R[i]
    w = ones(Int64, J)              # w[j]: selection weight for node C[j]
    # Step 2: Iterate for k̅ iterations until improvement
    for _ ∈ 1:k̅
        # Step 2.1: Randomly select a node
        j = sample(rng, J, OffsetWeights(w))
        c = C[j]
        # Step 2.2: Remove this node from its position between tail node nₜ and head node nₕ
        r = c.r
        nₜ = isequal(r.s, c.i) ? D[c.t] : C[c.t]
        nₕ = isequal(r.e, c.i) ? D[c.h] : C[c.h]
        removenode!(c, nₜ, nₕ, r, s)
        # Step 2.3: Iterate through all routes
        for (i,r) ∈ pairs(R)
            # Step 2.3.1: Iterate through all possible insertion positions
            v = V[r.o]
            d = D[v.o]
            nₛ = isclose(r) ? D[r.s] : C[r.s]
            nₑ = isclose(r) ? D[r.e] : C[r.e]
            nₜ = d
            nₕ = nₛ
            while true
                # Step 2.3.1.1: Insert customer node c between tail node nₜ and head node nₕ
                insertnode!(c, nₜ, nₕ, r, s)
                # Step 2.3.1.2: Compute insertion cost
                z′ = f(s)
                Δ  = z′ - z
                # Step 2.3.1.3: Revise least insertion cost in route r and the corresponding best insertion position in route r
                if Δ < x[i] x[i], p[i] = Δ, (nₜ.i, nₕ.i) end
                # Step 2.3.4: Remove node from its position between tail node nₜ and head node nₕ
                removenode!(c, nₜ, nₕ, r, s)
                if isequal(nₜ, nₑ) break end
                nₜ = nₕ
                nₕ = isequal(r.e, nₜ.i) ? D[nₜ.h] : C[nₜ.h]
            end
        end
        # Step 2.4: Move the node to its best position (this could be its original position as well)
        i = argmin(x)
        Δ = x[i]
        r = R[i]
        t = p[i][1]
        h = p[i][2]
        nₜ = t ≤ length(D) ? D[t] : C[t]
        nₕ = h ≤ length(D) ? D[h] : C[h]
        insertnode!(c, nₜ, nₕ, r, s)
        # Step 2.5: Revise vectors appropriately
        x .= Inf
        p .= ((0, 0), )
        w[j] = 0
        # Step 2.6: If the move results in reduction in objective function value, then go to step 3, else return to step 2.1
        Δ ≥ 0 ? continue : break
    end
    # Step 3: Return solution
    return s
end

# 2-opt
# Iteratively take 2 arcs and reconfigure them (total possible reconfigurations 2²-1 = 3) if the 
# reconfigure results in reduction in objective function value for k̅ iterations until improvement
function intraopt!(rng::AbstractRNG, k̅::Int64, s::Solution)
    z = f(s)
    D = s.D
    C = s.C
    V = s.V
    R = [r for v ∈ V for r ∈ v.R]
    w = isopen.(R)
    # Step 1: Iterate for k̅ iterations until improvement
    for _ ∈ 1:k̅
        # Step 1.1: Iteratively take 2 arcs from the same route
        # d → ... → n₁ → n₂ → n₃ → ... → n₄ → n₅ → n₆ → ... → d
        r = sample(rng, R, Weights(w))
        (i,j) = sample(rng, 1:r.n, 2)
        (i,j) = j < i ? (j,i) : (i,j)  
        k = 1
        c = C[r.s]
        n₂ = c
        n₅ = c
        while true
            if isequal(k, i) n₂ = c end
            if isequal(k, j) n₅ = c end
            if isequal(k, j) break end
            k += 1
            c = C[c.h]
        end
        n₁ = isequal(r.s, n₂.i) ? D[n₂.t] : C[n₂.t]
        n₃ = isequal(r.e, n₂.i) ? D[n₂.h] : C[n₂.h]
        n₄ = isequal(r.s, n₅.i) ? D[n₅.t] : C[n₅.t]
        n₆ = isequal(r.e, n₅.i) ? D[n₅.h] : C[n₅.h]
        if isequal(n₂, n₅) || isequal(n₁, n₅) continue end 
        # Step 1.2: Reconfigure
        # d → ... → n₁ → n₅ → n₄ → ... → n₃ → n₂ → n₆ → ... → d
        n  = n₂
        tₒ = n₁
        hₒ = n₃
        tₙ = n₅
        hₙ = n₆
        while true
            removenode!(n, tₒ, hₒ, r, s)
            insertnode!(n, tₙ, hₙ, r, s)
            hₙ = n
            n  = hₒ
            hₒ = isdepot(hₒ) ? C[r.s] : (isequal(r.e, hₒ.i) ? D[hₒ.h] : C[hₒ.h])
            if isequal(n, n₅) break end
        end
        # Step 1.3: Compute change in objective function value
        z′ = f(s)
        Δ  = z′ - z 
        # Step 1.4: If the reconfiguration results in reduction in objective function value then go to step 2, else go to step 1.5
        if Δ < 0 return s end
        # Step 1.5: Reconfigure back to the original state
        # d → ... → n₁ → n₂ → n₃ → ... → n₄ → n₅ → n₆ → ... → d
        n  = n₅
        tₒ = n₁
        hₒ = n₄
        tₙ = n₂
        hₙ = n₆
        while true
            removenode!(n, tₒ, hₒ, r, s)
            insertnode!(n, tₙ, hₙ, r, s)
            hₙ = n
            n  = hₒ
            hₒ = isdepot(hₒ) ? C[r.s] : (isequal(r.e, hₒ.i) ? D[hₒ.h] : C[hₒ.h])
            if isequal(n, n₂) break end
        end
    end
    # Step 2: Return solution
    return s
end
function interopt!(rng::AbstractRNG, k̅::Int64, s::Solution)
    z = f(s)
    D = s.D
    C = s.C
    V = s.V
    R = [r for v ∈ V for r ∈ v.R]
    w = isopen.(R)
    # Step 1: Iterate for k̅ iterations until improvement
    for _ ∈ 1:k̅
        # Step 1.1: Iteratively take 2 arcs from different routes
        # d₂ → ... → n₁ → n₂ → n₃ → ... → d₂ and d₅ → ... → n₄ → n₅ → n₆ → ... → d₅
        r₂, r₅ = sample(rng, R, Weights(w), 2)
        v₂, v₅ = V[r₂.o], V[r₅.o]
        d₂, d₅ = D[v₂.o], D[v₅.o]
        if isequal(r₂, r₅) continue end
        i = rand(rng, 1:r₂.n)
        k = 1
        c₂ = C[r₂.s]
        n₂ = c₂
        while true
            if isequal(k, i) n₂ = c₂ end
            if isequal(k, i) break end
            k += 1
            c₂ = C[c₂.h]
        end
        n₁ = isequal(r₂.s, n₂.i) ? D[n₂.t] : C[n₂.t]
        n₃ = isequal(r₂.e, n₂.i) ? D[n₂.h] : C[n₂.h]
        j = rand(rng, 1:r₅.n)
        k = 1
        c₅ = C[r₅.s]
        n₅ = c₅
        while true
            if isequal(k, j) n₅ = c₅ end
            if isequal(k, j) break end
            k += 1
            c₅ = C[c₅.h]
        end
        n₄ = isequal(r₅.s, n₅.i) ? D[n₅.t] : C[n₅.t]
        n₆ = isequal(r₅.e, n₅.i) ? D[n₅.h] : C[n₅.h]
        # Step 1.2: Reconfigure
        # d₂ → ... → n₁ → n₅ → n₆ → ...  → d₂ and d₅ → ... → n₄ → n₂ → n₃ → ... → d₅
        c₂ = n₂
        tₒ = n₁
        hₒ = n₃
        tₙ = n₄
        hₙ = n₅
        while true
            removenode!(c₂, tₒ, hₒ, r₂, s)
            insertnode!(c₂, tₙ, hₙ, r₅, s)
            if isequal(hₒ, d₂) break end
            tₙ = c₂ 
            c₂ = C[hₒ.i]
            hₒ = isequal(r₂.e, c₂.i) ? D[c₂.h] : C[c₂.h]
        end
        c₅ = n₅
        tₒ = c₂
        hₒ = n₆
        tₙ = n₁
        hₙ = d₂
        while true
            removenode!(c₅, tₒ, hₒ, r₅, s)
            insertnode!(c₅, tₙ, hₙ, r₂, s)
            if isequal(hₒ, d₅) break end
            tₙ = c₅
            c₅ = C[hₒ.i]
            hₒ = isequal(r₅.e, c₅.i) ? D[c₅.h] : C[c₅.h]
        end
        # Step 1.3: Compute change in objective function value
        z′ = f(s)
        Δ  = z′ - z 
        # Step 1.4: If the reconfiguration results in reduction in objective function value then go to step 2, else go to step 1.5
        if Δ < 0 break end
        # Step 1.5: Reconfigure back to the original state
        # d₂ → ... → n₁ → n₂ → n₃ → ... → d₂ and d₅ → ... → n₄ → n₅ → n₆ → ... → d₅
        c₂ = n₅
        tₒ = n₁
        hₒ = isequal(r₂.e, c₂.i) ? D[c₂.h] : C[c₂.h]
        tₙ = n₄
        hₙ = n₂
        while true
            removenode!(c₂, tₒ, hₒ, r₂, s)
            insertnode!(c₂, tₙ, hₙ, r₅, s)
            if isequal(hₒ, d₂) break end
            tₙ = c₂ 
            c₂ = C[hₒ.i]
            hₒ = isequal(r₂.e, c₂.i) ? D[c₂.h] : C[c₂.h]
        end
        c₅ = n₂
        tₒ = c₂
        hₒ = isequal(r₅.e, c₅.i) ? D[c₅.h] : C[c₅.h]
        tₙ = n₁
        hₙ = d₂
        while true
            removenode!(c₅, tₒ, hₒ, r₅, s)
            insertnode!(c₅, tₙ, hₙ, r₂, s)
            if isequal(hₒ, d₅) break end
            tₙ = c₅
            c₅ = C[hₒ.i]
            hₒ = isequal(r₅.e, c₅.i) ? D[c₅.h] : C[c₅.h]
        end
    end
    # Step 2: Return solution
    return s
end

# Split
# Iteratively split routes by moving a randomly selected depot node at best position if the
# split results in reduction in objective function value for k̅ iterations until improvement
function split!(rng::AbstractRNG, k̅::Int64, s::Solution)
    z = f(s)
    z̅ = z
    D = s.D
    C = s.C
    K = eachindex(D)
    w = isopen.(D)
    # Step 1: Iterate for k̅ iterations until improvement
    for _ ∈ 1:k̅
        # Step 1.1: Select a random depot node d
        k = sample(rng, K, Weights(w))
        d = D[k]
        V = d.V
        # Step 1.2: Iterate through every route originating from this depot node
        for v ∈ V
            R = v.R
            for r ∈ R
                # Step 1.2.1: Remove depot node d from its position in route r
                if isclose(r) continue end
                cₛ = C[r.s]
                cₑ = C[r.e]
                x = 0.
                p = (cₑ.i, cₛ.i)
                removenode!(d, cₑ, cₛ, r, s)
                # Step 1.2.2: Iterate through all possible positions in route r
                cₜ = cₛ
                cₕ = C[cₜ.h]
                while true
                    # Step 1.2.2.1: Insert depot node d between tail node nₜ and head node nₕ
                    insertnode!(d, cₜ, cₕ, r, s)
                    # Step 1.2.2.2: Compute change in objective function value
                    z′ = f(s) 
                    Δ  = z′ - z
                    # Step 1.2.2.3: Revise least insertion cost in route r and the corresponding best insertion position in route r
                    if Δ < x x, p = Δ, (cₜ.i, cₕ.i) end
                    # Step 1.2.2.4: Remove depot node d from its position between tail node nₜ and head node nₕ
                    removenode!(d, cₜ, cₕ, r, s)
                    if isequal(cₕ, cₑ) break end
                    cₜ = cₕ
                    cₕ = C[cₜ.h]
                end
                # Step 1.2.3: Move the depot node to its best position in route r (this could be its original position as well)
                t = p[1]
                h = p[2]
                cₜ = C[t]
                cₕ = C[h]
                insertnode!(d, cₜ, cₕ, r, s)
                z = f(s) 
            end
        end
        # Step 1.3: Revise vectors appropriately
        k = d.i
        w[k] = 0
        Δ = z - z̅
        # Step 1.4: If the overall change results in reduction in objective function value, then go to step 2, else return to step 1.1
        Δ ≥ 0 ? continue : break
    end
    # Step 2: Return solution
    return s
end

# Swap nodes
# Iteratively swap two randomly selected customer nodes if the swap results
# in reduction in objective function value for k̅ iterations until improvement
function swap!(rng::AbstractRNG, k̅::Int64, s::Solution)
    z = f(s)
    D = s.D
    C = s.C
    # Step 1: Iterate for k̅ iterations until improvement
    for _ ∈ 1:k̅
        # Step 1.1: Swap two randomly selected customer nodes
        # n₁ → n₂ → n₃ and n₄ → n₅ → n₆
        n₂, n₅ = sample(rng, C), sample(rng, C)
        if isequal(n₂, n₅) continue end
        r₂, r₅ = n₂.r, n₅.r
        n₁ = isequal(r₂.s, n₂.i) ? D[n₂.t] : C[n₂.t]
        n₃ = isequal(r₂.e, n₂.i) ? D[n₂.h] : C[n₂.h]
        n₄ = isequal(r₅.s, n₅.i) ? D[n₅.t] : C[n₅.t]
        n₆ = isequal(r₅.e, n₅.i) ? D[n₅.h] : C[n₅.h]
        # n₁ → n₂ (n₄) → n₃ (n₅) → n₆   ⇒   n₁ → n₃ (n₅) → n₂ (n₄) → n₆
        if isequal(n₃, n₅)
            removenode!(n₂, n₁, n₃, r₂, s)
            insertnode!(n₂, n₅, n₆, r₅, s)
        # n₄ → n₅ (n₁) → n₂ (n₆) → n₃   ⇒   n₄ → n₂ (n₆) → n₅ (n₁) → n₃   
        elseif isequal(n₂, n₆)
            removenode!(n₂, n₁, n₃, r₂, s)
            insertnode!(n₂, n₄, n₅, r₅, s)
        # n₁ → n₂ → n₃ and n₄ → n₅ → n₆ ⇒   n₁ → n₅ → n₃ and n₄ → n₂ → n₆
        else 
            removenode!(n₂, n₁, n₃, r₂, s)
            removenode!(n₅, n₄, n₆, r₅, s)
            insertnode!(n₅, n₁, n₃, r₂, s)
            insertnode!(n₂, n₄, n₆, r₅, s)
        end
        # Step 1.2: Compute change in objective function value
        z′ = f(s)
        Δ  = z′ - z 
        # Step 1.3: If the swap results in reduction in objective function value then go to step 2, else go to step 1.4
        if Δ < 0 break end
        # Step 1.4: Reswap the two customer nodes and go to step 1.1
        # n₁ → n₂ (n₄) → n₃ (n₅) → n₆   ⇒   n₁ → n₃ (n₅) → n₂ (n₄) → n₆
        if isequal(n₃, n₅)
            removenode!(n₂, n₅, n₆, r₅, s)
            insertnode!(n₂, n₁, n₃, r₂, s)
        # n₄ → n₅ (n₁) → n₂ (n₆) → n₃   ⇒   n₄ → n₂ (n₆) → n₅ (n₁) → n₃   
        elseif isequal(n₂, n₆)
            removenode!(n₂, n₄, n₅, r₅, s)
            insertnode!(n₂, n₁, n₃, r₂, s)
        # n₁ → n₂ → n₃ and n₄ → n₅ → n₆ ⇒   n₁ → n₅ → n₃ and n₄ → n₂ → n₆
        else 
            removenode!(n₅, n₁, n₃, r₂, s)
            removenode!(n₂, n₄, n₆, r₅, s)
            insertnode!(n₂, n₁, n₃, r₂, s)
            insertnode!(n₅, n₄, n₆, r₅, s)
        end
    end
    # Step 2: Return solution
    return s
end