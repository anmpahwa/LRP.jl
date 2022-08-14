"""
    localsearch!(rng::AbstractRNG, k̅::Int, s::Solution, χₒ::ObjectiveFunctionParameters, method::Symbol)

Return solution `s` performing local seach on the solution using given `method` for `k̅` iterations 
until improvement. `χₒ` includes the objective function parameters for objective function evaluation.

Available methods include,
- Move  : `move!`
- Swap  : `swap!`

Optionally specify a random number generator `rng` as the first argument (defaults to `Random.GLOBAL_RNG`).
"""
localsearch!(rng::AbstractRNG, k̅::Int, s::Solution, χₒ::ObjectiveFunctionParameters, method::Symbol)::Solution = getfield(LRP, method)(rng, k̅, s, χₒ)
localsearch!(k̅::Int, s::Solution, χₒ::ObjectiveFunctionParameters, method::Symbol) = localsearch!(Random.GLOBAL_RNG, k̅, s, χₒ, method)

# Move
# Iteratively move a randomly seceted customer node in its best position if the move 
# results in reduction in objective function value for k̅ iterations until improvement
function move!(rng::AbstractRNG, k̅::Int, s::Solution, χₒ::ObjectiveFunctionParameters)
    z = f(s, χₒ)
    D = s.D
    C = s.C
    V = s.V
    R = [r for v ∈ V for r ∈ v.R]
    # Step 1: Initialize
    I = length(R)
    x = fill(Inf, I)                # x[i]: insertion cost in route R[i]
    p = fill(Int64.((0, 0)), I)     # p[i]: best insertion postion in route R[i]
    w = ones(Int64, eachindex(C))   # w[j]: selection weight for node C[j]
    # Step 2: Iterate for k̅ iterations until improvement
    for _ ∈ 1:k̅
        # Step 2.1: Randomly select a node
        j = sample(rng, eachindex(C), OffsetWeights(w))
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
                z′ = f(s, χₒ)
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

# Swap nodes
# Iteratively swap two randomly selected customer nodes if the swap results
# in reduction in objective function value for k̅ iterations until improvement
function swap!(rng::AbstractRNG, k̅::Int, s::Solution, χₒ::ObjectiveFunctionParameters)
    z = f(s, χₒ)
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
        z′ = f(s, χₒ)
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