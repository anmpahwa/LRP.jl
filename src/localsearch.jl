"""
    localsearch!(rng::AbstractRNG, k̅, s::Solution, method::Symbol)

Return solution `s` performing local seach on the solution using 
given `method` for `k̅` iterations until improvement.

Available methods include,
- Move  : `move!`
- 2-Opt : `opt!`
- Swap  : `swap!`

Optionally specify a random number generator `rng` as the first argument 
(defaults to `Random.GLOBAL_RNG`).
"""
localsearch!(rng::AbstractRNG, k̅, s::Solution, method::Symbol)::Solution = getfield(LRP, method)(rng, k̅, s)
localsearch!(k̅, s::Solution, method::Symbol) = localsearch!(Random.GLOBAL_RNG, k̅, s, method)

# Move
# Iteratively move a randomly seceted node in its best position if the move 
# results in reduction in objective function value for k̅ iterations until improvement
function move!(rng::AbstractRNG, k̅, s::Solution)
    z = f(s)
    N = s.N
    # Step 1: Initialize
    I = length(N)
    p = (0, 0)                  # p   : best insertion postion
    x = Inf                     # x   : insertion cost at best position
    w = ones(Int64, I)          # w[i]: selection weight for node N[i]
    # Step 2: Iterate for k̅ iterations until improvement
    for _ ∈ 1:k̅
        # Step 2.1: Randomly select a node
        i  = sample(rng, 1:I, Weights(w))
        nₒ = N[i]
        # Step 2.2: Remove this node from its position between tail node nₜ and head node nₕ
        nₜ = N[nₒ.t]
        nₕ = N[nₒ.h]
        pₒ = (nₜ.i, nₕ.i)
        removenode!(nₒ, nₜ, nₕ, s)
        nₜ = nₕ
        nₕ = N[nₜ.h]
        # Step 2.3: Iterate through all possible insertion positions
        while true
            # Step 2.3.1: Insert node between tail node nₜ and head node nₕ
            insertnode!(nₒ, nₜ, nₕ, s)
            # Step 2.3.2: Compute insertion cost
            z′ = f(s)
            Δ  = z′ - z
            # Step 2.3.3: Revise least insertion cost and the corresponding best insertion position
            if Δ < x x, p = Δ, (nₜ.i, nₕ.i) end
            # Step 2.3.4: Remove node from its position between tail node nₜ and head node nₕ
            removenode!(nₒ, nₜ, nₕ, s)
            if isequal((nₜ.i, nₕ.i), pₒ) break end
            nₜ = nₕ
            nₕ = N[nₜ.h]
        end
        # Step 2.4: Move the node to its best position (this could be its original position as well)
        Δ  = x
        t  = p[1]
        h  = p[2]
        nₜ = N[t]
        nₕ = N[h]
        insertnode!(nₒ, nₜ, nₕ, s)
        # Step 2.5: Revise vectors appropriately
        p = (0, 0)
        x = Inf
        w[i] = 0
        # Step 2.6: If the move results in reduction in objective function value, then go to step 3, else return to step 2.1
        Δ ≥ 0 ? continue : break
    end
    # Step 3: Return solution
    return s
end

# 2-opt: 
# Iteratively takes 2 arcs and reconfigure them if the reconfigure 
# results in reduction in objective function value (total possible 
# reconfigurations 2²-1 = 3) for k̅ iterations until improvement
function opt!(rng::AbstractRNG, k̅, s::Solution)
    z = f(s)
    N = s.N
    # Step 1: Iterate for k̅ iterations until improvement
    for _ ∈ 1:k̅
        # Step 1.1: Reconfigure two randomly selected arcs
        # n₁ → n₂ → n₃ and n₄ → n₅ → n₆ 
        n₂, n₅ = sample(rng, N), sample(rng, N)
        n₁ = N[n₂.t]
        n₃ = N[n₂.h]
        n₄ = N[n₅.t]
        n₆ = N[n₅.h]
        if isequal(n₂, n₅) || isequal(n₁, n₅) continue end 
        q  = n₆
        nₒ = n₂
        p  = n₃
        while true
            removenode!(nₒ, n₁, p, s)
            insertnode!(nₒ, n₅, q, s)
            q  = nₒ
            nₒ = p
            p  = N[p.h]
            if isequal(nₒ, n₅) break end
        end
        # Step 1.2: Compute change in objective function value
        z′ = f(s)
        Δ  = z′ - z 
        # Step 1.3: If the reconfigure results in reduction in objective function value then go to step 2, else go to step 1.4
        if Δ < 0 break end
        # Step 1.4: Reconfigure the two arcs to original state and go to step 1.1
        q  = n₆
        nₒ = n₅
        p  = n₄
        while true
            removenode!(nₒ, n₁, p, s)
            insertnode!(nₒ, n₂, q, s)
            q  = nₒ
            nₒ = p
            p  = N[p.h]
            if isequal(nₒ, n₂) break end
        end
    end
    # Step 2: Return solution
    return s
end

# Swap nodes
# Iteratively swap two randomly selected nodes if the swap results in 
# reduction in objective function value for k̅ iterations until improvement
function swap!(rng::AbstractRNG, k̅, s::Solution)
    N = s.N
    z = f(s)
    # Step 1: Iterate for k̅ iterations until improvement
    for _ ∈ 1:k̅
        # Step 1.1: Swap two randomly selected nodes
        # n₁ → n₂ → n₃ and n₄ → n₅ → n₆
        n₂, n₅ = sample(rng, N), sample(rng, N)
        if isequal(n₂, n₅) continue end
        n₁ = N[n₂.t]
        n₃ = N[n₂.h]
        n₄ = N[n₅.t]
        n₆ = N[n₅.h]
        # n₁ → n₂ (n₄) → n₃ (n₅) → n₆   ⇒   n₁ → n₃ (n₅) → n₂ (n₄) → n₆
        if isequal(n₃, n₅)
            removenode!(n₂, n₁, n₃, s)
            insertnode!(n₂, n₅, n₆, s)
        # n₄ → n₅ (n₁) → n₂ (n₆) → n₃   ⇒   n₄ → n₂ (n₆) → n₅ (n₁) → n₃   
        elseif isequal(n₂, n₆)
            removenode!(n₂, n₁, n₃, s)
            insertnode!(n₂, n₄, n₅, s)
        # n₁ → n₂ → n₃ and n₄ → n₅ → n₆ ⇒   n₁ → n₅ → n₃ and n₄ → n₂ → n₆
        else 
            removenode!(n₂, n₁, n₃, s)
            removenode!(n₅, n₄, n₆, s)
            insertnode!(n₅, n₁, n₃, s)
            insertnode!(n₂, n₄, n₆, s)
        end
        # Step 1.2: Compute change in objective function value
        z′ = f(s)
        Δ  = z′ - z 
        # Step 1.3: If the swap results in reduction in objective function value then go to step 2, else go to step 1.4
        if Δ < 0 break end
        # Step 1.4: Reswap the two nodes and go to step 1.1
        # n₁ → n₂ (n₄) → n₃ (n₅) → n₆   ⇒   n₁ → n₃ (n₅) → n₂ (n₄) → n₆
        if isequal(n₃, n₅)
            removenode!(n₂, n₅, n₆, s)
            insertnode!(n₂, n₁, n₃, s)
        # n₄ → n₅ (n₁) → n₂ (n₆) → n₃   ⇒   n₄ → n₂ (n₆) → n₅ (n₁) → n₃   
        elseif isequal(n₂, n₆)
            removenode!(n₂, n₄, n₅, s)
            insertnode!(n₂, n₁, n₃, s)
        # n₁ → n₂ → n₃ and n₄ → n₅ → n₆ ⇒   n₁ → n₅ → n₃ and n₄ → n₂ → n₆
        else 
            removenode!(n₅, n₁, n₃, s)
            removenode!(n₂, n₄, n₆, s)
            insertnode!(n₂, n₁, n₃, s)
            insertnode!(n₅, n₄, n₆, s)
        end
    end
    # Step 2: Return solution
    return s
end


        