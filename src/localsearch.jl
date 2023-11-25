"""
    localsearch!(rng::AbstractRNG, k̅::Int64, s::Solution, method::Symbol)

Returns solution `s` after performing local seach on the solution using given `method` for `k̅` iterations.

Available methods include,
- intra-Opt     : `:intraopt!`
- inter-Opt     : `:interopt!`
- intra-Move    : `:intramove!`
- inter-Move    : `:intramove!`
- split         : `:split!`
- swap          : `:swap!`

Optionally specify a random number generator `rng` as the first argument (defaults to `Random.GLOBAL_RNG`).
"""
localsearch!(rng::AbstractRNG, k̅::Int64, s::Solution, method::Symbol)::Solution = isdefined(LRP, method) ? getfield(LRP, method)(rng, k̅, s) : getfield(Main, method)(rng, k̅, s)
localsearch!(k̅::Int64, s::Solution, method::Symbol) = localsearch!(Random.GLOBAL_RNG, k̅, s, method)



"""
    intraopt!(rng::AbstractRNG, k̅::Int64, s::Solution)

Returns solution `s` after iteratively taking 2 arcs from the same route 
and reconfiguring them (total possible reconfigurations 2²-1 = 3) if the 
reconfiguration results in a reduction in objective function value, repeating 
for `k̅` iterations.
"""
function intraopt!(rng::AbstractRNG, k̅::Int64, s::Solution)
    prelocalsearch!(s)
    z = f(s)
    D = s.D
    C = s.C
    R = [r for d ∈ D for v ∈ d.V for r ∈ v.R]
    W = isopt.(R)                   # W[i]: selection weight for route R[i]
    # Step 1: Iterate for k̅ iterations
    for _ ∈ 1:k̅
        # Step 1.1: Iteratively take 2 arcs from the same route
        # d → ... → n¹ → n² → n³ → ... → n⁴ → n⁵ → n⁶ → ... → d
        r = sample(rng, R, Weights(W))
        (i,j) = sample(rng, 1:r.n, 2)
        (i,j) = j < i ? (j,i) : (i,j)  
        k  = 1
        c  = C[r.iˢ]
        n² = c
        n⁵ = c
        while true
            if isequal(k, i) n² = c end
            if isequal(k, j) n⁵ = c end
            if isequal(k, j) break end
            k += 1
            c  = C[c.iʰ]
        end
        n¹ = isequal(r.iˢ, n².iⁿ) ? D[n².iᵗ] : C[n².iᵗ]
        n³ = isequal(r.iᵉ, n².iⁿ) ? D[n².iʰ] : C[n².iʰ]
        n⁴ = isequal(r.iˢ, n⁵.iⁿ) ? D[n⁵.iᵗ] : C[n⁵.iᵗ]
        n⁶ = isequal(r.iᵉ, n⁵.iⁿ) ? D[n⁵.iʰ] : C[n⁵.iʰ] 
        if isequal(n², n⁵) || isequal(n¹, n⁵) continue end 
        # Step 1.2: Reconfigure
        # d → ... → n¹ → n⁵ → n⁴ → ... → n³ → n² → n⁶ → ... → d
        n  = n²
        tᵒ = n¹
        hᵒ = n³
        tⁿ = n⁵
        hⁿ = n⁶
        while true
            removenode!(n, tᵒ, hᵒ, r, s)
            insertnode!(n, tⁿ, hⁿ, r, s)
            hⁿ = n
            n  = hᵒ
            hᵒ = isdepot(hᵒ) ? C[r.iˢ] : (isequal(r.iᵉ, hᵒ.iⁿ) ? D[hᵒ.iʰ] : C[hᵒ.iʰ])
            if isequal(n, n⁵) break end
        end
        # Step 1.3: Compute change in objective function value
        z′ = f(s)
        Δ  = z′ - z 
        # Step 1.4: If the reconfiguration results in reduction in objective function value then go to step 1, else go to step 1.5
        if Δ < 0 z = z′
        # Step 1.5: Reconfigure back to the original state
        else
            # d → ... → n¹ → n² → n³ → ... → n⁴ → n⁵ → n⁶ → ... → d
            n  = n⁵
            tᵒ = n¹
            hᵒ = n⁴
            tⁿ = n²
            hⁿ = n⁶
            while true
                removenode!(n, tᵒ, hᵒ, r, s)
                insertnode!(n, tⁿ, hⁿ, r, s)
                hⁿ = n
                n  = hᵒ
                hᵒ = isdepot(hᵒ) ? C[r.iˢ] : (isequal(r.iᵉ, hᵒ.iⁿ) ? D[hᵒ.iʰ] : C[hᵒ.iʰ])
                if isequal(n, n²) break end
            end
        end
    end
    postlocalsearch!(s)
    # Step 2: Return solution
    return s
end



"""
    interopt!(rng::AbstractRNG, k̅::Int64, s::Solution)

Returns solution `s` after iteratively taking 2 arcs from the different 
routes and reconfiguring them (total possible reconfigurations 2²-1 = 3) 
if the reconfiguration results in a reduction in objective function value, 
repeating for `k̅` iterations.
"""
function interopt!(rng::AbstractRNG, k̅::Int64, s::Solution)
    prelocalsearch!(s)
    z = f(s)
    D = s.D
    C = s.C
    R = [r for d ∈ D for v ∈ d.V for r ∈ v.R]
    W = isopt.(R)                   # W[i]: selection weight for route R[i]
    # Step 1: Iterate for k̅ iterations
    for _ ∈ 1:k̅
        # Step 1.1: Iteratively take 2 arcs from different routes
        # d² → ... → n¹ → n² → n³ → ... → d² and d⁵ → ... → n⁴ → n⁵ → n⁶ → ... → d⁵
        r² = sample(rng, R, Weights(W))
        W′ = [relatedness(r², r⁵, s) * (!isequal(r², r⁵) * isopt(r⁵)) for r⁵ ∈ R]
        r⁵ = sample(rng, R, Weights(W′))
        d² = D[r².iᵈ]
        d⁵ = D[r⁵.iᵈ]
        i  = rand(rng, 1:r².n)
        k  = 1
        c² = C[r².iˢ]
        n² = c²
        while true
            if isequal(k, i) n² = c² end
            if isequal(k, i) break end
            k += 1
            c² = C[c².iʰ]
        end
        n¹ = isequal(r².iˢ, n².iⁿ) ? D[n².iᵗ] : C[n².iᵗ]
        n³ = isequal(r².iᵉ, n².iⁿ) ? D[n².iʰ] : C[n².iʰ]
        j  = rand(rng, 1:r⁵.n)
        k  = 1
        c⁵ = C[r⁵.iˢ]
        n⁵ = c⁵
        while true
            if isequal(k, j) n⁵ = c⁵ end
            if isequal(k, j) break end
            k += 1
            c⁵ = C[c⁵.iʰ]
        end
        n⁴ = isequal(r⁵.iˢ, n⁵.iⁿ) ? D[n⁵.iᵗ] : C[n⁵.iᵗ]
        n⁶ = isequal(r⁵.iᵉ, n⁵.iⁿ) ? D[n⁵.iʰ] : C[n⁵.iʰ]
        # Step 1.2: Reconfigure
        # d² → ... → n¹ → n⁵ → n⁶ → ...  → d² and d⁵ → ... → n⁴ → n² → n³ → ... → d⁵
        c² = n²
        tᵒ = n¹
        hᵒ = n³
        tⁿ = n⁴
        hⁿ = n⁵
        while true
            removenode!(c², tᵒ, hᵒ, r², s)
            insertnode!(c², tⁿ, hⁿ, r⁵, s)
            if isequal(hᵒ, d²) break end
            tⁿ = c² 
            c² = C[hᵒ.iⁿ]
            hᵒ = isequal(r².iᵉ, c².iⁿ) ? D[c².iʰ] : C[c².iʰ]
        end
        c⁵ = n⁵
        tᵒ = c²
        hᵒ = n⁶
        tⁿ = n¹
        hⁿ = d²
        while true
            removenode!(c⁵, tᵒ, hᵒ, r⁵, s)
            insertnode!(c⁵, tⁿ, hⁿ, r², s)
            if isequal(hᵒ, d⁵) break end
            tⁿ = c⁵
            c⁵ = C[hᵒ.iⁿ]
            hᵒ = isequal(r⁵.iᵉ, c⁵.iⁿ) ? D[c⁵.iʰ] : C[c⁵.iʰ]
        end
        # Step 1.3: Compute change in objective function value
        z′ = f(s)
        Δ  = z′ - z 
        # Step 1.4: If the reconfiguration results in reduction in objective function value then go to step 1, else go to step 1.5
        if Δ < 0 z = z′
        # Step 1.5: Reconfigure back to the original state
        else
            # d² → ... → n¹ → n² → n³ → ... → d² and d⁵ → ... → n⁴ → n⁵ → n⁶ → ... → d⁵
            c² = n⁵
            tᵒ = n¹
            hᵒ = isequal(r².iᵉ, c².iⁿ) ? D[c².iʰ] : C[c².iʰ]
            tⁿ = n⁴
            hⁿ = n²
            while true
                removenode!(c², tᵒ, hᵒ, r², s)
                insertnode!(c², tⁿ, hⁿ, r⁵, s)
                if isequal(hᵒ, d²) break end
                tⁿ = c² 
                c² = C[hᵒ.iⁿ]
                hᵒ = isequal(r².iᵉ, c².iⁿ) ? D[c².iʰ] : C[c².iʰ]
            end
            c⁵ = n²
            tᵒ = c²
            hᵒ = isequal(r⁵.iᵉ, c⁵.iⁿ) ? D[c⁵.iʰ] : C[c⁵.iʰ]
            tⁿ = n¹
            hⁿ = d²
            while true
                removenode!(c⁵, tᵒ, hᵒ, r⁵, s)
                insertnode!(c⁵, tⁿ, hⁿ, r², s)
                if isequal(hᵒ, d⁵) break end
                tⁿ = c⁵
                c⁵ = C[hᵒ.iⁿ]
                hᵒ = isequal(r⁵.iᵉ, c⁵.iⁿ) ? D[c⁵.iʰ] : C[c⁵.iʰ]
            end
        end
    end
    postlocalsearch!(s)
    # Step 2: Return solution
    return s
end



"""
    intramove!(rng::AbstractRNG, k̅::Int64, s::Solution)

Returns solution `s` after moving a randomly selected customer node 
to its best position in the same route if the move results in a reduction 
in objective function value, repeating for `k̅` iterations.
"""
function intramove!(rng::AbstractRNG, k̅::Int64, s::Solution)
    prelocalsearch!(s)
    D = s.D
    C = s.C
    # Step 1: Initialize
    for _ ∈ 1:k̅
        z = f(s)
        # Step 1.1: Randomly select a customer node
        c = sample(rng, C)
        # Step 1.2: Remove this node from its position between tail node nᵗ and head node nʰ
        r  = c.r
        nᵗ = isequal(r.iˢ, c.iⁿ) ? D[c.iᵗ] : C[c.iᵗ]
        nʰ = isequal(r.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ] 
        removenode!(c, nᵗ, nʰ, r, s)
        # Step 1.3: Iterate through all position in the route
        x  = 0.
        p  = (nᵗ.iⁿ, nʰ.iⁿ)
        d  = s.D[r.iᵈ]
        nˢ = isopt(r) ? C[r.iˢ] : D[r.iˢ] 
        nᵉ = isopt(r) ? C[r.iᵉ] : D[r.iᵉ]
        nᵗ = d
        nʰ = nˢ
        while true
            # Step 1.3.1: Insert customer node c between tail node nᵗ and head node nʰ
            insertnode!(c, nᵗ, nʰ, r, s)
            # Step 1.3.2: Compute insertion cost
            z′ = f(s)
            Δ  = z′ - z
            # Step 1.3.3: Revise least insertion cost in route r and the corresponding best insertion position in route r
            if Δ < x x, p = Δ, (nᵗ.iⁿ, nʰ.iⁿ) end
            # Step 1.3.4: Remove node from its position between tail node nᵗ and head node nʰ
            removenode!(c, nᵗ, nʰ, r, s)
            if isequal(nᵗ, nᵉ) break end
            nᵗ = nʰ
            nʰ = isequal(r.iᵉ, nᵗ.iⁿ) ? D[nᵗ.iʰ] : C[nᵗ.iʰ]
        end
        # Step 1.4: Move the node to its best position (this could be its original position as well)
        iᵗ = p[1]
        iʰ = p[2]
        nᵗ = iᵗ ≤ length(D) ? D[iᵗ] : C[iᵗ]
        nʰ = iʰ ≤ length(D) ? D[iʰ] : C[iʰ]
        insertnode!(c, nᵗ, nʰ, r, s)
    end
    postlocalsearch!(s)
    # Step 2: Return solution
    return s
end



"""
    intermove!(rng::AbstractRNG, k̅::Int64, s::Solution)

Returns solution `s` after moving a randomly selected customer node 
to its best position in another route if the move results in a reduction 
in objective function value, repeating for `k̅` iterations.
"""
function intermove!(rng::AbstractRNG, k̅::Int64, s::Solution)
    prelocalsearch!(s)
    D = s.D
    C = s.C
    # Step 1: Initialize
    R = [r for d ∈ D for v ∈ d.V for r ∈ v.R]
    for _ ∈ 1:k̅
        z = f(s)
        # Step 1.1: Randomly select a customer node
        c  = sample(rng, C)
        # Step 1.2: Remove this node from its position between tail node nᵗ and head node nʰ
        r¹ = c.r
        nᵗ = isequal(r¹.iˢ, c.iⁿ) ? D[c.iᵗ] : C[c.iᵗ]
        nʰ = isequal(r¹.iᵉ, c.iⁿ) ? D[c.iʰ] : C[c.iʰ] 
        removenode!(c, nᵗ, nʰ, r¹, s)
        # Step 1.3: Select a random route
        W  = [!isequal(r¹, r²) for r² ∈ R]
        r² = sample(rng, R, Weights(W))
        # Step 1.4: Iterate through all position in the route
        x  = 0.
        p  = (nᵗ.iⁿ, nʰ.iⁿ)
        r  = r¹
        d  = s.D[r².iᵈ]
        nˢ = isopt(r²) ? C[r².iˢ] : D[r².iˢ] 
        nᵉ = isopt(r²) ? C[r².iᵉ] : D[r².iᵉ]
        nᵗ = d
        nʰ = nˢ
        while true
            # Step 1.4.1: Insert customer node c between tail node nᵗ and head node nʰ
            insertnode!(c, nᵗ, nʰ, r², s)
            # Step 1.4.2: Compute insertion cost
            z′ = f(s)
            Δ  = z′ - z
            # Step 1.4.3: Revise least insertion cost in route r and the corresponding best insertion position in route r
            if Δ < x x, p, r = Δ, (nᵗ.iⁿ, nʰ.iⁿ), r² end
            # Step 1.4.4: Remove node from its position between tail node nᵗ and head node nʰ
            removenode!(c, nᵗ, nʰ, r², s)
            if isequal(nᵗ, nᵉ) break end
            nᵗ = nʰ
            nʰ = isequal(r².iᵉ, nᵗ.iⁿ) ? D[nᵗ.iʰ] : C[nᵗ.iʰ]
        end
        # Step 1.5: Move the node to its best position (this could be its original position as well)
        iᵗ = p[1]
        iʰ = p[2]
        nᵗ = iᵗ ≤ length(D) ? D[iᵗ] : C[iᵗ]
        nʰ = iʰ ≤ length(D) ? D[iʰ] : C[iʰ]
        insertnode!(c, nᵗ, nʰ, r, s)
    end
    postlocalsearch!(s)
    # Step 2: Return solution
    return s
end



"""
    split!(rng::AbstractRNG, k̅::Int64, s::Solution)

Returns solution `s` after iteratively spliting routes by moving a randomly 
selected depot node at best position if the split results in reduction in 
objective function value, repeating for `k̅` iterations.
"""
function split!(rng::AbstractRNG, k̅::Int64, s::Solution)
    prelocalsearch!(s)
    D = s.D
    C = s.C
    W = isopt.(D)                   # W[iᵈ]: selection weight for depot node D[iᵈ]
    # Step 1: Iterate for k̅ iterations
    for _ ∈ 1:k̅
        z = f(s)
        # Step 1.1: Select a random depot node d
        i = sample(rng, eachindex(D), Weights(W))
        d = D[i]
        # Step 1.2: Select a random route originating from this depot node
        R = [r for v ∈ d.V for r ∈ v.R]
        r = sample(rng, R, Weights(isopt.(R)))
        # Step 1.3: Remove depot node d from its position in route r
        cˢ = C[r.iˢ]
        cᵉ = C[r.iᵉ]
        removenode!(d, cᵉ, cˢ, r, s)
        # Step 1.4: Iterate through all possible positions in route r
        x = 0.
        p = (cᵉ.iⁿ, cˢ.iⁿ)
        cᵗ = cˢ
        cʰ = C[cᵗ.iʰ]
        while true
            # Step 1.4.1: Insert depot node d between tail node nᵗ and head node nʰ
            insertnode!(d, cᵗ, cʰ, r, s)
            # Step 1.4.2: Compute change in objective function value
            z′ = f(s) 
            Δ  = z′ - z
            # Step 1.4.3: Revise least insertion cost in route r and the corresponding best insertion position in route r
            if Δ < x x, p = Δ, (cᵗ.iⁿ, cʰ.iⁿ) end
            # Step 1.4.4: Remove depot node d from its position between tail node nᵗ and head node nʰ
            removenode!(d, cᵗ, cʰ, r, s)
            if isequal(cʰ, cᵉ) break end
            cᵗ = cʰ
            cʰ = C[cᵗ.iʰ]
        end
        # Step 1.5: Move the depot node to its best position in route r (this could be its original position as well)
        iᵗ = p[1]
        iʰ = p[2]
        cᵗ = C[iᵗ]
        cʰ = C[iʰ]
        insertnode!(d, cᵗ, cʰ, r, s) 
    end
    postlocalsearch!(s)
    # Step 2: Return solution
    return s
end



"""
    swap!(rng::AbstractRNG, k̅::Int64, s::Solution)

Returns solution `s` after swapping two randomly selected 
customers if the swap results in a reduction in objective 
function value, repeating for `k̅` iterations.
"""
function swap!(rng::AbstractRNG, k̅::Int64, s::Solution)
    prelocalsearch!(s)
    z = f(s)
    D = s.D
    C = s.C
    # Step 1: Iterate for k̅ iterations
    for _ ∈ 1:k̅
        # Step 1.1: Swap two randomly selected customer nodes
        # n¹ → n² → n³ and n⁴ → n⁵ → n⁶
        n² = sample(rng, C)
        W  = [isequal(n², n⁵) ? 0. : relatedness(n², n⁵, s) for n⁵ ∈ C]
        n⁵ = sample(rng, C, OffsetWeights(W))
        r², r⁵ = n².r, n⁵.r
        n¹ = isequal(r².iˢ, n².iⁿ) ? D[n².iᵗ] : C[n².iᵗ]
        n³ = isequal(r².iᵉ, n².iⁿ) ? D[n².iʰ] : C[n².iʰ]
        n⁴ = isequal(r⁵.iˢ, n⁵.iⁿ) ? D[n⁵.iᵗ] : C[n⁵.iᵗ]
        n⁶ = isequal(r⁵.iᵉ, n⁵.iⁿ) ? D[n⁵.iʰ] : C[n⁵.iʰ]
        # n¹ → n² (n⁴) → n³ (n⁵) → n⁶   ⇒   n¹ → n³ (n⁵) → n² (n⁴) → n⁶
        if isequal(n³, n⁵)
            removenode!(n², n¹, n³, r², s)
            insertnode!(n², n⁵, n⁶, r⁵, s)
        # n⁴ → n⁵ (n¹) → n² (n⁶) → n³   ⇒   n⁴ → n² (n⁶) → n⁵ (n¹) → n³   
        elseif isequal(n², n⁶)
            removenode!(n², n¹, n³, r², s)
            insertnode!(n², n⁴, n⁵, r⁵, s)
        # n¹ → n² → n³ and n⁴ → n⁵ → n⁶ ⇒   n¹ → n⁵ → n³ and n⁴ → n² → n⁶
        else 
            removenode!(n², n¹, n³, r², s)
            removenode!(n⁵, n⁴, n⁶, r⁵, s)
            insertnode!(n⁵, n¹, n³, r², s)
            insertnode!(n², n⁴, n⁶, r⁵, s)
        end
        # Step 1.2: Compute change in objective function value
        z′ = f(s)
        Δ  = z′ - z
        # Step 1.3: If the swap results in reduction in objective function value then go to step 1, else go to step 1.4
        if Δ < 0 z = z′
        # Step 1.4: Reswap the two customer nodes and go to step 1.1
        else
            # n¹ → n² (n⁴) → n³ (n⁵) → n⁶   ⇒   n¹ → n³ (n⁵) → n² (n⁴) → n⁶
            if isequal(n³, n⁵)
                removenode!(n², n⁵, n⁶, r⁵, s)
                insertnode!(n², n¹, n³, r², s)
            # n⁴ → n⁵ (n¹) → n² (n⁶) → n³   ⇒   n⁴ → n² (n⁶) → n⁵ (n¹) → n³   
            elseif isequal(n², n⁶)
                removenode!(n², n⁴, n⁵, r⁵, s)
                insertnode!(n², n¹, n³, r², s)
            # n¹ → n² → n³ and n⁴ → n⁵ → n⁶ ⇒   n¹ → n⁵ → n³ and n⁴ → n² → n⁶
            else 
                removenode!(n⁵, n¹, n³, r², s)
                removenode!(n², n⁴, n⁶, r⁵, s)
                insertnode!(n², n¹, n³, r², s)
                insertnode!(n⁵, n⁴, n⁶, r⁵, s)
            end
        end
    end
    postlocalsearch!(s)
    # Step 2: Return solution
    return s
end