"""
    insert!([rng::AbstractRNG], s::Solution, method::Symbol)

Returns solution `s` after inserting open customer nodes to the solution using the given `method`.

Available methods include,
- Best Insertion           : `:best!`
- Precise Greedy Insertion : `:precise!`
- Perturb Greedy insertion : `:perturb!`
- Regret-two Insertion     : `:regret2!`
- Regret-three Insertion   : `:regret3!`

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
insert!(rng::AbstractRNG, s::Solution, method::Symbol)::Solution = isdefined(LRP, method) ? getfield(LRP, method)(rng, s) : getfield(Main, method)(rng, s)
insert!(s::Solution, method::Symbol) = insert!(Random.GLOBAL_RNG, s, method)



"""
    best!(rng::AbstractRNG, s::Solution)

Returns solution `s` after inserting randomly selected customer node 
at its best position until all open nodes have been inserted to the 
solution.
"""
function best!(rng::AbstractRNG, s::Solution)
    if all(isclose, s.C) return s end
    # Step 1: Initialize
    preinsert!(s)
    D = s.D
    C = s.C
    R = [r for d ∈ D for v ∈ d.V for r ∈ v.R if isactive(r)]
    L = [c for c ∈ C if isopen(c)]
    I = eachindex(L)
    J = eachindex(R)
    W = ones(Int, I)                        # W[j]  : selection weight for customer node L[i]
    X = ElasticMatrix(fill(Inf, (I,J)))     # X[i,j]: insertion cost of customer node L[i] at best position in route R[j]
    P = ElasticMatrix(fill((0, 0), (I,J)))  # P[i,j]: best insertion postion of customer node L[i] in route R[j]
    # Step 2: Iterate until all open customer nodes have been inserted into the route
    for _ ∈ I
        # Step 2.1: Randomly select an open customer node and iterate through all possible insertion positions in each route
        z = f(s)
        i = sample(rng, I, Weights(W))
        c = L[i]
        for (j,r) ∈ pairs(R)
            d  = s.D[r.iᵈ]
            nˢ = isopt(r) ? C[r.iˢ] : D[r.iˢ]
            nᵉ = isopt(r) ? C[r.iᵉ] : D[r.iᵉ]
            nᵗ = d
            nʰ = nˢ
            while true
                # Step 2.1.1: Insert customer node c between tail node nᵗ and head node nʰ in route r
                insertnode!(c, nᵗ, nʰ, r, s)
                # Step 2.1.2: Compute the insertion cost
                z′ = f(s)
                Δ  = z′ - z
                # Step 2.1.3: Revise least insertion cost in route r and the corresponding best insertion position in route r
                if Δ < X[i,j] X[i,j], P[i,j] = Δ, (nᵗ.iⁿ, nʰ.iⁿ) end
                # Step 2.1.4: Remove customer node c from its position between tail node nᵗ and head node nʰ
                removenode!(c, nᵗ, nʰ, r, s)
                if isequal(nᵗ, nᵉ) break end
                nᵗ = nʰ
                nʰ = isequal(r.iᵉ, nᵗ.iⁿ) ? D[nᵗ.iʰ] : C[nᵗ.iʰ]
            end
        end
        # Step 2.2: Insert the customer node at its best position
        j  = argmin(X[i,:])
        r  = R[j]
        d  = s.D[r.iᵈ]
        v  = d.V[r.iᵛ]
        iᵗ = P[i,j][1]
        iʰ = P[i,j][2]
        nᵗ = iᵗ ≤ lastindex(D) ? D[iᵗ] : C[iᵗ]
        nʰ = iʰ ≤ lastindex(D) ? D[iʰ] : C[iʰ]
        insertnode!(c, nᵗ, nʰ, r, s)
        # Step 2.3: Revise vectors appropriately
        W[i] = 0
        X[i,:] .= Inf
        P[i,:] .= ((0, 0), )
        # Step 2.4: Update solution appropriately 
        if addroute(r, s)
            r = Route(v, d)
            push!(v.R, r)
            push!(R, r)
            append!(X, fill(Inf, (I,1)))
            append!(P, fill((0, 0), (I,1)))
        end
        if addvehicle(v, s)
            v = Vehicle(v, d)
            r = Route(v, d)
            push!(d.V, v)
            push!(v.R, r) 
            push!(R, r)
            append!(X, fill(Inf, (I,1)))
            append!(P, fill((0, 0), (I,1)))
        end
    end
    postinsert!(s)
    # Step 3: Return solution
    return s
end



"""
    greedy!(rng::AbstractRNG, s::Solution; mode::Symbol)

Returns solution `s` after iteratively inserting customer nodes with least 
insertion cost until all open customer nodes have been added to the solution. 
Available modes include `:pcs` (precise estimation of insertion cost) and 
`:ptb` (perturbed estimation of insertion cost).
"""
function greedy!(rng::AbstractRNG, s::Solution; mode::Symbol)
    if all(isclose, s.C) return s end
    # Step 1: Initialize
    preinsert!(s)
    D = s.D
    C = s.C
    φ = isequal(mode, :ptb)
    R = [r for d ∈ D for v ∈ d.V for r ∈ v.R if isactive(r)]
    L = [c for c ∈ C if isopen(c)]
    I = eachindex(L)
    J = eachindex(R)
    X = ElasticMatrix(fill(Inf, (I,J)))     # X[i,j]: insertion cost of customer node L[i] at best position in route R[j]
    P = ElasticMatrix(fill((0, 0), (I,J)))  # P[i,j]: best insertion postion of customer node L[i] in route R[j]
    ϕ = ones(Int, J)                        # ϕ[j]  : binary weight for route R[j]
    # Step 2: Iterate until all open customer nodes have been inserted into the route
    for _ ∈ I
        # Step 2.1: Iterate through all open customer nodes and every possible insertion position in each route
        z = f(s)
        for (i,c) ∈ pairs(L)
            if !isopen(c) continue end
            for (j,r) ∈ pairs(R)
                if iszero(ϕ[j]) continue end
                d  = s.D[r.iᵈ]
                nˢ = isopt(r) ? C[r.iˢ] : D[r.iˢ]
                nᵉ = isopt(r) ? C[r.iᵉ] : D[r.iᵉ]
                nᵗ = d
                nʰ = nˢ
                while true
                    # Step 2.1.1: Insert customer node c between tail node nᵗ and head node nʰ in route r
                    insertnode!(c, nᵗ, nʰ, r, s)
                    # Step 2.1.2: Compute the insertion cost
                    z′ = f(s) * (1 + φ * rand(rng, Uniform(-0.2, 0.2)))
                    Δ  = z′ - z
                    # Step 2.1.3: Revise least insertion cost in route r and the corresponding best insertion position in route r
                    if Δ < X[i,j] X[i,j], P[i,j] = Δ, (nᵗ.iⁿ, nʰ.iⁿ) end
                    # Step 2.1.4: Remove customer node c from its position between tail node nᵗ and head node nʰ
                    removenode!(c, nᵗ, nʰ, r, s)
                    if isequal(nᵗ, nᵉ) break end
                    nᵗ = nʰ
                    nʰ = isequal(r.iᵉ, nᵗ.iⁿ) ? D[nᵗ.iʰ] : C[nᵗ.iʰ]
                end
            end
        end
        # Step 2.2: Insert the customer node with least insertion cost at its best position
        i,j= Tuple(argmin(X))
        c  = L[i]
        r  = R[j]
        d  = s.D[r.iᵈ]
        v  = d.V[r.iᵛ]
        iᵗ = P[i,j][1]
        iʰ = P[i,j][2]
        nᵗ = iᵗ ≤ lastindex(D) ? D[iᵗ] : C[iᵗ]
        nʰ = iʰ ≤ lastindex(D) ? D[iʰ] : C[iʰ]
        insertnode!(c, nᵗ, nʰ, r, s)
        # Step 2.3: Revise vectors appropriately
        ϕ .= 0
        X[i,:] .= Inf
        P[i,:] .= ((0, 0), )
        for (j,r) ∈ pairs(R) 
            φʳ = isequal(r, c.r)
            φᵛ = isequal(r.iᵛ, v.iᵛ) && isless(c.r.tⁱ, r.tⁱ) && isequal(s.φ, true)
            φᵈ = isequal(r.iᵈ, d.iⁿ) && !hasslack(d)
            φˢ = φʳ || φᵛ || φᵈ
            if isequal(φˢ, false) continue end
            X[:,j] .= Inf
            P[:,j] .= ((0, 0), )
            ϕ[j] = 1
        end
        # Step 2.4: Update solution appropriately     
        if addroute(r, s)
            r = Route(v, d)
            push!(v.R, r)
            push!(R, r)
            append!(X, fill(Inf, (I,1)))
            append!(P, fill((0, 0), (I,1)))
            push!(ϕ, 1)
        end
        if addvehicle(v, s)
            v = Vehicle(v, d)
            r = Route(v, d)
            push!(d.V, v)
            push!(v.R, r) 
            push!(R, r)
            append!(X, fill(Inf, (I,1)))
            append!(P, fill((0, 0), (I,1)))
            push!(ϕ, 1)
        end
    end
    postinsert!(s)
    # Step 3: Return solution
    return s
end
"""
    precise!(rng::AbstractRNG, s::Solution)

Returns solution `s` after iteratively inserting customer nodes with 
least insertion cost until all open customer nodes have been added 
to the solution. 
Estimates insertion cost precisely.
"""
precise!(rng::AbstractRNG, s::Solution) = greedy!(rng, s; mode=:pcs)
"""
    precise!(rng::AbstractRNG, s::Solution)

Returns solution `s` after iteratively inserting customer nodes with 
least insertion cost until all open customer nodes have been added to 
the solution. 
Estimates insertion cost with a perturbration.
"""
perturb!(rng::AbstractRNG, s::Solution) = greedy!(rng, s; mode=:ptb)



"""
    regretk!(rng::AbstractRNG, s::Solution, k̅::Int)

Returns solution `s` after iteratively adding customer nodes with 
highest regret-k cost at its best position until all open customer 
nodes have been added to the solution.
Note, regretk mechanism breaks any ties by inserting the 
node with the lowest insertion cost.
"""
function regretk!(rng::AbstractRNG, s::Solution, k̅::Int)
    if all(isclose, s.C) return s end
    # Step 1: Initialize
    preinsert!(s)
    D = s.D
    C = s.C
    R = [r for d ∈ D for v ∈ d.V for r ∈ v.R if isactive(r)]
    L = [c for c ∈ C if isopen(c)]
    I = eachindex(L)
    J = eachindex(R)
    X = ElasticMatrix(fill(Inf, (I,J)))     # X[i,j]  : insertion cost of customer node L[i] at best position in route R[j]
    Y = ElasticArray(fill(Inf, (k̅,I,J)))    # Y[k,i,j]: insertion cost of customer node L[i] at kᵗʰ best position in route R[j]
    Z = fill(0., I)                         # Z[i]    : regret-k cost of customer node L[i]
    P = ElasticMatrix(fill((0, 0), (I,J)))  # P[i,j]  : best insertion postion of customer node L[i] in route R[j]
    ϕ = ones(Int, J)                        # ϕ[j]    : binary weight for route R[j]
    # Step 2: Iterate until all open customer nodes have been inserted into the route
    for _ ∈ I
        # Step 2.1: Iterate through all open customer nodes and every route
        z = f(s)
        for (i,c) ∈ pairs(L)
            if !isopen(c) continue end
            for (j,r) ∈ pairs(R)
                # Step 2.1.1: Iterate through all possible insertion position in route r
                if iszero(ϕ[j]) continue end
                d  = s.D[r.iᵈ]
                nˢ = isopt(r) ? C[r.iˢ] : D[r.iˢ]
                nᵉ = isopt(r) ? C[r.iᵉ] : D[r.iᵉ]
                nᵗ = d
                nʰ = nˢ
                while true
                    # Step 2.1.1.1: Insert customer node c between tail node nᵗ and head node nʰ in route r
                    insertnode!(c, nᵗ, nʰ, r, s)
                    # Step 2.1.1.2: Compute the insertion cost
                    z′ = f(s)
                    Δ  = z′ - z
                    # Step 2.1.1.3: Revise least insertion cost in route r and the corresponding best insertion position in route r
                    if Δ < X[i,j] X[i,j], P[i,j] = Δ, (nᵗ.iⁿ, nʰ.iⁿ) end
                    # Step 2.1.1.4: Revise k least insertion costs in route r
                    k̲ = 1
                    for k ∈ 1:k̅ Δ < Y[k,i,j] ? break : k̲ += 1 end
                    for k ∈ k̅:-1:k̲ Y[k,i,j] = isequal(k, k̲) ? Δ : Y[k-1,i,j] end
                    # Step 2.1.1.5: Remove customer node c from its position between tail node nᵗ and head node nʰ in route r
                    removenode!(c, nᵗ, nʰ, r, s)
                    if isequal(nᵗ, nᵉ) break end
                    nᵗ = nʰ
                    nʰ = isequal(r.iᵉ, nᵗ.iⁿ) ? D[nᵗ.iʰ] : C[nᵗ.iʰ]
                end
            end
            # Step 2.1.2: Compute regret cost for customer node c
            U = sort(vec(Y[:,i,:]))
            for k ∈ 1:k̅ Z[i] += U[k] - U[1] end
        end
        # Step 2.2: Insert customer node with highest regret cost in its best position (break ties by inserting the node with the lowest insertion cost)
        I̲  = findall(isequal.(Z, maximum(Z)))
        i,j= Tuple(argmin(X[I̲,:]))
        i  = I̲[i]
        c  = L[i]
        r  = R[j]
        d  = s.D[r.iᵈ]
        v  = d.V[r.iᵛ]
        iᵗ = P[i,j][1]
        iʰ = P[i,j][2]
        nᵗ = iᵗ ≤ lastindex(D) ? D[iᵗ] : C[iᵗ]
        nʰ = iʰ ≤ lastindex(D) ? D[iʰ] : C[iʰ]
        insertnode!(c, nᵗ, nʰ, r, s)
        # Step 2.3: Revise vectors appropriately
        ϕ .= 0
        Z .= 0.
        X[i,:] .= Inf
        Y[:,i,:] .= Inf
        P[i,:] .= ((0, 0), )
        for (j,r) ∈ pairs(R) 
            φʳ = isequal(r, c.r)
            φᵛ = isequal(r.iᵛ, v.iᵛ) && isless(c.r.tⁱ, r.tⁱ) && isequal(s.φ, true)
            φᵈ = isequal(r.iᵈ, d.iⁿ) && !hasslack(d)
            φˢ = φʳ || φᵛ || φᵈ
            if isequal(φˢ, false) continue end
            X[:,j] .= Inf
            Y[:,:,j] .= Inf
            P[:,j] .= ((0, 0), )
            ϕ[j] = 1  
        end
        # Step 2.4: Update solution appropriately     
        if addroute(r, s)
            r = Route(v, d)
            push!(v.R, r)
            push!(R, r)
            append!(X, fill(Inf, (I,1)))
            append!(Y, fill(Inf, (I,1,k̅)))
            append!(P, fill((0, 0), (I,1)))
            push!(ϕ, 1)
        end
        if addvehicle(v, s)
            v = Vehicle(v, d)
            r = Route(v, d)
            push!(d.V, v)
            push!(v.R, r) 
            push!(R, r)
            append!(X, fill(Inf, (I,1)))
            append!(Y, fill(Inf, (k̅,I,1)))
            append!(P, fill((0, 0), (I,1)))
            push!(ϕ, 1)
        end
    end
    postinsert!(s)
    # Step 3: Return solution
    return s
end
"""
    regret2!(rng::AbstractRNG, s::Solution)

Returns solution `s` after iteratively adding customer nodes with 
highest regret-2 cost at its best position until all open customer 
nodes have been added to the solution.
Note, regret2 mechanism breaks any ties by inserting the 
node with the lowest insertion cost.
"""
regret2!(rng::AbstractRNG, s::Solution) = regretk!(rng, s, 2)
"""
    regret3!(rng::AbstractRNG, s::Solution)

Returns solution `s` after iteratively adding customer nodes with 
highest regret-3 cost at its best position until all open customer 
nodes have been added to the solution.
Note, regret3 mechanism breaks any ties by inserting the 
node with the lowest insertion cost.
"""
regret3!(rng::AbstractRNG, s::Solution) = regretk!(rng, s, 3)