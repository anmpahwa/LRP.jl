"""
    insert!([rng], s::Solution, method::Symbol)

Return solution inserting open customer nodes to the solution `s` using the given `method`.

Available methods include,
- Best Insertion    : `:best!`
- Greedy Insertion  : `:greedy!`
- Regret Insertion  : `:regret₂insert!`, `:regret₃insert!`

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
insert!(rng::AbstractRNG, s::Solution, method::Symbol)::Solution = getfield(LRP, method)(rng, s)
insert!(s::Solution, method::Symbol) = insert!(Random.GLOBAL_RNG, s, method)

# Best insertion
# Iteratively insert randomly selected customer node at its best position until all open customer nodes have been added to the solution
function best!(rng::AbstractRNG, s::Solution)
    if all(isclose, s.C) return s end
    D = s.D
    C = s.C
    V = s.V
    R = s.R
    # Step 1: Initialize
    L = [c for c ∈ C if isopen(c)]
    I = eachindex(L)
    J = eachindex(R)
    x = ElasticMatrix(fill(Inf, (I,J)))     # x[i,j]: insertion cost of customer node L[i] at best position in route R[j]
    p = ElasticMatrix(fill((0, 0), (I,J)))  # p[i,j]: best insertion postion of customer node L[i] in route R[j]
    w = ones(Int64, I)                      # w[j]  : selection weight for customer node L[i]
    ϕ = ones(Int64, J)                      # ϕ[j]  : selection weight for route R[j]
    # Step 2: Iterate until all open customer nodes have been inserted into the route
    for _ ∈ I
        # Step 2.1: Iterate through all open customer nodes and every possible insertion position in each route
        zᵒ = f(s)
        for (i,c) ∈ pairs(L)
            if !isopen(c) continue end
            for (j,r) ∈ pairs(R)
                if iszero(ϕ[j]) continue end
                v  = V[r.o]
                d  = D[v.o]
                nₛ = isopt(r) ? C[r.s] : D[r.s]
                nₑ = isopt(r) ? C[r.e] : D[r.e]
                nₜ = d
                nₕ = nₛ
                while true
                    # Step 2.1.1: Insert customer node c between tail node nₜ and head node nₕ in route r
                    insertnode!(c, nₜ, nₕ, r, s)
                    # Step 2.1.2: Compute the insertion cost
                    z⁺ = f(s)
                    Δ  = (z⁺ + (d.πᶠ/d.q + v.πᶠ/v.q)) - (zᵒ - (d.πᶠ + v.πᶠ))
                    # Step 2.1.3: Revise least insertion cost in route r and the corresponding best insertion position in route r
                    if Δ < x[i,j] x[i,j], p[i,j] = Δ, (nₜ.i, nₕ.i) end
                    # Step 2.1.4: Remove customer node c from its position between tail node nₜ and head node nₕ
                    removenode!(c, nₜ, nₕ, r, s)
                    if isequal(nₜ, nₑ) break end
                    nₜ = nₕ
                    nₕ = isequal(r.e, nₜ.i) ? D[nₜ.h] : C[nₜ.h]
                end
            end
        end
        # Step 2.2: Randomly select a customer node to insert at its best position
        i = sample(rng, I, Weights(w))
        j = argmin(x[i,:])
        c = L[i]
        r = R[j]
        v = V[r.o]
        t, h = p[i,j]
        nₜ = t ≤ length(D) ? D[t] : C[t]
        nₕ = h ≤ length(D) ? D[h] : C[h]
        insertnode!(c, nₜ, nₕ, r, s)
        # Step 2.3: Revise vectors appropriately
        x[i,:] .= Inf
        p[i,:] .= ((0, 0), )
        w[i] = 0
        ϕ .= 0
        for (j,r) ∈ pairs(R) 
            if isequal(r.o, v.i) 
                x[:,j] .= Inf
                p[:,j] .= ((0, 0), )
                ϕ[j] = 1 
            end 
        end
        for (k,v) ∈ pairs(V)
            if any(!isopt, v.R) continue end
            d = D[v.o]
            r = Route(rand(rng, 1:M), v, d)
            push!(v.R, r) 
            push!(R, r)
            append!(x, fill(Inf, (I,1)))
            append!(p, fill((0, 0), (I,1)))
            push!(ϕ, 1)
        end
    end
    return s
end

# Greedy insertion
# Iteratively insert customer nodes with least insertion cost until all open customer nodes have been added to the solution
function greedy!(rng::AbstractRNG, s::Solution)
    if all(isclose, s.C) return s end
    D = s.D
    C = s.C
    V = s.V
    R = s.R
    # Step 1: Initialize
    L = [c for c ∈ C if isopen(c)]
    I = eachindex(L)
    J = eachindex(R)
    x = ElasticMatrix(fill(Inf, (I,J)))     # x[i,j]: insertion cost of customer node L[i] at best position in route R[j]
    p = ElasticMatrix(fill((0, 0), (I,J)))  # p[i,j]: best insertion postion of customer node L[i] in route R[j]
    ϕ = ones(Int64, J)                      # ϕ[j]  : selection weight for route R[j]
    # Step 2: Iterate until all open customer nodes have been inserted into the route
    for _ ∈ I
        # Step 2.1: Iterate through all open customer nodes and every possible insertion position in each route
        zᵒ = f(s)
        for (i,c) ∈ pairs(L)
            if !isopen(c) continue end
            for (j,r) ∈ pairs(R)
                if iszero(ϕ[j]) continue end
                v  = V[r.o]
                d  = D[v.o]
                nₛ = isopt(r) ? C[r.s] : D[r.s]
                nₑ = isopt(r) ? C[r.e] : D[r.e]
                nₜ = d
                nₕ = nₛ
                while true
                    # Step 2.1.1: Insert customer node c between tail node nₜ and head node nₕ in route r
                    insertnode!(c, nₜ, nₕ, r, s)
                    # Step 2.1.2: Compute the insertion cost
                    z⁺ = f(s)
                    Δ  = (z⁺ + (d.πᶠ/d.q + v.πᶠ/v.q)) - (zᵒ - (d.πᶠ + v.πᶠ))
                    # Step 2.1.3: Revise least insertion cost in route r and the corresponding best insertion position in route r
                    if Δ < x[i,j] x[i,j], p[i,j] = Δ, (nₜ.i, nₕ.i) end
                    # Step 2.1.4: Remove customer node c from its position between tail node nₜ and head node nₕ
                    removenode!(c, nₜ, nₕ, r, s)
                    if isequal(nₜ, nₑ) break end
                    nₜ = nₕ
                    nₕ = isequal(r.e, nₜ.i) ? D[nₜ.h] : C[nₜ.h]
                end
            end
        end
        # Step 2.2: Randomly select a customer node to insert at its best position        
        i,j = Tuple(argmin(x))
        c = L[i]
        r = R[j]
        v = V[r.o]
        t, h = p[i,j]
        nₜ = t ≤ length(D) ? D[t] : C[t]
        nₕ = h ≤ length(D) ? D[h] : C[h]
        insertnode!(c, nₜ, nₕ, r, s)
        # Step 2.3: Revise vectors appropriately
        x[i,:] .= Inf
        p[i,:] .= ((0, 0), )
        ϕ .= 0
        for (j,r) ∈ pairs(R) 
            if isequal(r.o, v.i) 
                x[:,j] .= Inf
                p[:,j] .= ((0, 0), )
                ϕ[j] = 1 
            end 
        end
        for (k,v) ∈ pairs(V)
            if any(!isopt, v.R) continue end
            d = D[v.o]
            r = Route(rand(rng, 1:M), v, d)
            push!(v.R, r) 
            push!(R, r)
            append!(x, fill(Inf, (I,1)))
            append!(p, fill((0, 0), (I,1)))
            push!(ϕ, 1)
        end
    end
    return s
end

# Regret-N insertion
# Iteratively add customer nodes with highest regret cost at its best position until all open customer nodes have been added to the solution
function regretₙinsert!(rng::AbstractRNG, N::Int64, s::Solution)
    if all(isclose, s.C) return s end
    D = s.D
    C = s.C
    V = s.V
    R = s.R
    # Step 1: Initialize
    L = [c for c ∈ C if isopen(c)]
    I = eachindex(L)
    J = eachindex(R)
    x = ElasticMatrix(fill(Inf, (I,J)))     # x[i,j]: insertion cost of customer node L[i] at best position in route R[j]
    p = ElasticMatrix(fill((0, 0), (I,J)))  # p[i,j]: best insertion postion of customer node L[i] in route R[j]
    y = fill(Inf, (I,N))                    # y[i,n]: insertion cost of customer node L[i] at nᵗʰ best position
    w = zeros(Int64, (I,N))                 # w[i,n]: route index of customer node L[j] at nᵗʰ best position
    z = fill(-Inf, I)                       # z[i]  : regret-N cost of customer node L[i]
    ϕ = ones(Int64, J)                      # ϕ[j]  : selection weight for route R[j]
    # Step 2: Iterate until all open customer nodes have been inserted into the route
    for _ ∈ I
        # Step 2.1: Iterate through all open customer nodes and every route
        zᵒ = f(s)
        for (i,c) ∈ pairs(L)
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
        c = L[i]
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
        for (i,c) ∈ pairs(L)
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
            d = D[v.o]
            r = Route(rand(rng, 1:M), v, d)
            push!(v.R, r) 
            push!(R, r)
            append!(x, fill(Inf, (I,1)))
            append!(p, fill((0, 0), (I,1)))
            push!(ϕ, 1)
        end
    end
    # Step 3: Return initial solution
    return s
end
regret₂insert!(rng::AbstractRNG, s::Solution) = regretₙinsert!(rng, Int64(2), s)
regret₃insert!(rng::AbstractRNG, s::Solution) = regretₙinsert!(rng, Int64(3), s)
