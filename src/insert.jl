"""
    insert!([rng], s::Solution, method::Symbol)

Return solution inserting open customer nodes to the solution `s` using the given `method`.

Available methods include,
- Best Insertion    : `:best!`
- Greedy Insertion  : `:greedy!`
- Regret Insertion  : `:regret2!`, `:regret3!`

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
    # Step 1: Initialize
    preinsertion(s)
    R = [r for d ∈ D for v ∈ d.V for r ∈ v.R]
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
        zᵒ = f(s; fixed=false)
        for (i,c) ∈ pairs(L)
            if !isopen(c) continue end
            for (j,r) ∈ pairs(R)
                if iszero(ϕ[j]) continue end
                d = s.D[r.iᵈ]
                v = d.V[r.iᵛ]
                m = length(v.R)
                nˢ = isopt(r) ? C[r.iˢ] : D[r.iˢ]
                nᵉ = isopt(r) ? C[r.iᵉ] : D[r.iᵉ]
                nᵗ = d
                nʰ = nˢ
                while true
                    # Step 2.1.1: Insert customer node c between tail node nᵗ and head node nʰ in route r
                    insertnode!(c, nᵗ, nʰ, r, s)
                    # Step 2.1.2: Compute the insertion cost
                    z⁺ = f(s; fixed=false) + c.q * (d.πᶠ/d.q + v.πᶠ/(m * v.q))
                    Δ  = z⁺ - zᵒ
                    # Step 2.1.3: Revise least insertion cost in route r and the corresponding best insertion position in route r
                    if Δ < x[i,j] x[i,j], p[i,j] = Δ, (nᵗ.iⁿ, nʰ.iⁿ) end
                    # Step 2.1.4: Remove customer node c from its position between tail node nᵗ and head node nʰ
                    removenode!(c, nᵗ, nʰ, r, s)
                    if isequal(nᵗ, nᵉ) break end
                    nᵗ = nʰ
                    nʰ = isequal(r.iᵉ, nᵗ.iⁿ) ? D[nᵗ.iʰ] : C[nᵗ.iʰ]
                end
            end
        end
        # Step 2.2: Randomly select a customer node to insert at its best position
        i = sample(rng, I, Weights(w))
        j = argmin(x[i,:])
        c = L[i]
        r = R[j]
        d = s.D[r.iᵈ]
        v = d.V[r.iᵛ]
        iᵗ = p[i,j][1]
        iʰ = p[i,j][2]
        nᵗ = iᵗ ≤ length(D) ? D[iᵗ] : C[iᵗ]
        nʰ = iʰ ≤ length(D) ? D[iʰ] : C[iʰ]
        insertnode!(c, nᵗ, nʰ, r, s)
        # Step 2.3: Revise vectors appropriately
        x[i,:] .= Inf
        p[i,:] .= ((0, 0), )
        w[i] = 0
        ϕ .= 0
        for (j,r) ∈ pairs(R) 
            if !isequal(r.iᵛ, v.iᵛ) continue end
            x[:,j] .= Inf
            p[:,j] .= ((0, 0), )
            ϕ[j] = 1  
        end
        if addroute(r, s)
            r = Route(v, d)
            push!(v.R, r) 
            push!(R, r)
            append!(x, fill(Inf, (I,1)))
            append!(p, fill((0, 0), (I,1)))
            push!(ϕ, 1)
        end
        if addvehicle(v, s)
            v = Vehicle(v, d)
            r = Route(v, d)
            push!(d.V, v)
            push!(v.R, r) 
            push!(R, r)
            append!(x, fill(Inf, (I,1)))
            append!(p, fill((0, 0), (I,1)))
            push!(ϕ, 1)
        end
    end
    postinsertion(s)
    # Step 3: Return solution
    return s
end

# Greedy insertion
# Iteratively insert customer nodes with least insertion cost until all open customer nodes have been added to the solution
function greedy!(rng::AbstractRNG, s::Solution)
    if all(isclose, s.C) return s end
    D = s.D
    C = s.C
    # Step 1: Initialize
    preinsertion(s)
    R = [r for d ∈ D for v ∈ d.V for r ∈ v.R]
    L = [c for c ∈ C if isopen(c)]
    I = eachindex(L)
    J = eachindex(R)
    x = ElasticMatrix(fill(Inf, (I,J)))     # x[i,j]: insertion cost of customer node L[i] at best position in route R[j]
    p = ElasticMatrix(fill((0, 0), (I,J)))  # p[i,j]: best insertion postion of customer node L[i] in route R[j]
    ϕ = ones(Int64, J)                      # ϕ[j]  : selection weight for route R[j]
    # Step 2: Iterate until all open customer nodes have been inserted into the route
    for _ ∈ I
        # Step 2.1: Iterate through all open customer nodes and every possible insertion position in each route
        zᵒ = f(s; fixed=false)
        for (i,c) ∈ pairs(L)
            if !isopen(c) continue end
            for (j,r) ∈ pairs(R)
                if iszero(ϕ[j]) continue end
                d = s.D[r.iᵈ]
                v = d.V[r.iᵛ]
                m = length(v.R)
                nˢ = isopt(r) ? C[r.iˢ] : D[r.iˢ]
                nᵉ = isopt(r) ? C[r.iᵉ] : D[r.iᵉ]
                nᵗ = d
                nʰ = nˢ
                while true
                    # Step 2.1.1: Insert customer node c between tail node nᵗ and head node nʰ in route r
                    insertnode!(c, nᵗ, nʰ, r, s)
                    # Step 2.1.2: Compute the insertion cost
                    z⁺ = f(s; fixed=false) + c.q * (d.πᶠ/d.q + v.πᶠ/(m * v.q))
                    Δ  = z⁺ - zᵒ
                    # Step 2.1.3: Revise least insertion cost in route r and the corresponding best insertion position in route r
                    if Δ < x[i,j] x[i,j], p[i,j] = Δ, (nᵗ.iⁿ, nʰ.iⁿ) end
                    # Step 2.1.4: Remove customer node c from its position between tail node nᵗ and head node nʰ
                    removenode!(c, nᵗ, nʰ, r, s)
                    if isequal(nᵗ, nᵉ) break end
                    nᵗ = nʰ
                    nʰ = isequal(r.iᵉ, nᵗ.iⁿ) ? D[nᵗ.iʰ] : C[nᵗ.iʰ]
                end
            end
        end
        # Step 2.2: Randomly select a customer node to insert at its best position        
        i,j = Tuple(argmin(x))
        c = L[i]
        r = R[j]
        d = s.D[r.iᵈ]
        v = d.V[r.iᵛ]
        iᵗ = p[i,j][1]
        iʰ = p[i,j][2]
        nᵗ = iᵗ ≤ length(D) ? D[iᵗ] : C[iᵗ]
        nʰ = iʰ ≤ length(D) ? D[iʰ] : C[iʰ]
        insertnode!(c, nᵗ, nʰ, r, s)
        # Step 2.3: Revise vectors appropriately
        x[i,:] .= Inf
        p[i,:] .= ((0, 0), )
        ϕ .= 0
        for (j,r) ∈ pairs(R) 
            if !isequal(r.iᵛ, v.iᵛ) continue end
            x[:,j] .= Inf
            p[:,j] .= ((0, 0), )
            ϕ[j] = 1  
        end
        if addroute(r, s)
            r = Route(v, d)
            push!(v.R, r) 
            push!(R, r)
            append!(x, fill(Inf, (I,1)))
            append!(p, fill((0, 0), (I,1)))
            push!(ϕ, 1)
        end
        if addvehicle(v, s)
            v = Vehicle(v, d)
            r = Route(v, d)
            push!(d.V, v)
            push!(v.R, r) 
            push!(R, r)
            append!(x, fill(Inf, (I,1)))
            append!(p, fill((0, 0), (I,1)))
            push!(ϕ, 1)
        end
    end
    postinsertion(s)
    # Step 3: Return solution
    return s
end

# Regret-N insertion
# Iteratively add customer nodes with highest regret cost at its best position until all open customer nodes have been added to the solution
function regretN!(rng::AbstractRNG, N::Int64, s::Solution)
    if all(isclose, s.C) return s end
    D = s.D
    C = s.C
    # Step 1: Initialize
    preinsertion(s)
    R = [r for d ∈ D for v ∈ d.V for r ∈ v.R]
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
        zᵒ = f(s; fixed=false)
        for (i,c) ∈ pairs(L)
            if !isopen(c) continue end
            for (j,r) ∈ pairs(R)
                # Step 2.1.1: Iterate through all possible insertion position in route r
                if iszero(ϕ[j]) continue end
                d = s.D[r.iᵈ]
                v = d.V[r.iᵛ]
                m = length(v.R)
                nˢ = isopt(r) ? C[r.iˢ] : D[r.iˢ]
                nᵉ = isopt(r) ? C[r.iᵉ] : D[r.iᵉ]
                nᵗ = d
                nʰ = nˢ
                while true
                    # Step 2.1.1.1: Insert customer node c between tail node nᵗ and head node nʰ in route r
                    insertnode!(c, nᵗ, nʰ, r, s)
                    # Step 2.1.1.2: Compute the insertion cost
                    z⁺ = f(s; fixed=false) + c.q * (d.πᶠ/d.q + v.πᶠ/(m * v.q))
                    Δ  = z⁺ - zᵒ
                    # Step 2.1.1.3: Revise least insertion cost in route r and the corresponding best insertion position in route r
                    if Δ < x[i,j] x[i,j], p[i,j] = Δ, (nᵗ.iⁿ, nʰ.iⁿ) end
                    # Step 2.1.1.4: Revise N least insertion costs
                    n̲ = 1
                    for n ∈ 1:N 
                        n̲ = n
                        if Δ < y[i,n] break end
                    end
                    for n ∈ N:-1:n̲ 
                        y[i,n] = isequal(n, n̲) ? Δ : y[i,n-1]::Float64
                        w[i,n] = isequal(n, n̲) ? r.iʳ : w[i,n-1]::Int64
                    end
                    # Step 2.1.1.5: Remove customer node c from its position between tail node nᵗ and head node nʰ in route r
                    removenode!(c, nᵗ, nʰ, r, s)
                    if isequal(nᵗ, nᵉ) break end
                    nᵗ = nʰ
                    nʰ = isequal(r.iᵉ, nᵗ.iⁿ) ? D[nᵗ.iʰ] : C[nᵗ.iʰ]
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
        d = s.D[r.iᵈ]
        v = d.V[r.iᵛ]
        iᵗ = p[i,j][1]
        iʰ = p[i,j][2]
        nᵗ = iᵗ ≤ length(D) ? D[iᵗ] : C[iᵗ]
        nʰ = iʰ ≤ length(D) ? D[iʰ] : C[iʰ]
        insertnode!(c, nᵗ, nʰ, r, s)
        # Step 2.3: Revise vectors appropriately
        x[i,:] .= Inf
        p[i,:] .= ((0, 0), )
        y[i,:] .= Inf
        w[i,:] .= 0
        z .= -Inf 
        for (i,c) ∈ pairs(L)
            for n ∈ 1:N
                if iszero(w[i,n]) break end
                k = findfirst(x -> x.iʳ == w[i,n], R)
                r = R[k]
                if isequal(r.iᵛ, v.iᵛ) y[i,n], w[i,n] = Inf, 0 end
            end
            ix = sortperm(y[i,:])
            y[i,:] .= y[i,ix]
            w[i,:] .= w[i,ix]
        end
        ϕ .= 0
        for (j,r) ∈ pairs(R) 
            if !isequal(r.iᵛ, v.iᵛ) continue end
            x[:,j] .= Inf
            p[:,j] .= ((0, 0), )
            ϕ[j] = 1  
        end
        if addroute(r, s)
            r = Route(v, d)
            push!(v.R, r) 
            push!(R, r)
            append!(x, fill(Inf, (I,1)))
            append!(p, fill((0, 0), (I,1)))
            push!(ϕ, 1)
        end
        if addvehicle(v, s)
            v = Vehicle(v, d)
            r = Route(v, d)
            push!(d.V, v)
            push!(v.R, r) 
            push!(R, r)
            append!(x, fill(Inf, (I,1)))
            append!(p, fill((0, 0), (I,1)))
            push!(ϕ, 1)
        end
    end
    postinsertion(s)
    # Step 3: Return solution
    return s
end
regret2!(rng::AbstractRNG, s::Solution) = regretN!(rng, Int64(2), s)
regret3!(rng::AbstractRNG, s::Solution) = regretN!(rng, Int64(3), s)