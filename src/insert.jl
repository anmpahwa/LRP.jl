"""
    insert!([rng], s::Solution, χₒ::ObjectiveFunctionParameters, method::Symbol)

Return solution inserting open customer nodes to the solution `s` using the given `method`.
`χₒ` includes the objective function parameters for objective function evaluation.

Available methods include,
- Best Insertion    : `best!`
- Greedy Insertion  : `greedy!`
- Regret Insertion  : `regret₂insert!`, `regret₃insert!`

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
insert!(rng::AbstractRNG, s::Solution, χₒ::ObjectiveFunctionParameters, method::Symbol)::Solution = getfield(LRP, method)(rng, s, χₒ)
insert!(s::Solution, χₒ::ObjectiveFunctionParameters, method::Symbol) = insert!(Random.GLOBAL_RNG, s, χₒ, method)

# Best insertion
# Iteratively insert randomly selected customer node at its best position until all open customer nodes have been added to the solution
function best!(rng::AbstractRNG, s::Solution, χₒ::ObjectiveFunctionParameters)
    D = s.D
    C = s.C
    V = s.V
    R = [r for v ∈ V for r ∈ v.R]
    L = [c for c ∈ C if isopen(c)]
    # Step 1: Initialize
    I = length(R)
    J = length(L)
    K = length(V)
    p = fill((0, 0), (I,J)) # p[i,j]: best insertion postion of customer node L[j] in route R[i]
    x = fill(Inf, (I,J))    # x[i,j]: insertion cost of customer node L[j] at best position in route R[i]
    w = ones(Int64, J)      # w[j]  : selection weight for customer node L[j]
    ϕ = ones(Int64, K)      # ϕ[k]  : selection weight for vehicle V[k]
    # Step 2: Iterate until all open customer nodes have been inserted into the route
    for _ ∈ 1:J
        # Step 2.1: Iterate through all open customer nodes and every possible insertion position in each route
        z = f(s, χₒ)
        for (j,c) ∈ pairs(L)
            if isclose(c) continue end
            for (i,r) ∈ pairs(R)
                if iszero(ϕ[r.o]) continue end
                v = V[r.o]
                d = D[v.o]
                nₛ = isclose(r) ? D[r.s] : C[r.s]
                nₑ = isclose(r) ? D[r.e] : C[r.e]
                nₜ = d
                nₕ = nₛ
                while true
                    # Step 2.1.1: Insert customer node c between tail node nₜ and head node nₕ in route r
                    insertnode!(c, nₜ, nₕ, r, s)
                    # Step 2.1.2: Compute the insertion cost
                    z⁺ = f(s, χₒ)
                    Δ  = z⁺ - z
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
        j = sample(rng, 1:J, Weights(w))
        i = argmin(x[:,j])
        r = R[i]
        c = L[j]
        v = V[r.o]
        t, h = p[i,j]
        nₜ = t ≤ length(D) ? D[t] : C[t]
        nₕ = h ≤ length(D) ? D[h] : C[h]
        insertnode!(c, nₜ, nₕ, r, s)
        # Step 2.3: Revise vectors appropriately
        for r ∈ R if isequal(r.o, v.i) p[i,:] .= ((0, 0), ) end end
        for r ∈ R if isequal(r.o, v.i) x[i,:] .= Inf end end
        w[j] = 0
        for (k,v) ∈ pairs(V) ϕ[k] = isequal(r.o, v.i) ? 1 : 0 end
    end
    return s
end

# Greedy insertion
# Iteratively insert customer nodes with least insertion cost until all open customer nodes have been added to the solution
function greedy!(rng::AbstractRNG, s::Solution, χₒ::ObjectiveFunctionParameters)
    D = s.D
    C = s.C
    V = s.V
    R = [r for v ∈ V for r ∈ v.R]
    L = [c for c ∈ C if isopen(c)]
    # Step 1: Initialize
    I = length(R)
    J = length(L)
    K = length(V)
    p = fill((0, 0), (I,J)) # p[i,j] : best insertion postion of customer node L[j] in route R[i]
    x = fill(Inf, (I,J))    # x[i,j]: insertion cost of customer node L[j] at best position route R[i]
    ϕ = ones(Int64, K)      # ϕ[k]  : selection weight for vehicle V[k]
    # Step 2: Iterate until all open customer nodes have been inserted into the route
    for _ ∈ 1:J
        # Step 2.1: Iterate through all open customer nodes and every possible insertion position in each route
        z = f(s, χₒ)
        for (j,c) ∈ pairs(L)
            if isclose(c) continue end
            for (i,r) ∈ pairs(R)
                if iszero(ϕ[r.o]) continue end
                v = V[r.o]
                d = D[v.o]
                nₛ = isclose(r) ? D[r.s] : C[r.s]
                nₑ = isclose(r) ? D[r.e] : C[r.e]
                nₜ = d
                nₕ = nₛ
                while true
                    # Step 2.1.1: Insert customer node c between tail node nₜ and head node nₕ in route r
                    insertnode!(c, nₜ, nₕ, r, s)
                    # Step 2.1.2: Compute the insertion cost
                    z⁺ = f(s, χₒ)
                    Δ  = z⁺ - z
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
        # Step 2.2: Insert node with least insertion cost at its best position
        i,j = Tuple(argmin(x))
        r = R[i]
        c = L[j]
        t, h = p[i,j]
        nₜ = t ≤ length(D) ? D[t] : C[t]
        nₕ = h ≤ length(D) ? D[h] : C[h]
        insertnode!(c, nₜ, nₕ, r, s)
        v = V[r.o]
        for r ∈ R if isequal(r.o, v.i) p[i,:] .= ((0, 0), ) end end
        for r ∈ R if isequal(r.o, v.i) x[i,:] .= Inf end end
        x[:,j] .= Inf
        for (k,v) ∈ pairs(V) ϕ[k] = isequal(r.o, v.i) ? 1 : 0 end
    end
    return s
end

# Regret-N insertion
# Iteratively add customer nodes with highest regret cost at its best position until all open customer nodes have been added to the solution
function regretₙinsert!(rng::AbstractRNG, N::Integer, s::Solution, χₒ::ObjectiveFunctionParameters)
    D = s.D
    C = s.C
    V = s.V
    R = [r for v ∈ V for r ∈ v.R]
    L = [c for c ∈ C if isopen(c)]
    # Step 1: Initialize
    I = length(R)
    J = length(L)
    K = length(V)
    xᵢ = fill(Inf, (I,J))               # xᵢ[i,j]: insertion cost of customer node C[j] at best position in route R[i]
    pᵢ = fill((0, 0), (I,J))            # pᵢ[i,j]: best insertion postion of customer node C[j] in route R[i]
    xₙ = fill(Inf, (N,J))               # xₙ[i,n]: insertion cost of customer node C[j] at nᵗʰ best position
    rₙ = zeros(Int64, (N,J))            # rₙ[i,j]: N best insertion route index of customer node C[j]
    xᵣ = fill(-Inf, J)                  # x[j]   : regret-N cost of customer node C[j]
    ϕ  = ones(Int64, K)                 # ϕ[k]   : selection weight for vehicle V[k]
    # Step 2: Iterate until all open customer nodes have been inserted into the route
    for _ ∈ 1:J
        # Step 2.1: Iterate through all open customer nodes and every route
        z = f(s, χₒ)
        for (j,c) ∈ pairs(L)
            if isclose(c) continue end
            for (i,r) ∈ pairs(R)
                # Step 2.1.1: Iterate through all possible insertion position in route r
                if iszero(ϕ[r.o]) continue end
                v = V[r.o]
                d = D[v.o]
                nₛ = isclose(r) ? D[r.s] : C[r.s]
                nₑ = isclose(r) ? D[r.e] : C[r.e]
                nₜ = d
                nₕ = nₛ
                while true
                    # Step 2.1.1.1: Insert customer node c between tail node nₜ and head node nₕ in route r, and compute the insertion cost
                    insertnode!(c, nₜ, nₕ, r, s)
                    # Step 2.1.1.2: Compute the insertion cost
                    z⁺ = f(s, χₒ)
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
        c = L[j]
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
        for j ∈ eachindex(L) 
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
    return s
end
regret₂insert!(rng::AbstractRNG, s::Solution, χₒ::ObjectiveFunctionParameters) = regretₙinsert!(rng, 2, s, χₒ)
regret₃insert!(rng::AbstractRNG, s::Solution, χₒ::ObjectiveFunctionParameters) = regretₙinsert!(rng, 3, s, χₒ)
