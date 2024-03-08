"""
    build(instance::String; dir=joinpath(dirname(@__DIR__), "instances"))
    
Returns a tuple of depot nodes, customer nodes, and arcs for the `instance`.

Note, `dir` locates the the folder containing instance files as sub-folders,
as follows,

    <dir>
    |-<instance>
        |-arcs.csv
        |-customer_nodes.csv
        |-depot_nodes.csv
        |-vehicles.csv
"""
function build(instance::String; dir=joinpath(dirname(@__DIR__), "instances"))
    # Depot nodes
    df = DataFrame(CSV.File(joinpath(dir, "$instance/depot_nodes.csv")))
    D  = Vector{DepotNode}(undef, nrow(df))
    for k ∈ 1:nrow(df)
        iⁿ = df[k,1]
        x  = df[k,2]
        y  = df[k,3]
        qᵈ = df[k,4]
        tˢ = df[k,5]
        tᵉ = df[k,6]
        πᵒ = df[k,7]
        πᶠ = df[k,8]
        φ  = df[k,9]
        τ  = Inf
        n  = 0
        q  = 0.
        l  = 0.
        d  = DepotNode(iⁿ, x, y, qᵈ, tˢ, tᵉ, πᵒ, πᶠ, φ, τ, n, q, l, Vehicle[])
        D[iⁿ] = d
    end
    # Vehicles
    df = DataFrame(CSV.File(joinpath(dir, "$instance/vehicles.csv")))
    for k ∈ 1:nrow(df)
        d  = D[df[k,3]]
        iᵛ = df[k,1]
        jᵛ = df[k,2]
        iᵈ = df[k,3]
        qᵛ = df[k,4]
        lᵛ = df[k,5]
        sᵛ = df[k,6]
        τᶠ = df[k,7]
        τᵈ = df[k,8]
        τᶜ = df[k,9]
        τʷ = df[k,10]
        r̅  = df[k,11]
        πᵈ = df[k,12]
        πᵗ = df[k,13]
        πᶠ = df[k,14]
        tˢ = d.tˢ
        tᵉ = d.tˢ
        τ  = Inf
        n  = 0
        q  = 0.
        l  = 0.
        v  = Vehicle(iᵛ, jᵛ, iᵈ, qᵛ, lᵛ, sᵛ, τᶠ, τᵈ, τᶜ, τʷ, r̅, πᵈ, πᵗ, πᶠ, tˢ, tᵉ, τ, n, q, l, Route[])
        push!(d.V, v)
    end
    # Customer nodes
    df = DataFrame(CSV.File(joinpath(dir, "$instance/customer_nodes.csv")))
    I  = (df[1,1]:df[nrow(df),1])::UnitRange{Int64}
    C  = OffsetVector{CustomerNode}(undef, I)
    for k ∈ 1:nrow(df)
        iⁿ = df[k,1]
        x  = df[k,2]
        y  = df[k,3]
        qᶜ = df[k,4]
        τᶜ = df[k,5]
        tʳ = df[k,6]
        tᵉ = df[k,7]
        tˡ = df[k,8]
        iᵗ = 0
        iʰ = 0
        tˢ = 0.
        tᵃ = 0.
        tᵈ = 0.
        r  = NullRoute
        c  = CustomerNode(iⁿ, x, y, qᶜ, τᶜ, tʳ, tᵉ, tˡ, iᵗ, iʰ, tˢ, tᵃ, tᵈ, r)
        C[iⁿ] = c
    end
    # Arcs
    df = DataFrame(CSV.File(joinpath(dir, "$instance/arcs.csv"), header=false))
    A  = Dict{Tuple{Int,Int},Arc}()
    n  = lastindex(C)
    for iᵗ ∈ 1:n
        for iʰ ∈ 1:n
            l = df[iᵗ,iʰ] 
            a = Arc(iᵗ, iʰ, l)
            A[(iᵗ,iʰ)] = a
        end
    end
    G  = (D, C, A)
    return G
end



"""
    cluster(rng::AbstractRNG, k::Int, instance::String; dir=joinpath(dirname(@__DIR__), "instances"))

Returns `Solution` created using `k`-means clustering algorithm.
Here, each cluster is assigned to the nearest depot node. Each 
customer in a cluster is then best inserted into the assigned 
depot node until this depot is capacitated. Any remaining customer 
nodes are finally inserted using best insertion method over the 
entire solution.

Note, `dir` locates the the folder containing instance files as sub-folders,
as follows,

    <dir>
    |-<instance>
        |-arcs.csv
        |-customer_nodes.csv
        |-depot_nodes.csv
        |-vehicles.csv
"""
function cluster(rng::AbstractRNG, k::Int, instance::String; dir=joinpath(dirname(@__DIR__), "instances"))
    # Step 1: Initialize
    G = build(instance; dir=dir)
    s = Solution(G...)
    preinitialize!(s)
    D = s.D
    C = s.C
    # Step 2: Clustering
    N = zeros(4, eachindex(C))
    for (iⁿ,c) ∈ pairs(C) N[:,iⁿ] = [c.x, c.y, c.tᵉ, c.tˡ] end
    K = kmeans(N.parent, k; rng=rng)
    A = OffsetVector(K.assignments, eachindex(C))
    M = K.centers
    # Step 3: Assignment
    Y = zeros(Int, nclusters(K))  
    for k ∈ 1:nclusters(K)
        xᵏ = M[1,k]
        yᵏ = M[2,k]
        lᵏ = Inf
        for j ∈ eachindex(D)
            xʲ = D[j].x
            yʲ = D[j].y
            lᵏʲ  = sqrt((xᵏ-xʲ)^2 + (yᵏ-yʲ)^2)
            Y[k] = lᵏʲ < lᵏ ? j : Y[k] 
            lᵏ   = lᵏʲ < lᵏ ? lᵏʲ : lᵏ
        end
    end
    # Step 4: Add customers from each cluster to the assigned depot
    for k ∈ 1:nclusters(K)
        d = D[Y[k]]
        R = [r for v ∈ d.V for r ∈ v.R]
        L = filter(c -> isequal(A[c.iⁿ], k), C)
        if isempty(L) continue end
        I = eachindex(L)
        J = eachindex(R)
        W = ones(Int, I)                        # W[j]  : selection weight for customer node L[i]
        X = ElasticMatrix(fill(Inf, (I,J)))     # X[i,j]: insertion cost of customer node L[i] at best position in route R[j]
        P = ElasticMatrix(fill((0, 0), (I,J)))  # P[i,j]: best insertion postion of customer node L[i] in route R[j]
        # Step 4.1: Iterate through all open customer nodes until the depot is capacitated
        for _ ∈ I
            if !hasslack(d) break end
            # Step 4.1.1: Randomly select an open customer nodes and iterate through all possible insertion positions in each route
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
                    # Step 4.1.1.1: Insert customer node c between tail node nᵗ and head node nʰ in route r
                    insertnode!(c, nᵗ, nʰ, r, s)
                    # Step 4.1.1.2: Compute the insertion cost
                    z′ = f(s)
                    Δ  = z′ - z
                    # Step 4.1.1.3: Revise least insertion cost in route r and the corresponding best insertion position in route r
                    if Δ < X[i,j] X[i,j], P[i,j] = Δ, (nᵗ.iⁿ, nʰ.iⁿ) end
                    # Step 4.1.1.4: Remove customer node c from its position between tail node nᵗ and head node nʰ
                    removenode!(c, nᵗ, nʰ, r, s)
                    if isequal(nᵗ, nᵉ) break end
                    nᵗ = nʰ
                    nʰ = isequal(r.iᵉ, nᵗ.iⁿ) ? D[nᵗ.iʰ] : C[nᵗ.iʰ]
                end
            end
            # Step 4.1.2: Randomly select a customer node to insert at its best position
            j  = argmin(X[i,:])
            r  = R[j]
            d  = s.D[r.iᵈ]
            v  = d.V[r.iᵛ]
            iᵗ = P[i,j][1]
            iʰ = P[i,j][2]
            nᵗ = iᵗ ≤ lastindex(D) ? D[iᵗ] : C[iᵗ]
            nʰ = iʰ ≤ lastindex(D) ? D[iʰ] : C[iʰ]
            insertnode!(c, nᵗ, nʰ, r, s)
            # Step 4.1.3: Revise vectors appropriately
            W[i] = 0
            # Step 4.1.4: Update solution appropriately     
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
    end
    # Step 5: Insert any remaining open customer nodes 
    best!(rng, s)
    postinitialize!(s)
    # Step 6: Return initial solution
    return s
end



"""
    initialize([rng::AbstractRNG], instance::String; method=:local, dir=joinpath(dirname(@__DIR__), "instances"))

Returns initial VRP `Solution` developed using iterated clustering method. 
If the `method` is set to `:local` search, the number of clusters are 
increased iteratively for at most as many iterations as the number of 
depot nodes. Else if the `method` is set to `:global` search, the number 
of clusters are increased iteratively for at least as many iterations as 
the number of depot nodes and at most the number of customer nodes until a 
feasible solution is found. Finally, the solution with the least objective 
function value is returned as the initial solution.

Note, `dir` locates the the folder containing instance files as sub-folders,
as follows,

    <dir>
    |-<instance>
        |-arcs.csv
        |-customer_nodes.csv
        |-depot_nodes.csv
        |-vehicles.csv

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
function initialize(rng::AbstractRNG, instance::String; method=:local, dir=joinpath(dirname(@__DIR__), "instances"))
    # Step 1. Initialize
    G = build(instance; dir=dir)
    s = Solution(G...)
    z = Inf
    # Step 2. Iteratively increase the number of clusters
    k = 0
    k̲ = length(s.D)
    k̅ = length(s.C)
    ϕ = isequal(method, :local)
    while k < k̅
        k += 1
        s′ = cluster(rng, k, instance; dir=dir)
        z′ = f(s′)
        # Step 2.1. Update solution
        if z′ < z
            z = z′ 
            s = deepcopy(s′)
        end
        # Step 2.2. Check for break conditions
        ϕ = ϕ || isfeasible(s′)
        k < k̲ ? continue : (ϕ ? break : continue)
    end
    # Step 3. Return solution
    return s
end
initialize(instance::String; method=:local, dir=joinpath(dirname(@__DIR__), "instances")) = initialize(Random.GLOBAL_RNG, instance; method=:method, dir=dir)