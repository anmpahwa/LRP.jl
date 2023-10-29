"""
    build(instance::String)
    
Returns a tuple of depot nodes, customer nodes, and arcs for the `instance`.
"""
function build(instance::String)
    # Depot nodes
    df = DataFrame(CSV.File(joinpath(dirname(@__DIR__), "instances/$instance/depot_nodes.csv")))
    D  = Vector{DepotNode}(undef, nrow(df))
    for k ∈ 1:nrow(df)
        iⁿ = df[k,1]
        jⁿ = df[k,2]
        x  = df[k,3]
        y  = df[k,4]
        q  = df[k,5]
        pˡ = df[k,6]
        pᵘ = df[k,7]
        tˢ = df[k,8]
        tᵉ = df[k,9]
        πᵒ = df[k,10]
        πᶠ = df[k,11]
        φ  = df[k,12]
        d  = DepotNode(iⁿ, jⁿ, x, y, q, pˡ, pᵘ, tˢ, tᵉ, πᵒ, πᶠ, φ, Vehicle[])
        D[iⁿ] = d
    end

    # Customer nodes
    df = DataFrame(CSV.File(joinpath(dirname(@__DIR__), "instances/$instance/customer_nodes.csv")))
    I  = (df[1,1]:df[nrow(df),1])
    C  = OffsetVector{CustomerNode}(undef, I)
    for k ∈ 1:nrow(df)
        iⁿ = df[k,1]
        iʳ = 0
        iᵛ = 0
        iᵈ = 0
        x  = df[k,2]
        y  = df[k,3]
        q  = df[k,4]
        τᶜ = df[k,5]
        tᵉ = df[k,6]
        tˡ = df[k,7]
        iᵗ = 0
        iʰ = 0
        tᵃ = 0.
        tᵈ = 0.
        c  = CustomerNode(iⁿ, iʳ, iᵛ, iᵈ, x, y, q, τᶜ, tᵉ, tˡ, iᵗ, iʰ, tᵃ, tᵈ, NullRoute)
        C[iⁿ] = c
    end

    # Arcs
    df = DataFrame(CSV.File(joinpath(dirname(@__DIR__), "instances/$instance/arcs.csv"), header=false))
    A  = Dict{Tuple{Int64,Int64},Arc}()
    N  = length(D)+length(C)
    for iᵗ ∈ 1:N
        for iʰ ∈ 1:N
            l = df[iᵗ,iʰ] 
            a = Arc(iᵗ, iʰ, l)
            A[(iᵗ,iʰ)] = a
        end
    end

    # Vehicles
    df = DataFrame(CSV.File(joinpath(dirname(@__DIR__), "instances/$instance/vehicles.csv")))
    for k ∈ 1:nrow(df)
        iᵛ = df[k,1]
        jᵛ = df[k,2]
        iᵈ = df[k,3]
        q  = df[k,4]
        l  = df[k,5]
        s  = df[k,6]
        τᶠ = df[k,7]
        τᵈ = df[k,8]
        τᶜ = df[k,9]
        τʷ = df[k,10]
        r̅  = df[k,11]
        πᵈ = df[k,12]
        πᵗ = df[k,13]
        πᶠ = df[k,14]
        d  = D[iᵈ]
        v  = Vehicle(iᵛ, jᵛ, iᵈ, q, l, s, τᶠ, τᵈ, τᶜ, τʷ, r̅, d.tˢ, d.tˢ, πᵈ, πᵗ, πᶠ, Route[])
        push!(d.V, v)
    end
    V  = [v for d ∈ D for v ∈ d.V]

    φᵈ = iszero(getproperty.(D, :tˢ)) && iszero(getproperty.(D, :tᵉ))
    φᶜ = iszero(getproperty.(C, :tᵉ)) && iszero(getproperty.(C, :tˡ))
    φᵛ = iszero(getproperty.(V, :τʷ)) && iszero(getproperty.(V, :πᵗ))
    global φᵀ = !(φᵈ && φᶜ && φᵛ)::Bool
    
    G  = (D, C, A)
    return G
end
