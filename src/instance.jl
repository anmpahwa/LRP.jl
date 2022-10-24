# Builds instance as a graph with set of depot nodes, customer nodes, and arcs.
function build(instance; root=joinpath(dirname(@__DIR__), "instances"))
    # Depot nodes
    file = joinpath(root, "$instance/depot_nodes.csv")
    csv = CSV.File(file, types=[Int64, Int64, Float64, Float64, Int64, Float64, Float64, Float64, Float64, Float64, Float64])
    df = DataFrame(csv)
    D = Vector{DepotNode}(undef, nrow(df))
    for k ∈ 1:nrow(df)
        iⁿ = df[k,1]::Int64
        jⁿ = df[k,2]::Int64
        x  = df[k,3]::Float64
        y  = df[k,4]::Float64
        q  = df[k,5]::Int64
        pˡ = df[k,6]::Float64
        pᵘ = df[k,7]::Float64
        tˢ = df[k,8]::Float64
        tᵉ = df[k,9]::Float64
        πᵒ = df[k,10]::Float64
        πᶠ = df[k,11]::Float64
        d  = DepotNode(iⁿ, jⁿ, x, y, q, pˡ, pᵘ, tˢ, tᵉ, πᵒ, πᶠ, Vehicle[])
        D[iⁿ] = d
    end
    φᴱ = Int64(!isone(length(unique(getproperty.(D, :jⁿ)))))::Int64

    # Customer nodes
    file = joinpath(root, "$instance/customer_nodes.csv")
    csv = CSV.File(file, types=[Int64, Float64, Float64, Int64, Float64, Float64])
    df = DataFrame(csv)
    ix = (df[1,1]:df[nrow(df),1])::UnitRange{Int64}
    C = OffsetVector{CustomerNode}(undef, ix)
    for k ∈ 1:nrow(df)
        iⁿ = df[k,1]::Int64
        x  = df[k,2]::Float64
        y  = df[k,3]::Float64
        q  = df[k,4]::Int64
        tᵉ = df[k,5]::Float64
        tˡ = df[k,6]::Float64
        iᵗ = 0
        iʰ = 0
        tᵃ = 0
        tᵈ = 0
        c  = CustomerNode(iⁿ, x, y, q, tᵉ, tˡ, iᵗ, iʰ, tᵃ, tᵈ, NullRoute)
        C[iⁿ] = c
    end
    # Arcs
    A = Dict{Tuple{Int64,Int64},Arc}()
    N = length(D)+length(C)
    for iᵗ ∈ 1:N
        nᵗ = iᵗ ≤ length(D) ? D[iᵗ] : C[iᵗ]
        xᵗ = nᵗ.x
        yᵗ = nᵗ.y
        for iʰ ∈ 1:N
            nʰ = iʰ ≤ length(D) ? D[iʰ] : C[iʰ]
            xʰ = nʰ.x
            yʰ = nʰ.y
            l  = sqrt((xʰ - xᵗ)^2 + (yʰ - yᵗ)^2)
            a  = Arc(iᵗ, iʰ, l)
            A[(iᵗ,iʰ)] = a
        end
    end
    # Vehicles
    file = joinpath(root, "$instance/vehicles.csv")
    csv = CSV.File(file, types=[Int64, Int64, Int64, Int64, Int64, Int64, Float64, Float64, Float64, Int64, Int64, Float64, Float64])
    df = DataFrame(csv)
    for k ∈ 1:nrow(df)
        iᵛ = df[k,1]::Int64
        jᵛ = df[k,2]::Int64
        iᵈ = df[k,3]::Int64
        q  = df[k,4]::Int64
        l  = df[k,5]::Int64
        s  = df[k,6]::Int64
        τᶠ = df[k,7]::Float64
        τᵈ = df[k,8]::Float64
        τᶜ = df[k,9]::Float64
        τʷ = df[k,10]::Int64
        r̅  = df[k,11]::Int64
        πᵒ = df[k,12]::Float64
        πᶠ = df[k,13]::Float64
        d  = D[iᵈ]
        v  = Vehicle(iᵛ, jᵛ, iᵈ, q, l, s, τᶠ, τᵈ, τᶜ, τʷ, r̅, d.tˢ, d.tˢ, πᵒ, πᶠ, Route[])
        push!(d.V, v)
    end
    V  = [v for d ∈ D for v ∈ d.V]
    φᵀ = Int64(!(iszero(getproperty.(D, :tˢ)) && iszero(getproperty.(D, :tᵉ)) && iszero(getproperty.(C, :tᵉ)) && iszero(getproperty.(C, :tˡ)) && iszero(getproperty.(V, :τʷ))))::Int64
    G  = (D, C, A, φᴱ, φᵀ)
    return G
end