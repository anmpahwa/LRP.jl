# Builds instance as a graph with set of depot nodes, customer nodes, and arcs.
function build(instance)
    # Depot nodes
    file = joinpath(dirname(@__DIR__), "instances/$instance/depot_nodes.csv")
    csv = CSV.File(file, types=[Int64, Float64, Float64, Int64, Float64, Float64])
    df = DataFrame(csv)
    D = Vector{DepotNode}(undef, nrow(df))
    for k ∈ 1:nrow(df)
        iⁿ = df[k,1]::Int64
        x  = df[k,2]::Float64
        y  = df[k,3]::Float64
        q  = df[k,4]::Int64 
        πᵒ = df[k,5]::Float64
        πᶠ = df[k,6]::Float64
        d  = DepotNode(iⁿ, x, y, q, πᵒ, πᶠ, Vehicle[])
        D[iⁿ] = d
    end
    # Customer nodes
    file = joinpath(dirname(@__DIR__), "instances/$instance/customer_nodes.csv")
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
        tᵃ = Inf
        tᵈ = Inf
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
    file = joinpath(dirname(@__DIR__), "instances/$instance/vehicles.csv")
    csv = CSV.File(file, types=[Int64, Int64, Int64, Int64, Float64, Float64, Float64, Float64])
    df = DataFrame(csv)
    for k ∈ 1:nrow(df)
        iᵛ = df[k,1]::Int64
        iᵈ = df[k,2]::Int64
        q  = df[k,3]::Int64
        s  = df[k,4]::Int64
        τᵈ = df[k,5]::Float64
        τᶜ = df[k,6]::Float64
        πᵒ = df[k,7]::Float64
        πᶠ = df[k,8]::Float64
        v  = Vehicle(iᵛ, iᵈ, q, s, τᵈ, τᶜ, πᵒ, πᶠ, Route[])
        d  = D[iᵈ]
        push!(d.V, v)
    end
    G = (D, C, A)
    return G
end