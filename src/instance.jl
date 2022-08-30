# Builds instance as a graph with set of depot nodes, customer nodes, arcs, and vehicles.
function build(instance)
    # Depot nodes
    file = joinpath(dirname(@__DIR__), "instances/$instance/depot_nodes.csv")
    csv = CSV.File(file, types=[Int64, Float64, Float64, Int64, Float64, Float64])
    df = DataFrame(csv)
    D = Vector{DepotNode}(undef, nrow(df))
    for k ∈ 1:nrow(df)
        i  = df[k,1]::Int64
        x  = df[k,2]::Float64
        y  = df[k,3]::Float64
        q  = df[k,4]::Int64 
        πᵒ = df[k,5]::Float64
        πᶠ = df[k,6]::Float64
        d  = DepotNode(i, x, y, q, πᵒ, πᶠ, Vehicle[])
        D[i] = d
    end
    # Customer nodes
    file = joinpath(dirname(@__DIR__), "instances/$instance/customer_nodes.csv")
    csv = CSV.File(file, types=[Int64, Float64, Float64, Int64])
    df = DataFrame(csv)
    ix = (df[1,1]:df[nrow(df),1])::UnitRange{Int64}
    C = OffsetVector{CustomerNode}(undef, ix)
    for k ∈ 1:nrow(df)
        i = df[k,1]::Int64
        x = df[k,2]::Float64
        y = df[k,3]::Float64
        q = df[k,4]::Int64
        iₜ = 0
        iₕ = 0
        c = CustomerNode(i, x, y, q, iₜ, iₕ, NullRoute)
        C[i] = c
    end
    # Arcs
    A = Dict{Tuple{Int64,Int64},Arc}()
    N = length(D)+length(C)
    for iₜ ∈ 1:N
        nₜ = iₜ ≤ length(D) ? D[iₜ] : C[iₜ]
        xₜ = nₜ.x
        yₜ = nₜ.y
        for iₕ ∈ 1:N
            nₕ = iₕ ≤ length(D) ? D[iₕ] : C[iₕ]
            xₕ = nₕ.x
            yₕ = nₕ.y
            l  = sqrt((xₕ - xₜ)^2 + (yₕ - yₜ)^2)
            a  = Arc(iₜ, iₕ, l)
            A[(iₜ,iₕ)] = a
        end
    end
    # Vehicles
    file = joinpath(dirname(@__DIR__), "instances/$instance/vehicles.csv")
    csv = CSV.File(file, types=[Int64, Int64, Int64, Float64, Float64])
    df = DataFrame(csv)
    V = Vector{Vehicle}(undef, nrow(df))
    for k ∈ 1:nrow(df)
        i  = df[k,1]::Int64
        o  = df[k,2]::Int64
        q  = df[k,3]::Int64
        πᵒ = df[k,4]::Float64
        πᶠ = df[k,5]::Float64
        v  = Vehicle(i, o, q, πᵒ, πᶠ, Route[])
        V[k] = v
        push!(D[o].V, v)
    end
    G = (D, C, A, V)
    return G
end