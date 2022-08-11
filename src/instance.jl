"""
    build(instance)

Builds graph as set of depot nodes, customer nodes, and arcs for `instance`.
"""
function build(instance)
    # Vehicle
    file = joinpath(dirname(@__DIR__), "instances/$instance/vehicles.csv")
    csv = CSV.File(file, types=[Int64, Int64, Int64, Float64, Float64, Float64, Float64])
    df = DataFrame(csv)
    V = Vector{Vehicle}(undef, nrow(df))
    for k ∈ 1:nrow(df)
        i  = df[k,1]::Int64
        o  = df[k,2]::Int64
        q  = df[k,3]::Int64
        πᵐ = df[k,4]::Float64
        πʷ = df[k,5]::Float64
        πᶠ = df[k,6]::Float64
        πᵛ = df[k,7]::Float64
        v  = Vehicle(i, o, q, πᵐ, πʷ, πᶠ, πᵛ, Route[])
        V[k] = v
    end
    # Depot nodes
    file = joinpath(dirname(@__DIR__), "instances/$instance/depot_nodes.csv")
    csv = CSV.File(file, types=[Int64, Float64, Float64, Int64, Float64])
    df = DataFrame(csv)
    D = Vector{DepotNode}(undef, nrow(df))
    for k ∈ 1:nrow(df)
        i  = df[k,1]::Int64
        x  = df[k,2]::Float64
        y  = df[k,3]::Float64
        q  = df[k,4]::Int64 
        πᵈ = df[k,5]::Float64
        d  = DepotNode(i, x, y, q, πᵈ, [v for v ∈ V if isequal(v.o, i)])
        D[i] = d
    end
    # Customer nodes
    file = joinpath(dirname(@__DIR__), "instances/$instance/customer_nodes.csv")
    csv = CSV.File(file, types=[Int64, Float64, Float64, Int64])
    df = DataFrame(csv)
    ix = (df[begin,1]:df[end,1])::UnitRange{Int64}
    rₒ = Route()
    C = OffsetVector{CustomerNode}(undef, ix)
    for k ∈ 1:nrow(df)
        i = df[k,1]::Int64
        x = df[k,2]::Float64
        y = df[k,3]::Float64
        q = df[k,4]::Int64
        t = 0
        h = 0
        c = CustomerNode(i, x, y, q, t, h, rₒ)
        C[i] = c
    end
    # Arcs
    file = joinpath(dirname(@__DIR__), "instances/$instance/arcs.csv")
    csv = CSV.File(file, types=[Int64, Int64, Float64, Float64, Float64])
    df = DataFrame(csv)
    A = Dict{Tuple{Int64,Int64},Arc}()
    for k ∈ 1:nrow(df)
        i = df[k,1]::Int64
        j = df[k,2]::Int64
        l = df[k,3]::Float64
        t = df[k,4]::Float64
        f = df[k,5]::Float64
        a = Arc(i, j, l, t, f)
        A[(i,j)] = a       
    end
    G = (D, C, A, V)
    return G
end