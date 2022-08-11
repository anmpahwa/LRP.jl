"""
    build(instance)

Builds graph as set of nodes and arcs for `instance`.
"""
function build(instance)
    f = joinpath(dirname(@__DIR__), "$instances/$instance.csv")
    csv = CSV.File(f, types=[Int64, Float64, Float64])
    df = DataFrame(csv)
    K = nrow(df)
    # Nodes
    N = Vector{Node}(undef, K)
    for k ∈ 1:K
        i = df[k,1]::Int64
        x = df[k,2]::Float64
        y = df[k,3]::Float64
        n = Node(i, x, y)
        N[i] = n
    end
    # Arcs
    A = Dict{Tuple{Int64,Int64},Arc}()
    for p ∈ 1:K
        i  = df[p,1]::Int64
        xᵢ = df[p,2]::Float64
        yᵢ = df[p,3]::Float64
        for q ∈ 1:K
            j  = df[q,1]::Int64
            xⱼ = df[q,2]::Float64
            yⱼ = df[q,3]::Float64
            c  = sqrt((xⱼ - xᵢ)^2 + (yⱼ - yᵢ)^2)
            a  = Arc(i, j, c)
            A[(i,j)] = a
        end
    end
    G = (N, A)
    return G
end