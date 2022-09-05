@doc """
    Arc(iₜ::Int64, iₕ::Int64, l::Float64)

An `Arc` is a connection between tail   
node with index `iᵗ` and head node with 
index `iʰ` with length `l`.
"""
struct Arc
    iᵗ::Int64                                                                       # Tail node index
    iʰ::Int64                                                                       # Head node index
    l::Float64                                                                      # Arc length
end

@doc """
    Route(iʳ::Int64, iᵛ::Int64, iᵈ::Int64, iˢ::Int64, iᵉ::Int64, n::Int64, q::Int64, l::Float64)

A `Route` is a connection between nodes, with index `iʳ`, vehicle index `iᵛ`, depot
node index `iᵈ`, start node index `iˢ`, end node index `iᵉ`, number of customers `n`
, load `q`, and length `l`.
"""
mutable struct Route
    iʳ::Int64                                                                       # Route index
    iᵛ::Int64                                                                       # Vehicle index
    iᵈ::Int64                                                                       # Depot node index
    iˢ::Int64                                                                       # Route start node index
    iᵉ::Int64                                                                       # Route end node index
    n::Int64                                                                        # Route size (number of customers)
    q::Int64                                                                        # Route load
    l::Float64                                                                      # Route length
end
    
@doc """
    Vehicle(iᵛ::Int64, o::Int64, q::Int64, πᵒ::Float64, πᶠ::Float64, R::Vector{Route})

A `Vehicle` is a mode of delivery with index `iᵛ` from depot node index `iᵈ`, 
capacity `q`, operational cost `πᵒ`, fixed cost `πᶠ`, and set of routes `R`.
"""
struct Vehicle
    iᵛ::Int64                                                                       # Vehicle index
    iᵈ::Int64                                                                       # Depot node index
    q::Int64                                                                        # Vehicle capacity
    πᵒ::Float64                                                                     # Operational cost
    πᶠ::Float64                                                                     # Fixed cost
    R::Vector{Route}                                                                # Vector of vehicle routes
end

@doc """
    Node

A `Node` is a point on the graph.
"""
abstract type Node end

@doc """
    DepotNode(iⁿ::Int64, x::Float64, y::Float64, q::Float64, πᵒ::Float64, πᶠ::Float64, V::Vector{Vehicle})

A `DepotNode` is a source point on the graph at `(x,y)` with index `iⁿ`, capacity 
`q`, operational cost `πᵒ`, fixed cost `πᶠ`, and fleet of vehicles `V`.
"""
struct DepotNode <: Node
    iⁿ::Int64                                                                        # Depot node index
    x::Float64                                                                      # Depot node location on the x-axis
    y::Float64                                                                      # Depot node location in the y-axis
    q::Int64                                                                        # Depot capacity
    πᵒ::Float64                                                                     # Operational cost
    πᶠ::Float64                                                                     # Fixed cost
    V::Vector{Vehicle}                                                              # Vector of depot vehicles
end

@doc """
    CustomerNode(iⁿ::Int64, x::Float64, y::Float64, q::Float64, iᵗ::Int64, iʰ::Int64, r::Route)

A `CustomerNode` is a sink point on the graph at `(x,y)` with index `iⁿ`, demand `q`
, tail node index `iᵗ`, head node index `iʰ`, on route `r`.
"""
mutable struct CustomerNode <: Node
    iⁿ::Int64                                                                       # Customer node index
    x::Float64                                                                      # Customer node location on the x-axis
    y::Float64                                                                      # Customer node location in the y-axis
    q::Int64                                                                        # Customer demand
    iᵗ::Int64                                                                       # Tail (predecessor) node index
    iʰ::Int64                                                                       # Head (successor) node index
    r::Route                                                                        # Route visiting the customer
end

@doc """
    Solution(D::Vector{DepotNode}, C::Vector{CustomerNode}, A::Dict{Tuple{Int64,Int64}, Arc})

A Solution is a graph with depot nodes `D`, customer nodes `C`, and arcs `A`.
"""
struct Solution
    D::Vector{DepotNode}                                                            # Vector of depot nodes
    C::OffsetVector{CustomerNode, Vector{CustomerNode}}                             # Vector of customer nodes
    A::Dict{Tuple{Int64,Int64}, Arc}                                                # Set of arcs
end

# is operational
isopt(r::Route) = (r.n ≥ 1)                                                         # A route is defined operational if it serves at least one customer
isopt(v::Vehicle) = any(isopt, v.R)                                                 # A vehicle is defined operational if any of its routes is operational
isopt(d::DepotNode) = any(isopt, d.V)                                               # A depot is defined operational if any of its vehicles is operational
isopen(c::CustomerNode) = isequal(c.r, NullRoute)                                   # A customer is defined open if it is not being served by any vehicle-route

# is close
isclose(d::DepotNode) = !isopt(d)                                                   # A depot node is defined closed if it is non-operational
isclose(c::CustomerNode) = !isopen(c)                                               # A customer node is defined closed it is being served by any vehicle-route

# is equal
Base.isequal(p::Route, q::Route) = isequal(p.iʳ, q.iʳ)
Base.isequal(p::Vehicle, q::Vehicle) = isequal(p.iᵛ, q.iᵛ)
Base.isequal(p::Node, q::Node) = isequal(p.iⁿ, q.iⁿ)

# Node type
isdepot(n::Node) = typeof(n) == DepotNode
iscustomer(n::Node) = typeof(n) == CustomerNode

# Null route 
const NullRoute = Route(0, 0, 0, 0, 0, 0, 0, Inf)

# Empty (closed) route traversed by vehicle v from depot d                   
function Route(v::Vehicle, d::DepotNode)
    iʳ = length(v.R) + 1
    iᵛ = v.iᵛ
    iᵈ = d.iⁿ
    return Route(iʳ, iᵛ, iᵈ, iᵈ, iᵈ, 0, 0, 0)
end

# Hash solution
Base.hash(s::Solution) = hash(vectorize(s))