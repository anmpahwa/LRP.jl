@doc """
    Arc(iₜ::Int64, iₕ::Int64, l::Float64)

An `Arc` is a connection between tail   
node with index `iₜ` and head node with 
index `iₕ` with length `l`.
"""
struct Arc
    iₜ::Int64                                                                       # Tail node index
    iₕ::Int64                                                                       # Head node index
    l::Float64                                                                      # Arc length
end

@doc """
    Route(i::Int64, o::Int64, iₛ::Int64, iₑ::Int64, n::Int64, q::Int64, l::Float64)

A `Route` is a connection between nodes, with index `i`, origin vehicle index `o`, 
start node index `iₛ`, end node index `iₑ`, number of customers `n`, load `q`, and 
length `l`.
"""
mutable struct Route
    i::Int64                                                                        # Route index
    o::Int64                                                                        # Route (origin) vehicle index
    iₛ::Int64                                                                       # Route start node index
    iₑ::Int64                                                                       # Route end node index
    n::Int64                                                                        # Route size (number of customers)
    q::Int64                                                                        # Route load
    l::Float64                                                                      # Route length
end
    
@doc """
    Vehicle(i::Int64, o::Int64, q::Int64, πᵒ::Float64, πᶠ::Float64, R::Vector{Route})

A `Vehicle` is a mode of delivery with index `i`, origin depot node index `o`, 
capacity `q`, operational cost `πᵒ`, fixed cost `πᶠ`, and set of routes `R`.
"""
struct Vehicle
    i::Int64                                                                        # Vehicle index
    o::Int64                                                                        # Vehicle origin (depot node) index
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
    DepotNode(i::Int64, x::Float64, y::Float64, q::Float64, πᵒ::Float64, πᶠ::Float64, V::Vector{Vehicle})

A `DepotNode` is a source point on the graph at `(x,y)` with index `i`, capacity `q`, operational cost 
`πᵒ`, fixed cost `πᶠ`, and fleet of vehicles `V`.
"""
struct DepotNode <: Node
    i::Int64                                                                        # Depot node index
    x::Float64                                                                      # Depot node location on the x-axis
    y::Float64                                                                      # Depot node location in the y-axis
    q::Int64                                                                        # Depot capacity
    πᵒ::Float64                                                                     # Operational cost
    πᶠ::Float64                                                                     # Fixed cost
    V::Vector{Vehicle}                                                              # Vector of depot vehicles
end

@doc """
    CustomerNode(i::Int64, x::Float64, y::Float64, q::Float64, iₜ::Int64, iₕ::Int64, r::Route)

A `CustomerNode` is a sink point on the graph at `(x,y)` with index `i`, demand `q`, tail node 
index `iₜ`, head node index `iₕ`, on route `r`.
"""
mutable struct CustomerNode <: Node
    i::Int64                                                                        # Customer node index
    x::Float64                                                                      # Customer node location on the x-axis
    y::Float64                                                                      # Customer node location in the y-axis
    q::Int64                                                                        # Customer demand
    iₜ::Int64                                                                       # Tail (predecessor) node index
    iₕ::Int64                                                                       # Head (successor) node index
    r::Route                                                                        # Route visiting the customer
end

@doc """
    Solution(D::Vector{DepotNode}, C::Vector{CustomerNode}, A::Dict{Tuple{Int64,Int64}, Arc}, V::Vector{Vehicle})

A Solution is a graph with depot nodes `D`, customer nodes `C`, arcs `A`, and vehicles `V`.
"""
struct Solution
    D::Vector{DepotNode}                                                            # Vector of depot nodes
    C::OffsetVector{CustomerNode, Vector{CustomerNode}}                             # Vector of customer nodes
    A::Dict{Tuple{Int64,Int64}, Arc}                                                # Set of arcs
    V::Vector{Vehicle}                                                              # Set of vehicles
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
Base.isequal(p::Route, q::Route) = isequal(p.i, q.i)
Base.isequal(p::Vehicle, q::Vehicle) = isequal(p.i, q.i)
Base.isequal(p::Node, q::Node) = isequal(p.i, q.i)

# Node type
isdepot(n::Node) = typeof(n) == DepotNode
iscustomer(n::Node) = typeof(n) == CustomerNode

# Null route 
const NullRoute = Route(0, 0, 0, 0, 0, 0, Inf)

# Empty (closed) route traversed by vehicle v from depot d                   
Route(i, v::Vehicle, d::DepotNode) = Route(i, v.i, d.i, d.i, 0, 0, 0)

# Hash solution
Base.hash(s::Solution) = hash(vectorize(s))