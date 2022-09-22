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
    Route(iʳ::Int64, iᵛ::Int64, iᵈ::Int64, iˢ::Int64, iᵉ::Int64, tˢ::Float64, tᵉ::Float64, n::Int64, q::Int64, l::Float64)

A `Route` is a connection between nodes, with index `iʳ`, vehicle index `iᵛ`, depot
node index `iᵈ`, start node index `iˢ`, end node index `iᵉ`, start time `tₛ`, end 
time `tₑ`, number of customers `n`, load `q`, and length `l`.
"""
mutable struct Route
    iʳ::Int64                                                                       # Route index
    iᵛ::Int64                                                                       # Vehicle index
    iᵈ::Int64                                                                       # Depot node index
    iˢ::Int64                                                                       # Route start node index
    iᵉ::Int64                                                                       # Route end node index
    tˢ::Float64                                                                     # Route start time (departure time from the depot node)
    tᵉ::Float64                                                                     # Route end time (subsequent arrival time at the depot node)
    n::Int64                                                                        # Route size (number of customers)
    q::Int64                                                                        # Route load
    l::Float64                                                                      # Route length
end
    
@doc """
    Vehicle(iᵛ::Int64, jᵛ::Int64, iᵈ::Int64, q::Int64, l::Int64, s::Int64, τᶠ::Float64, τᵈ::Float64, τᶜ::Float64, πᵒ::Float64, πᶠ::Float64, r̅::Int64, w::Int64, tˢ::Float64, tᵉ::Float64, R::Vector{Route})

A `Vehicle` is a mode of delivery with index `iᵛ`, vehicle type index `jᵛ`, depot 
node index `iᵈ`, capacity `q`, range `l`, speed `s`, refueling time `τᶠ`, service 
time `τᵈ` at depot node (per unit demand), service time `τᶜ` at customer node, 
operational cost `πₒ` per unit distance traveled, fixed cost `πᶠ`, maximum number
of vehicle routes permitted `r̅`, working hours `w`, start time `tˢ`, end time `tᵉ`, 
and set of routes `R`.
"""
mutable struct Vehicle
    iᵛ::Int64                                                                       # Vehicle index
    jᵛ::Int64                                                                       # Vehicle type index
    iᵈ::Int64                                                                       # Depot node index
    q::Int64                                                                        # Vehicle capacity
    l::Int64                                                                        # Vehicle range
    s::Int64                                                                        # Vehicle speed
    τᶠ::Float64                                                                     # Re-fueling time
    τᵈ::Float64                                                                     # Depot node service time per unit demand
    τᶜ::Float64                                                                     # Customer node service time
    πᵒ::Float64                                                                     # Operational cost
    πᶠ::Float64                                                                     # Fixed cost
    r̅::Int64                                                                        # Maximum number of vehicle routes permitted
    w::Int64                                                                        # Working hours
    tˢ::Float64                                                                     # Vehicle start time (initial departure time from the depot node)
    tᵉ::Float64                                                                     # Vehicle end time (final arrival time at the depot node)
    R::Vector{Route}                                                                # Vector of vehicle routes
end

@doc """
    Node

A `Node` is a point on the graph.
"""
abstract type Node end

@doc """
    DepotNode(i::Int64, x::Float64, y::Float64, q::Float64, V::Vector{Vehicle}, ϕ::Int64)

A `DepotNode` is a source point on the graph at `(x,y)` with index `iⁿ`, capacity 
`q`, operational cost `πᵒ` per package, fixed cost `πᶠ`, and fleet of vehicles `V`.
"""
struct DepotNode <: Node
    iⁿ::Int64                                                                       # Depot node index
    x::Float64                                                                      # Depot node location on the x-axis
    y::Float64                                                                      # Depot node location in the y-axis
    q::Int64                                                                        # Depot capacity
    πᵒ::Float64                                                                     # Operational cost
    πᶠ::Float64                                                                     # Fixed cost
    V::Vector{Vehicle}                                                              # Vector of depot vehicles
end

@doc """
    CustomerNode(iⁿ::Int64, x::Float64, y::Float64, q::Float64, tᵉ::Float64, tˡ::Float64, iᵗ::Int64, iʰ::Int64, tᵃ::Float64, tᵈ::Float64, r::Route)

A `CustomerNode` is a sink point on the graph at `(x,y)` with index `iⁿ`, demand `q`, 
earliest service time `tᵉ`, latest service time `tˡ`, tail node index `iᵗ`, head node 
index `iʰ`, arrival time `tᵃ`, departure time `tᵈ`, on route `r`.
"""
mutable struct CustomerNode <: Node
    iⁿ::Int64                                                                       # Customer node index
    x::Float64                                                                      # Customer node location on the x-axis
    y::Float64                                                                      # Customer node location in the y-axis
    q::Int64                                                                        # Customer demand
    tᵉ::Float64                                                                     # Customer node earliest service time
    tˡ::Float64                                                                     # Customer node latest service time
    iᵗ::Int64                                                                       # Tail (predecessor) node index
    iʰ::Int64                                                                       # Head (successor) node index
    tᵃ::Float64                                                                     # Customer node arrival time
    tᵈ::Float64                                                                     # Customer node departure time
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

# is identical
isidentical(v¹::Vehicle, v²::Vehicle) = isequal(v¹.jᵛ, v².jᵛ)

# Node type
isdepot(n::Node) = typeof(n) == DepotNode
iscustomer(n::Node) = typeof(n) == CustomerNode

# Null route
const NullRoute = Route(0, 0, 0, 0, 0, Inf, Inf, 0, 0, Inf)

# Create a non-operational route traversed by vehicle v from depot d
function Route(v::Vehicle, d::DepotNode)
    iʳ = length(v.R) + 1
    iᵛ = v.iᵛ
    iᵈ = d.iⁿ
    r  = Route(iʳ, iᵛ, iᵈ, iᵈ, iᵈ, 0., 0., 0, 0, 0)
    return r
end            

# Create a non-operational vehicle cloning vehicle v at depot node d
function Vehicle(v::Vehicle, d::DepotNode)
    iᵛ = length(d.V) + 1
    v  = Vehicle(iᵛ, v.jᵛ, v.iᵈ, v.q, v.l, v.s, v.τᶠ, v.τᵈ, v.τᶜ, v.πᵒ, v.πᶠ, v.r̅, v.w, 0., 0., Route[])
    return v
end

# Hash solution
Base.hash(s::Solution) = hash(vectorize(s))