@doc """
    Arc(iᵗ::Int64, iʰ::Int64, l::Float64)

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
    Route(iʳ::Int64, iᵛ::Int64, iᵈ::Int64, iˢ::Int64, iᵉ::Int64, θⁱ::Float64, θˢ::Float64, θᵉ::Float64, tⁱ::Float64, tˢ::Float64, tᵉ::Float64, τ::Float64, n::Int64, q::Int64, l::Float64, φ::Int64)

A `Route` is a connection between nodes, with index `iʳ`, vehicle index `iᵛ`, depot
node index `iᵈ`, start node index `iˢ`, end node index `iᵉ`, vehicle tank status 
`θⁱ`, `θˢ`, and `θᵉ` at route initiaition `tⁱ`, start `tˢ`, and end time `tᵉ`, 
repsectively, slack time `τ`, number of customers `n`, load `q`, length `l`, and
departure status `φ`.
"""
mutable struct Route
    iʳ::Int64                                                                       # Route index
    iᵛ::Int64                                                                       # Vehicle index
    iᵈ::Int64                                                                       # Depot node index
    iˢ::Int64                                                                       # Route start node index
    iᵉ::Int64                                                                       # Route end node index
    θⁱ::Float64                                                                     # Vehicle tank status at the initiation time
    θˢ::Float64                                                                     # Vehicle tank status at the start time
    θᵉ::Float64                                                                     # Vehicle tank status at the end time
    tⁱ::Float64                                                                     # Route initiation time
    tˢ::Float64                                                                     # Route start time
    tᵉ::Float64                                                                     # Route end time
    τ::Float64                                                                      # Route slack time
    n::Int64                                                                        # Route size (number of customers)
    q::Int64                                                                        # Route load
    l::Float64                                                                      # Route length
    φ::Int64                                                                        # Route departure status
end
    
@doc """
    Vehicle(iᵛ::Int64, jᵛ::Int64, iᵈ::Int64, q::Int64, l::Int64, s::Int64, τᶠ::Float64, τᵈ::Float64, τᶜ::Float64, r̅::Int64, τʷ::Int64, tˢ::Float64, tᵉ::Float64, πᵒ::Float64, πᶠ::Float64, R::Vector{Route})

A `Vehicle` is a mode of delivery with index `iᵛ`, vehicle type index `jᵛ`, depot 
node index `iᵈ`, capacity `q`, range `l`, speed `s`, refueling time `τᶠ`, service 
time `τᵈ` at depot node (per unit demand), service time `τᶜ` at customer node, 
maximum number of vehicle routes permitted `r̅`, working-hours `τʷ`, initial 
departure time `tˢ`, final arrival time `tᵉ`, operational cost `πₒ` per unit 
distance traveled, fixed cost `πᶠ`,  and set of routes `R`.
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
    τʷ::Int64                                                                       # Vehicle working-hours duration
    r̅::Int64                                                                        # Maximum number of vehicle routes permitted
    tˢ::Float64                                                                     # Vehicle start time (initial departure time from the depot node)
    tᵉ::Float64                                                                     # Vehicle end time (final arrival time at the depot node)
    πᵒ::Float64                                                                     # Vehicle operational cost
    πᶠ::Float64                                                                     # Vehicle fixed cost
    R::Vector{Route}                                                                # Vector of vehicle routes
end

@doc """
    Node

A `Node` is a point on the graph.
"""
abstract type Node end

@doc """
    DepotNode(iⁿ::Int64, jⁿ::Int64, x::Float64, y::Float64, q::Float64, pˡ::Float64, pᵘ::Float64, tˢ::Float64, tᵉ::Float64, πᵒ::Float64, πᶠ::Float64, V::Vector{Vehicle})

A `DepotNode` is a source point on the graph at `(x,y)` with index `iⁿ` in echelon
`jⁿ`, capacity `q`, lower threshold `pˡ` and upper threshold `pᵘ` on share of 
customers handled, working-hours start time `tˢ` and end tme  `tᵉ`,  operational 
cost  `πᵒ` per package, fixed cost `πᶠ`,  and fleet of vehicles `V`.
"""
struct DepotNode <: Node
    iⁿ::Int64                                                                       # Depot node index
    jⁿ::Int64                                                                       # Depot echelon
    x::Float64                                                                      # Depot node location on the x-axis
    y::Float64                                                                      # Depot node location in the y-axis
    q::Int64                                                                        # Depot capacity
    pˡ::Float64                                                                     # Lower threshold on share of customers handled
    pᵘ::Float64                                                                     # Upper threshold on share of customers handled
    tˢ::Float64                                                                     # Depot working-hours start time
    tᵉ::Float64                                                                     # Depot working-hours end time
    πᵒ::Float64                                                                     # Depot operational cost
    πᶠ::Float64                                                                     # Depot fixed cost
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
    Solution(D::Vector{DepotNode}, C::Vector{CustomerNode}, A::Dict{Tuple{Int64,Int64}, Arc}, V::Vector{Vehicle}, φᴱ::Int64, φᵀ::Int64)

A Solution is a graph with depot nodes `D`, customer nodes `C`, arcs `A`, and vehicles `V`.
"""
struct Solution
    D::Vector{DepotNode}                                                            # Vector of depot nodes
    C::OffsetVector{CustomerNode, Vector{CustomerNode}}                             # Vector of customer nodes
    A::Dict{Tuple{Int64,Int64}, Arc}                                                # Set of arcs
    φᴱ::Int64                                                                       # Binary (Internal use)
    φᵀ::Int64                                                                       # Binary (Internal use)
end

# is active
isactive(r::Route) = iszero(r.φ)                                                    # A route is said to be active if the vehicle hasn't departed yet

# is operational
isopt(r::Route) = (r.n ≥ 1)                                                         # A route is defined operational if it serves at least one customer
isopt(v::Vehicle) = any(isopt, v.R)                                                 # A vehicle is defined operational if any of its routes is operational
isopt(d::DepotNode) = any(isopt, d.V)                                               # A depot is defined operational if any of its vehicles is operational
isopen(c::CustomerNode) = isequal(c.r, NullRoute)                                   # A customer is defined open if it is not being served by any vehicle-route

# is close
isclose(d::DepotNode) = !isopt(d)                                                   # A depot node is defined closed if it is non-operational
isclose(c::CustomerNode) = !isopen(c)                                               # A customer node is defined closed it is being served by any vehicle-route

# is equal
Base.isequal(p::Route, q::Route) = isequal(p.iʳ, q.iʳ) && isequal(p.iᵛ, q.iᵛ) && isequal(p.iᵈ, q.iᵈ)
Base.isequal(p::Vehicle, q::Vehicle) = isequal(p.iᵛ, q.iᵛ) && isequal(p.iᵈ, q.iᵈ)
Base.isequal(p::Node, q::Node) = isequal(p.iⁿ, q.iⁿ)

# is identical
isidentical(v¹::Vehicle, v²::Vehicle) = isequal(v¹.jᵛ, v².jᵛ)

# Node type
isdepot(n::Node) = typeof(n) == DepotNode
iscustomer(n::Node) = typeof(n) == CustomerNode

# Null route
const NullRoute = Route(0, 0, 0, 0, 0, 0., 0., 0., Inf, Inf, Inf, 0., 0, 0, Inf, 0)

# Create a non-operational route traversed by vehicle v from depot d
function Route(v::Vehicle, d::DepotNode)
    iʳ = length(v.R) + 1
    iᵛ = v.iᵛ
    iᵈ = d.iⁿ
    iˢ = iᵈ
    iᵉ = iᵈ
    θⁱ = isone(iʳ) ? 1.0 : v.R[iʳ-1].θᵉ
    θˢ = θⁱ
    θᵉ = θˢ
    tⁱ = d.tˢ
    tˢ = d.tˢ
    tᵉ = d.tˢ
    τ  = Inf
    n  = 0 
    q  = 0
    l  = 0.
    φ  = 0
    r  = Route(iʳ, iᵛ, iᵈ, iˢ, iᵉ, θⁱ, θˢ, θᵉ, tⁱ, tˢ, tᵉ, τ, n, q, l, φ)
    return r
end            

# Create a non-operational vehicle cloning vehicle v at depot node d
function Vehicle(v::Vehicle, d::DepotNode)
    iᵛ = length(d.V) + 1
    jᵛ = v.jᵛ
    iᵈ = v.iᵈ
    q  = v.q
    l  = v.l
    s  = v.s
    τᶠ = v.τᶠ
    τᵈ = v.τᵈ
    τᶜ = v.τᶜ
    τʷ = v.τʷ
    r̅  = v.r̅
    tˢ = d.tˢ
    tᵉ = d.tˢ
    πᵒ = v.πᵒ
    πᶠ = v.πᶠ
    R  = Route[]
    v  = Vehicle(iᵛ, jᵛ, iᵈ, q, l, s, τᶠ, τᵈ, τᶜ, τʷ, r̅, tˢ, tᵉ, πᵒ, πᶠ, R)
    return v
end

# Hash solution
Base.hash(s::Solution) = hash(vectorize(s))