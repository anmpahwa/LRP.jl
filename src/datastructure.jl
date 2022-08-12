@doc """
    Arc(i::Int64, j::Int64, l::Float64, t::Float64, f::Float64)

An `Arc` is a connection between node `i` and `j` with  arc length `l`, travel 
time `t`, and fuel use `f`.
"""
struct Arc
    i::Int64                                                                        # Tail node index
    j::Int64                                                                        # Head node index
    l::Float64                                                                      # Arc length
    t::Float64                                                                      # Arc travel time
    f::Float64                                                                      # Arc fuel use
end

@doc """
    Route(i::Int64, o::Int64, s::Int64, e::Int64, n::Int64, q::Int64, l::Float64, t::Float64, f::Float64, c::Float64)

A `Route` is a ... with index value `i`, origin vehicle index `o`, start node 
index `s`, end node index `e`, number of customers `n`, load `q`, length `l`, 
travel time `t`, fuel use `f`, and cost `c`.
"""
mutable struct Route
    i::Int64                                                                        # Route index
    o::Int64                                                                        # Route origin (vehicle) index
    s::Int64                                                                        # Route start node index
    e::Int64                                                                        # Route end node index
    n::Int64                                                                        # Route size (number of customers)
    q::Int64                                                                        # Route load
    l::Float64                                                                      # Route length
    t::Float64                                                                      # Route travel time
    f::Float64                                                                      # Route fuel use
    c::Float64                                                                      # Route cost
end
    
# TODO: Make vehicles completely heterogenous (add vehicle efficiency and speed)

@doc """
    Vehicle(i::Int64, o::Int64, q::Int64, πᵐ::Float64, πʷ::Float64, πᶠ::Float64, πᵛ::Float64, R::Vector{Route})

A `Vehicle` is a delivery vehicle with index value `i`, origin depot node index 
`o`, vehicle capacity `q`, maintenance cost `πᵐ`, driver's wage `πʷ`, fuel cost 
`πᶠ`, rental cost `πᵛ`, and set of routes `R`.
"""
struct Vehicle
    i::Int64                                                                        # Vehicle index
    o::Int64                                                                        # Vehicle origin (depot node) index
    q::Int64                                                                        # Vehicle capacity
    πᵐ::Float64                                                                     # Maintenance cost ($ per unit distance)
    πʷ::Float64                                                                     # Driver wage cost ($ per unit time)
    πᶠ::Float64                                                                     # Fuel/Energy cost ($ per unit energy)
    πᵛ::Float64                                                                     # Vehicle rental cost ($)
    R::Vector{Route}                                                                # Vector of vehicle routes
end

@doc """
    Node

A `Node` is a point on the graph.
"""
abstract type Node end

@doc """
    DepotNode(i::Int64, x::Float64, y::Float64, q::Float64, πᵈ::Float64, V::Vector{Vehicle})

A `DepotNode` is a source point on the graph at `(x,y)` with index value `i`, 
capacity `q`, rental cost `πᵈ`, and fleet of vehicles `V`.
"""
struct DepotNode <: Node
    i::Int64                                                                        # Depot node index
    x::Float64                                                                      # Depot node location on the x-axis
    y::Float64                                                                      # Depot node location in the y-axis
    q::Int64                                                                        # Depot capacity
    πᵈ::Float64                                                                     # Depot rental cost ($)
    V::Vector{Vehicle}                                                              # Vector of depot vehicles
end

@doc """
    CustomerNode(i::Int64, x::Float64, y::Float64, q::Float64, t::Int64, h::Int64, r::Route)

A `CustomerNode` is a sink point on the graph at `(x,y)` with index value `i`, 
demand `q`, tail node index `t`, head node index `h`, and on route `r`.
"""
mutable struct CustomerNode <: Node
    i::Int64                                                                        # Customer node index
    x::Float64                                                                      # Customer node location on the x-axis
    y::Float64                                                                      # Customer node location in the y-axis
    q::Int64                                                                        # Customer demand
    t::Int64                                                                        # Tail (predecessor) node index
    h::Int64                                                                        # Head (successor) node index
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
    V::Vector{Vehicle}                                                              # Set of vehicles
end

# Route definitions
Route() = Route(0, 0, 0, 0, 0, 0, Inf, Inf, Inf, Inf)                               # Null route      
Route(i, v::Vehicle, d::DepotNode) = Route(i, v.i, d.i, d.i, 0, 0, 0, 0, 0, 0)      # Empty (closed) route traversed by vehicle v from depot d

# isequal
Base.isequal(p::Route, q::Route) = isequal(p.i, q.i)
Base.isequal(p::Vehicle, q::Vehicle) = isequal(p.i,q.i)
Base.isequal(p::Node, q::Node) = isequal(p.i,q.i)

# isclose
isclose(r::Route) = iszero(r.n)                                                     # A route is defined closed if it serves no customer
isclose(v::Vehicle) = all(isclose, v.R)                                             # A vehicle is defined closed if all its routes are closed
isclose(d::DepotNode) = all(isclose, d.V)                                           # A depot is defined closed if all its vehicles are closed
isclose(c::CustomerNode) = !(iszero(c.t) & iszero(c.h))                             # A customer is defined closed if its tail node and head node index is non-zero
isopen(x) = !isclose(x)

# node type
isdepot(n::Node) = typeof(n) == DepotNode
iscustomer(n::Node) = typeof(n) == CustomerNode

# Objective function evaluation
"""
    f(r::Route, χₒ::ObjectiveFunctionParameters)

Objective function evaluation for route `r` using objective function 
parameters `χₒ`.
"""
f(r::Route, χ::ObjectiveFunctionParameters) = r.c
"""
    f(v::Vehicle, χₒ::ObjectiveFunctionParameters)

Objective function evaluation for vehicle `v` using objective function 
parameters `χₒ`.
"""
function f(v::Vehicle, χₒ::ObjectiveFunctionParameters)
    if isclose(v) return 0. end
    α = χₒ.v
    z = v.πᵛ
    q = 0
    for r ∈ v.R 
        z += f(r, χₒ)
        q += r.q
    end
    if q > v.q z += α * (q - v.q) end
    return z
end
"""
    f(d::DepotNode, χₒ::ObjectiveFunctionParameters)

Objective function evaluation for depot node `d` using objective function 
parameters `χₒ`.
"""
function f(d::DepotNode, χₒ::ObjectiveFunctionParameters)
    if isclose(d) return 0. end
    α = χₒ.d
    z = d.πᵈ
    q = 0
    for v ∈ d.V
        z += f(v, χₒ)
        for r ∈ v.R q += r.q end 
    end
    if q > d.q z += α * (q - d.q) end
    return z
end
"""
    f(s::Solution, χₒ::ObjectiveFunctionParameters)

Objective function evaluation for solution `s` using objective function 
parameters `χₒ`.
"""
function f(s::Solution, χₒ::ObjectiveFunctionParameters)
    z = 0.
    for d ∈ s.D return z += f(d, χₒ) end
    return z
end
"""
    infeasible(s::Solution)

Returns true if node service constraint, node flow constraint, and
sub-tour elimination constraint are not violated.
"""
function isfeasible(s::Solution)
    D = s.D
    C = s.C
    V = s.V
    # Customer node service and flow constraints
    x = zeros(Int64, eachindex(C))
    for d ∈ D
        for v ∈ V
            R = v.R
            for r ∈ R
                if isclose(r) continue end
                nₛ = C[r.s]
                nₑ = C[r.e]
                nₒ = nₛ
                while true
                    k = nₒ.i
                    x[k] += 1
                    if isequal(nₒ, nₑ) break end
                    nₒ = C[nₒ.h]
                end
            end
        end
    end
    if any(!isone, x) return false end
    # Vehicle capacity constraint
    for v ∈ V
        q = 0
        R = v.R
        for r ∈ R q += r.q end
        if q > v.q return false end
    end
    return true
end