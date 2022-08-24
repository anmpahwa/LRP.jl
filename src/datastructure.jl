@doc """
    Arc(i::Int64, j::Int64, l::Float64, t::Float64, f::Float64)

An `Arc` is a connection between node `i` and `j` with length `l`.
"""
struct Arc
    i::Int64                                                                        # Tail node index
    j::Int64                                                                        # Head node index
    l::Float64                                                                      # Arc length
end

@doc """
    Route(i::Int64, o::Int64, s::Int64, e::Int64, n::Int64, q::Int64, l::Float64, t::Float64, f::Float64, c::Float64)

A `Route` is a ... with index value `i`, origin vehicle index `o`, start node 
index `s`, end node index `e`, number of customers `n`, load `q`, and length `l`.
"""
mutable struct Route
    i::Int64                                                                        # Route index
    o::Int64                                                                        # Route origin (vehicle) index
    s::Int64                                                                        # Route start node index
    e::Int64                                                                        # Route end node index
    n::Int64                                                                        # Route size (number of customers)
    q::Int64                                                                        # Route load
    l::Float64                                                                      # Route length
end
    
@doc """
    Vehicle(i::Int64, o::Int64, q::Int64, πᵐ::Float64, πʷ::Float64, πᶠ::Float64, πᵛ::Float64, R::Vector{Route})

A `Vehicle` is a mode of delivery with index value `i`, origin depot node index 
`o`, capacity `q`, operational cost `πᵒ`, fixed cost `πᶠ`, and set of routes `R`.
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
    DepotNode(i::Int64, x::Float64, y::Float64, q::Float64, πᵈ::Float64, V::Vector{Vehicle})

A `DepotNode` is a source point on the graph at `(x,y)` with index value `i`, 
capacity `q`, operational cost `πᵒ`, fixed cost `πᶠ`, and fleet of vehicles `V`.
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

A Solution is a graph with depot nodes `D`, customer nodes `C`, arcs `A`, and vehicles `V`.
"""
struct Solution
    D::Vector{DepotNode}                                                            # Vector of depot nodes
    C::OffsetVector{CustomerNode, Vector{CustomerNode}}                             # Vector of customer nodes
    A::Dict{Tuple{Int64,Int64}, Arc}                                                # Set of arcs
    V::Vector{Vehicle}                                                              # Set of vehicles
    R::Vector{Route}                                                                # Set of routes
end

# Route definitions
const NullRoute = Route(0, 0, 0, 0, 0, 0, Inf)                                      # Null route      
Route(i, v::Vehicle, d::DepotNode) = Route(i, v.i, d.i, d.i, 0, 0, 0)               # Empty (closed) route traversed by vehicle v from depot d

# is equal
Base.isequal(p::Route, q::Route) = isequal(p.i, q.i)
Base.isequal(p::Vehicle, q::Vehicle) = isequal(p.i,q.i)
Base.isequal(p::Node, q::Node) = isequal(p.i,q.i)

# is operational
isopt(r::Route) = (r.n ≥ 1)                                                         # A route is defined operational if it serves at least one customer
isopt(v::Vehicle) = any(isopt, v.R)                                                 # A vehicle is defined operational if any of its routes is operational
isopt(d::DepotNode) = any(isopt, d.V)                                               # A depot is defined operational if any of its vehicles is operational
isopen(c::CustomerNode) = isequal(c.r, NullRoute)                                   # A customer is defined open if it is not being served by any vehicle-route

# is close
isclose(d::DepotNode) = !isopt(d)
isclose(c::CustomerNode) = !isopen(c)

# Node type
isdepot(n::Node) = typeof(n) == DepotNode
iscustomer(n::Node) = typeof(n) == CustomerNode

# Objective function evaluation
"""
    f(s::Solution)

Objective function evaluation for solution `s`.
"""
function f(s::Solution)
    z = 0.
    for d ∈ s.D
        if !isopt(d) continue end 
        qᵈ = 0
        z += d.πᶠ
        for v ∈ d.V 
            if !isopt(v) continue end
            z += v.πᶠ
            for r ∈ v.R 
                if !isopt(r) continue end
                qʳ  = r.q
                qᵈ += qʳ
                z  += r.l * v.πᵒ
                z  += z * (qʳ > v.q) * (qʳ - v.q)
            end
        end
        z += qᵈ * d.πᵒ
        z += z * (qᵈ > d.q) * (qᵈ - d.q)
    end
    return z 
end

# Solution feasibility
"""
    infeasible(s::Solution)

Returns true if node service constraint, node flow constraint, and
sub-tour elimination constraint are not violated.
"""
function isfeasible(s::Solution)
    D = s.D
    C = s.C
    # Customer node service and flow constraints
    x = zeros(Int64, eachindex(C))
    for d ∈ D
        if !isopt(d) continue end 
        V = d.V
        for v ∈ V
            if !isopt(v) continue end 
            R = v.R
            for r ∈ R
                if !isopt(r) continue end
                cₛ = C[r.s]
                cₑ = C[r.e]
                c  = cₛ
                while true
                    k = c.i
                    x[k] += 1
                    if isequal(c, cₑ) break end
                    c = C[c.h]
                end
            end
        end
    end
    if any(!isone, x) return false end
    # Capacity constraints
    for d ∈ D
        if !isopt(d) continue end 
        qᵈ = 0
        for v ∈ d.V
            if !isopt(v) continue end 
            for r ∈ v.R 
                if !isopt(r) continue end
                qʳ  = r.q
                qᵈ += qʳ
                if qʳ > v.q return false end
            end
        end
        if qᵈ > d.q return false end
    end
    return true
end

# Hash solution
Base.hash(s::Solution) = hash(vectorize(s))