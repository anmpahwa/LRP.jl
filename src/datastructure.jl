"""
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



"""
    Route(iʳ::Int64, iᵛ::Int64, iᵈ::Int64, x::Float64, y::Float64, iˢ::Int64, iᵉ::Int64, θⁱ::Float64, θˢ::Float64, θᵉ::Float64, tⁱ::Float64, tˢ::Float64, tᵉ::Float64, τ::Float64, n::Int64, q::Float64, l::Float64)

A `Route` is a connection between nodes, with index `iʳ`, vehicle index `iᵛ`, depot
node index `iᵈ`, centroid coordinate `x, y` , start node index `iˢ`, end node index 
`iᵉ`, vehicle tank status `θⁱ`, `θˢ`, and `θᵉ` at route initiaition `tⁱ`, start `tˢ`, 
and end time `tᵉ`, repsectively, slack time `τ`, customers served `n`, demand served 
`q`, and route length `l`.
"""
mutable struct Route
    iʳ::Int64                                                                       # Route index
    iᵛ::Int64                                                                       # Vehicle index
    iᵈ::Int64                                                                       # Depot node index
    x::Float64                                                                      # Centroid x-coordinate
    y::Float64                                                                      # Centroid y-coordinate
    iˢ::Int64                                                                       # Route start node index
    iᵉ::Int64                                                                       # Route end node index
    θⁱ::Float64                                                                     # Vehicle tank status at the initiation time
    θˢ::Float64                                                                     # Vehicle tank status at the start time
    θᵉ::Float64                                                                     # Vehicle tank status at the end time
    tⁱ::Float64                                                                     # Route initiation time
    tˢ::Float64                                                                     # Route start time
    tᵉ::Float64                                                                     # Route end time
    τ::Float64                                                                      # Route slack time
    n::Int64                                                                        # Route total customers served
    q::Float64                                                                      # Route total demand served
    l::Float64                                                                      # Route length
end



"""
    Vehicle(iᵛ::Int64, jᵛ::Int64, iᵈ::Int64, qᵛ::Float64, lᵛ::Float64, sᵛ::Float64, τᶠ::Float64, τᵈ::Float64, τᶜ::Float64, r̅::Int64, τʷ::Float64, tˢ::Float64, tᵉ::Float64, τ::Float64, n::Int64, q::Float64, l::Float64, πᵈ::Float64, πᵗ::Float64, πᶠ::Float64, R::Vector{Route})

A `Vehicle` is a mode of delivery with index `iᵛ`, vehicle type index `jᵛ`, depot 
node index `iᵈ`, capacity `qᵛ`, range `lᵛ`, speed `sᵛ`, refueling time `τᶠ`, 
service time `τᵈ` at depot node (per unit demand), parking time `τᶜ` at customer 
node, maximum number of vehicle routes permitted `r̅`, working-hours `τʷ`, initial 
departure time `tˢ`, final arrival time `tᵉ`, slack time `τ`, customers served `n`, 
demand served `q`, total route length `l`, operational cost `πᵈ` per unit distance 
and `πᵗ` per unit time, fixed cost `πᶠ`, and set of routes `R`.
"""
mutable struct Vehicle
    iᵛ::Int64                                                                       # Vehicle index
    jᵛ::Int64                                                                       # Vehicle type index
    iᵈ::Int64                                                                       # Depot node index
    qᵛ::Float64                                                                     # Vehicle capacity
    lᵛ::Float64                                                                     # Vehicle range
    sᵛ::Float64                                                                     # Vehicle speed
    τᶠ::Float64                                                                     # Re-fueling time
    τᵈ::Float64                                                                     # Depot node service time per unit demand
    τᶜ::Float64                                                                     # Parking time at customer stop
    τʷ::Float64                                                                     # Vehicle working-hours duration
    r̅::Int64                                                                        # Maximum number of vehicle routes permitted
    tˢ::Float64                                                                     # Vehicle start time (initial departure time from the depot node)
    tᵉ::Float64                                                                     # Vehicle end time (final arrival time at the depot node)
    τ::Float64                                                                      # Vehicle slack time
    n::Int64                                                                        # Vehicle total customers served
    q::Float64                                                                      # Vehicle total demand served
    l::Float64                                                                      # Vehicle total route length
    πᵈ::Float64                                                                     # Vehicle operational cost (distance based)
    πᵗ::Float64                                                                     # Vehicle operational cost (time based)
    πᶠ::Float64                                                                     # Vehicle fixed cost
    R::Vector{Route}                                                                # Vector of vehicle routes
end



"""
    Node

A `Node` is a point on the graph.
"""
abstract type Node end
"""
    DepotNode(iⁿ::Int64, jⁿ::Int64, x::Float64, y::Float64, q::Float64, pˡ::Float64, pᵘ::Float64, tˢ::Float64, tᵉ::Float64, τ::Float64, n::Int64, q::Float64, l::Float64, πᵒ::Float64, πᶠ::Float64, V::Vector{Vehicle})

A `DepotNode` is a source point on the graph at `(x,y)` with index `iⁿ` in echelon
`jⁿ`, capacity `q`, lower threshold `pˡ` and upper threshold `pᵘ` on share of 
customers handled, working-hours start time `tˢ` and end tme `tᵉ`, slack time `τ`, 
customers served `n`, demand served `q`, total route length `l`, operational cost 
`πᵒ` per package, fixed cost `πᶠ`, and fleet of vehicles `V`.
"""
mutable struct DepotNode <: Node
    iⁿ::Int64                                                                       # Depot node index
    jⁿ::Int64                                                                       # Depot echelon
    x::Float64                                                                      # Depot node location on the x-axis
    y::Float64                                                                      # Depot node location in the y-axis
    qᵈ::Float64                                                                     # Depot capacity
    pˡ::Float64                                                                     # Lower threshold on share of customers handled
    pᵘ::Float64                                                                     # Upper threshold on share of customers handled
    tˢ::Float64                                                                     # Depot working-hours start time
    tᵉ::Float64                                                                     # Depot working-hours end time
    τ::Float64                                                                      # Vehicle slack time
    n::Int64                                                                        # Vehicle total customers served
    q::Float64                                                                      # Vehicle total demand served
    l::Float64                                                                      # Vehicle total route length
    πᵒ::Float64                                                                     # Depot operational cost
    πᶠ::Float64                                                                     # Depot fixed cost
    V::Vector{Vehicle}                                                              # Vector of depot vehicles
end
"""
    CustomerNode(iⁿ::Int64, iʳ::Int64, iᵛ::Int64, iᵈ::Int64, x::Float64, y::Float64, q::Float64, τᶜ::Float64, tᵉ::Float64, tˡ::Float64, iᵗ::Int64, iʰ::Int64, tᵃ::Float64, tᵈ::Float64, r::Route)

A `CustomerNode` is a sink point on the graph at `(x,y)` with index `iⁿ`, demand `q`, 
customer service time `τᶜ`, earliest service time `tᵉ`, latest service time `tˡ`, 
tail node index `iᵗ`, head node index `iʰ`, arrival time `tᵃ`, departure time `tᵈ`, 
on route `r` with route index `iʳ`, vehicle index `iᵛ`, depot node index `iᵈ`.
"""
mutable struct CustomerNode <: Node
    iⁿ::Int64                                                                       # Customer node index
    iʳ::Int64                                                                       # Route index
    iᵛ::Int64                                                                       # Vehicle index
    iᵈ::Int64                                                                       # Depot node index
    x::Float64                                                                      # Customer node location on the x-axis
    y::Float64                                                                      # Customer node location in the y-axis
    q::Float64                                                                      # Customer demand
    τᶜ::Float64                                                                     # Customer service time
    tᵉ::Float64                                                                     # Customer node earliest service time
    tˡ::Float64                                                                     # Customer node latest service time
    iᵗ::Int64                                                                       # Tail (predecessor) node index
    iʰ::Int64                                                                       # Head (successor) node index
    tᵃ::Float64                                                                     # Customer node arrival time
    tᵈ::Float64                                                                     # Customer node departure time
    r::Route                                                                        # Route visiting the customer
end



"""
    Solution(D::Vector{DepotNode}, C::Vector{CustomerNode}, A::Dict{Tuple{Int64,Int64}, Arc}, V::Vector{Vehicle})

A Solution is a graph with depot nodes `D`, customer nodes `C`, arcs `A`, and vehicles `V`.
"""
struct Solution
    D::Vector{DepotNode}                                                            # Vector of depot nodes
    C::OffsetVector{CustomerNode, Vector{CustomerNode}}                             # Vector of customer nodes
    A::Dict{Tuple{Int64,Int64}, Arc}                                                # Set of arcs
end