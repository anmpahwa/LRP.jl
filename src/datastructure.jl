"""
    Arc(iᵗ::Int, iʰ::Int, l::Float64)

An `Arc` is a connection between a tail node with index `iᵗ` and a head node with 
index `iʰ` spanning length `l`.
"""
struct Arc
    iᵗ::Int                                                                         # Tail node index
    iʰ::Int                                                                         # Head node index
    l::Float64                                                                      # Length
end



"""
    Route(iʳ::Int, iᵛ::Int, iᵈ::Int, x::Float64, y::Float64, iˢ::Int, iᵉ::Int, θⁱ::Float64, θˢ::Float64, θᵉ::Float64, tⁱ::Float64, tˢ::Float64, tᵉ::Float64, τ::Float64, n::Int, q::Float64, l::Float64, φ::Int64)

A `Route` is a connection between nodes, with index `iʳ`, vehicle index `iᵛ`, depot
node index `iᵈ`, centroid coordinates `(x, y)` , start node index `iˢ`, end node 
index `iᵉ`, vehicle tank status `θⁱ`, `θˢ`, and `θᵉ` at route initiaition time `tⁱ`, 
start time `tˢ`, and end time `tᵉ`, repsectively, slack time `τ`, customers served 
`n`, demand served `q`, length `l`, and initiaition status `φ` (usage in dynamic 
execution and simulation)
"""
mutable struct Route
    iʳ::Int                                                                         # Route index
    iᵛ::Int                                                                         # Vehicle index
    iᵈ::Int                                                                         # Depot node index
    x::Float64                                                                      # Centroid x-coordinate
    y::Float64                                                                      # Centroid y-coordinate
    iˢ::Int                                                                         # Start node index
    iᵉ::Int                                                                         # End node index
    θⁱ::Float64                                                                     # Vehicle tank status at the initiation time
    θˢ::Float64                                                                     # Vehicle tank status at the start time
    θᵉ::Float64                                                                     # Vehicle tank status at the end time
    tⁱ::Float64                                                                     # Initiation time
    tˢ::Float64                                                                     # Start time
    tᵉ::Float64                                                                     # End time
    τ::Float64                                                                      # Slack time
    n::Int                                                                          # Customers served
    q::Float64                                                                      # Demand served
    l::Float64                                                                      # Length
    φ::Int64                                                                        # Route initialization status
end



"""
    Vehicle(iᵛ::Int, jᵛ::Int, iᵈ::Int, qᵛ::Float64, lᵛ::Float64, sᵛ::Float64, τᶠ::Float64, τᵈ::Float64, τᶜ::Float64, τʷ::Float64, r̅::Int, πᵈ::Float64, πᵗ::Float64, πᶠ::Float64, tˢ::Float64, tᵉ::Float64, τ::Float64, n::Int, q::Float64, l::Float64, R::Vector{Route})

A `Vehicle` is a mode of delivery with index `iᵛ`, type index `jᵛ`, depot node index 
`iᵈ`, capacity `qᵛ`, range `lᵛ`, speed `sᵛ`, re-fueling time at the depot ndoe `τᶠ`, 
service time per package at the depot node `τᵈ`, parking time at a customer node 
`τᶜ`, driver working-hours `τʷ`, driver work-load (maximum vehicle-routes) `r̅`, 
operational cost `πᵈ` per unit distance and `πᵗ` per unit time, fixed cost `πᶠ`, 
initial departure time `tˢ`, final arrival time `tᵉ`, slack time `τ`, customers 
served `n`, demand served `q`, route length `l`, and set of routes `R`.
"""
mutable struct Vehicle
    iᵛ::Int                                                                         # Vehicle index
    jᵛ::Int                                                                         # Vehicle type index
    iᵈ::Int                                                                         # Depot node index
    qᵛ::Float64                                                                     # Capacity
    lᵛ::Float64                                                                     # Range
    sᵛ::Float64                                                                     # Speed
    τᶠ::Float64                                                                     # Re-fueling time at the depot node
    τᵈ::Float64                                                                     # service time per package at the depot node
    τᶜ::Float64                                                                     # Parking time at customer node
    τʷ::Float64                                                                     # Driver working-hours duration
    r̅::Int                                                                          # Driver work-load (maximum vehicle-routes)
    πᵈ::Float64                                                                     # Distance-based operational cost
    πᵗ::Float64                                                                     # Time-based operational cost
    πᶠ::Float64                                                                     # Fixed cost
    tˢ::Float64                                                                     # Start time (initial departure time from the depot node)
    tᵉ::Float64                                                                     # End time (final arrival time at the depot node)
    τ::Float64                                                                      # Slack time
    n::Int                                                                          # Customers served
    q::Float64                                                                      # Demand served
    l::Float64                                                                      # Route length
    R::Vector{Route}                                                                # Vector of vehicle routes
end



"""
    Node

A `Node` is a point on the graph.
"""
abstract type Node end
"""
    DepotNode(iⁿ::Int, x::Float64, y::Float64, qᵈ::Float64, tˢ::Float64, tᵉ::Float64, πᵒ::Float64, πᶠ::Float64, φ::Int, τ::Float64, n::Int, q::Float64, l::Float64, V::Vector{Vehicle})

A `DepotNode` is a source point on the graph at `(x,y)` with index `iⁿ`, capacity 
`qᵈ`, working-hours start time `tˢ` and end time `tᵉ`, operational cost `πᵒ` per 
package, fixed cost `πᶠ`, operations mandate `φ`, slack time `τ`, customers served 
`n`, demand served `q`, route length `l`, and fleet of vehicles `V`.
"""
mutable struct DepotNode <: Node
    iⁿ::Int                                                                         # Depot node index
    x::Float64                                                                      # Location on the x-axis
    y::Float64                                                                      # Location in the y-axis
    qᵈ::Float64                                                                     # Capacity
    tˢ::Float64                                                                     # Working-hours start time
    tᵉ::Float64                                                                     # Working-hours end time
    πᵒ::Float64                                                                     # Operational cost
    πᶠ::Float64                                                                     # Fixed cost
    φ::Int                                                                          # Operations mandate
    τ::Float64                                                                      # Depot slack time
    n::Int                                                                          # Customers served
    q::Float64                                                                      # Demand served
    l::Float64                                                                      # Route length
    V::Vector{Vehicle}                                                              # Vector of depot vehicles
end
"""
    CustomerNode(iⁿ::Int, x::Float64, y::Float64, qᶜ::Float64, τᶜ::Float64, tᵉ::Float64, tˡ::Float64, iᵗ::Int, iʰ::Int, tᵃ::Float64, tᵈ::Float64, r::Route)

A `CustomerNode` is a sink point on the graph at `(x,y)` with index `iⁿ`, demand 
`qᶜ`, customer service time (duration) `τᶜ`, earliest service time `tᵉ`, latest 
service time `tˡ`, tail node index `iᵗ`, head node index `iʰ`, arrival time `tᵃ`, 
and departure time `tᵈ`, serviced on route `r`.
"""
mutable struct CustomerNode <: Node
    iⁿ::Int                                                                         # Customer node index
    x::Float64                                                                      # Location on the x-axis
    y::Float64                                                                      # Location on the y-axis
    qᶜ::Float64                                                                     # Demand
    τᶜ::Float64                                                                     # Service time (duration)
    tᵉ::Float64                                                                     # Earliest service time
    tˡ::Float64                                                                     # Latest service time
    iᵗ::Int                                                                         # Tail (predecessor) node index
    iʰ::Int                                                                         # Head (successor) node index
    tᵃ::Float64                                                                     # Vehicle arrival time
    tᵈ::Float64                                                                     # Vehicle departure time
    r::Route                                                                        # Route visiting the customer
end



"""
    Solution(D::Vector{DepotNode}, C::OffsetVector{CustomerNode, Vector{CustomerNode}}, A::Dict{Tuple{Int,Int}, Arc}, πᶠ::Float64, πᵒ::Float64, πᵖ::Float64, φ::Bool)

A `Solution` is a graph with depot nodes `D`, customer nodes `C`, arcs `A`, fixed 
cost `πᶠ`, operational cost `πᵒ`, and penalty `πᵖ`. Note, `φ` is an internal binary 
trigger for en-route parameter evaluation.
"""
mutable struct Solution
    D::Vector{DepotNode}                                                            # Vector of depot nodes
    C::OffsetVector{CustomerNode, Vector{CustomerNode}}                             # Vector of customer nodes
    A::Dict{Tuple{Int,Int}, Arc}                                                    # Set of arcs
    πᶠ::Float64                                                                     # Fixed cost
    πᵒ::Float64                                                                     # Opertaional cost
    πᵖ::Float64                                                                     # Penalty
    φ::Bool                                                                         # En-route evaluation binary
end