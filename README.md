# Location Routing Problem (LRP)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://github.com/anmol1104/LRP.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/anmol1104/LRP.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/anmol1104/LRP.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/anmol1104/LRP.jl)

capacitated location routing problem with time-windows with heterogeneous fleet of multi-route delivery vehicles

Given, a graph `G = (D, C, A)` with set of depots `D` with capacity `d.q`, lower threshold `d.pˡ` and upper threshold `d.pᵘ` on share of customers handled, working-hours start time `d.tˢ` and end tme  `d.tᵉ`,  operational cost  `d.πᵒ` per package, fixed cost `d.πᶠ`, mandated depot use `d.φ`, and fleet of vehicles `d.V` with capacity `v.q`, range `v.l`, speed `v.s`, refueling time `v.τᶠ`, depot node service time `v.τᵈ` (per unit demand), customer node parking time `v.τᶜ`, driver working hours `v.τʷ`, maximum number of vehicle routes permitted `v.r̅`, operational cost `v.πᵈ` per unit distance and `v.πᵗ` per unit time, fixed cost `v.πᶠ`, and  for every vehicle `v ∈ d.V`, for every depot `d ∈ D`; set of customer nodes `C` with demand `c.q`, service time `c.τᶜ`, delivery time-window `[c.tᵉ,c.tˡ]` for every customer `c ∈ C`; set of arcs `A` with length `l` for every arc `(i,j) ∈ A`; the objective is to develop least cost routes from select depot nodes using select vehicles such that every customer node is visited exactly once while also accounting for depot capacity, vehicle capacity, vehicle range, driver working-hours, and customers' time-windows.

This package uses Adaptive Large Neighborhood Search (ALNS) algorithm to find an optimal solution for the Location Routing Problem given ALNS optimization parameters,
- j     :   Number of segments in the ALNS
- k     :   Number of segments to reset ALNS
- n     :   Number of iterations in an ALNS segment
- m     :   Number of iterations in a local search
- Ψᵣ    :   Vector of removal operators
- Ψᵢ    :   Vector of insertion operators
- Ψₗ    :   Vector of local search operators
- σ₁    :   Score for a new best solution
- σ₂    :   Score for a new better solution
- σ₃    :   Score for a new worse but accepted solution
- μ̲     :   Minimum removal fraction
- C̲     :   Minimum customer nodes removed
- μ̅     :   Maximum removal fraction
- C̅     :   Maximum customer nodes removed
- ω̅     :   Initial temperature deviation parameter
- τ̅     :   Initial temperatureprobability parameter
- ω̲     :   Final temperature deviation parameter
- τ̲     :   Final temperature probability parameter
- φ     :   Local search trigger
- θ     :   Cooling rate
- ρ     :   Reaction factor

and an initial solution developed using iterated clustering method,

The ALNS metaheuristic iteratively removes a set of nodes using,
- Random Customer Node Removal  : `:randomcustomer!`
- Random Route Removal          : `:randomroute!`
- Random Vehicle Removal        : `:randomvehicle!`
- Random Depot Removal          : `:randomdepot!` 
- Related Customer Node Removal : `:relatedcustomer!`
- Related Route removal         : `:relatedroute!`
- Related Vehicle Removal       : `:relatedvehicle!`
- Related Depot Removal         : `:relateddepot!`
- Worst Customer Node Removal   : `:worstcustomer!`
- Worst Route Removal           : `:worstroute!`
- Worst Vehicle Removal         : `:worstvehicle!`
- Worst Depot Removal           : `:worstdepot!`

and consequently inserts removed nodes using,
- Best Insertion           : `:best!`
- Precise Greedy Insertion : `:precise!`
- Perturb Greedy insertion : `:perturb!`
- Regret-two Insertion     : `:regret2!`
- Regret-three Insertion   : `:regret3!`

In every few iterations, the ALNS metaheuristic performs local search with,
- intra-move    : `:intramove!`
- inter-move    : `:intermove!`
- intra-swap    : `:intraswap!`
- inter-swap    : `:interswap!`
- intra-opt     : `:intraopt!`
- inter-opt     : `:interopt!`
- swapdepot     : `:swapdepot!`

See benchmark.jl for usage

Additional initialization, removal, insertion, and local search methods can be defined.