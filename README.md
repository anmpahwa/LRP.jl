[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Build Status](https://github.com/anmol1104/LRP.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/anmol1104/LRP.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/anmol1104/LRP.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/anmol1104/LRP.jl)

# Location Routing Problem (LRP)

capacitated location routing problem with time-windows with heterogeneous fleet of multi-route delivery vehicles

Given, a graph `G = (D, C, A)` with 
set of depots `D` with capacity `d.q`, operational cost `d.πᵒ`, fixed cost `d.πᶠ`, and
fleet of vehicles `d.V` with capacity `v.q`, range `v.l`, speed `v.s`, refueling time `v.τᶠ`, depot node service time `v.τᵈ` (per unit demand), customer node service time `v.τᶜ`, operational cost `v.πᵒ`, fixed cost `v.πᶠ`, and driver working hours `v.w` for every vehicle `v ∈ d.V`, 
for every depot `d ∈ D`;
set of customer nodes `C` with demand `c.q`, delivery time-window `[c.tᵉ,c.tˡ]` for every customer `c ∈ C`;
set of arcs `A` with length `l` for every arc `(i,j) ∈ A`; and 
the objective is to develop least cost routes from select depot nodes using select vehicles such that every customer node is visited exactly once while also accounting for depot capacity, vehicle capacity, vehicle range, driver working-hours, and customers' time-windows.

This package uses Adaptive Large Neighborhood Search (ALNS) algorithm to find an optimal solution for the Locatio Routing Problem given ALNS optimization 
parameters,
- `k̲`     :   Number of ALNS iterations triggering operator probability update (segment size)
- `l̲`     :   Number of ALNS iterations triggering local search
- `l̅`     :   Number of local search iterations
- `k̅`     :   Number of ALNS iterations
- `Ψᵣ`    :   Vector of removal operators
- `Ψᵢ`    :   Vector of insertion operators
- `Ψₗ`    :   Vector of local search operators
- `σ₁`    :   Score for a new best solution
- `σ₂`    :   Score for a new better solution
- `σ₃`    :   Score for a new worse but accepted solution
- `ω`     :   Start tempertature control threshold 
- `τ`     :   Start tempertature control probability
- `𝜃`     :   Cooling rate
- `C̲`     :   Minimum customer nodes removal
- `C̅`     :   Maximum customer nodes removal
- `μ̲`     :   Minimum removal fraction
- `μ̅`     :   Maximum removal fraction
- `ρ`     :   Reaction factor

and an initial solution developed using one of the following methods,
- K-means Clustering Intialization  : `:cluster`
- Random Initialization             : `:random`

The ALNS metaheuristic iteratively removes a set of nodes using,
- Random Node Removal       : `:randomnode!`
- Random Route Removal      : `:randomroute!`
- Random Vehicle Removal    : `:randomvehicle!`
- Random Depot Removal      : `:randomdepot!` 
- Related Node Removal      : `:relatednode!`
- Related Route removal     : `:relatedroute!`
- Related Vehicle Removal   : `:relatedvehicle!`
- Related Depot Removal     : `:relateddepot!`
- Worst Node Removal        : `:worstnode!`
- Worst Route Removal       : `:worstroute!`
- Worst Vehicle Removal     : `:worstvehicle!`
- Worst Depot Removal       : `:worstdepot!`

and consequently inserts removed nodes using,
- Precise Best Insertion    : `:bestprecise!`
- Perturb Best Insertion    : `:bestperturb!`
- Precise Greedy Insertion  : `:greedyprecise!`
- Perturb Greedy Insertion  : `:greedyperturb!`
- Regret-two Insertion      : `:regret2!`
- Regret-three Insertion    : `:regret3!`

In every few iterations, the ALNS metaheuristic performs local search with,
- Move          : `:move!`
- Opt           : `:opt!`
- Split         : `:split!`
- Swap-Customer : `:swapcustomers!`
- Swap-Depot    : `:swapdepots!`

See example.jl for usage

Additional initialization, removal, insertion, and local search methods can be defined.