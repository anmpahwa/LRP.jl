[![Build Status](https://github.com/anmol1104/LRP.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/anmol1104/LRP.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/anmol1104/LRP.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/anmol1104/LRP.jl)

# Location Routing Problem (LRP)

Given, a graph `G = (D, C, A, V)` with 
set of depots `D` with capacity `d.q`, fleet `d.V`, operational cost `d.πᵒ`, and fixed cost `d.πᶠ` for every depot `d ∈ D`; 
set of customer nodes `C` with demand `c.q` for every customer `c ∈ C`;
set of arcs `A = {(i,j); i,j ∈ N={D∪C}}` with length `l` for every arc `(i,j) ∈ A`; and 
set of vehicles `V` with capacity `v.q`, operational cost `v.πᵒ`, and fixed cost `v.πᶠ` for every vehicle `v ∈ V`, 
the objective is to develop least cost routes from select depot nodes using select vehicles such that every customer node is visited exactly once while also accounting for depot and vehicle capacities.  

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
- Clarke and Wright Savings Algorithm   : `:cw`
- Nearest Neighborhood Algorithm        : `:nn`
- Random Initialization                 : `:random`
- Regret N Insertion                    : `:regret₂init`, `:regret₃init`

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
- Best Insertion    : `best!`
- Greedy Insertion  : `greedy!`
- Regret Insertion  : `regret₂insert!`, `regret₃insert!`

In every few iterations, the ALNS metaheuristic performs local search with,
- Move      : `:move!`
- Inter-Opt : `:interopt!`
- Intra-Opt : `:intraopt!`
- Split     : `:split!`
- Swap      : `:swap!`

See example.jl for usage