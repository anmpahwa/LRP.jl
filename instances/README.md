# Instance files

Each instance folder contains the following files

- arcs.csv: matrix of arc lengths

- customer_nodes.csv: list of customer nodes in the network
    - in: customer node index
    - x: location on the x-axis
    - y: location on the y-axis
    - qc: demand
    - tc: service time (duration)
    - te: earliest service time
    - tl: latest service time

- depot_nodes.csv: list of depot nodes in the network
    - in: depot node index
    - x: location on the x-axis
    - y: location on the y-axis
    - qd: capacity
    - ts: working-hours start time
    - te: working-hours end time
    - co: operational cost per package
    - cf: fixed cost
    - phi: operations mandate

- vehicles.csv: list of vehicles at the depot nodes
    - iv: vehicle index
    - jv: vehicle type index
    - id: depot node index
    - qv: capacity
    - lv: range
    - sv: speed
    - tf: re-fueling time at the depot node
    - td: service time per package at the depot node
    - tc: parking time at a customer node
    - tw: driver working-hours
    - r: driver work-load (maximum vehicle-routes)
    - cd: distance-based operational cost
    - ct: time-based operational cost
    - cf: fixed cost
Note, at least one vehicle of every type at each depot node must be detailed. 
The size of the fleet at the depot nodes is optimized by this tool. 

