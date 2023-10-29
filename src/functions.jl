"""
    isopt(r::Route)

Returns true if route `r` is operational.
A route is defined operational if it serves at least one customer.
"""
isopt(r::Route) = (r.n ≥ 1)
"""
    isopt(v::Vehicle)

Returns true if vehicle `v` is operational.
A vehicle is defined operational if any of its routes is operational.
"""
isopt(v::Vehicle) = any(isopt, v.R)
"""
    isopt(d::DepotNode)
    
Returns true if depot node `d` is operational.
A depot node is defined operational if any of its vehicles is operational.
"""
isopt(d::DepotNode) = any(isopt, d.V)



"""
    isopen(c::CustomerNode)
    
Returns true if customer node `c` is open.
A customer node is defined open if it is not being served by any vehicle-route.
"""
isopen(c::CustomerNode) = isequal(c.r, NullRoute)



"""
    isclose(d::DepotNode)

Returns true if depot node `d` is not operational.
A depot node is defined operational if any of its vehicles is operational.
"""  
isclose(d::DepotNode)    = !isopt(d)
"""
    isclose(c::CustomerNode)

Returns true if customer node `c` is not open.
A customer node is defined open if it is not being served by any vehicle-route.
"""  
isclose(c::CustomerNode) = !isopen(c)



"""
    isequal(p::Route, q::Route)

Return true if route `p` equals route `q`.
Two routes are the equal if their indices (`iᵈ`, `iᵛ`, `iʳ`) match.
"""
Base.isequal(p::Route, q::Route)     = isequal(p.iᵈ, q.iᵈ) && isequal(p.iᵛ, q.iᵛ) && isequal(p.iʳ, q.iʳ)
"""
    isequal(p::Vehicle, q::Vehicle)

Return true if vehicle `p` equals vehicle `q`.
Two vehicles are equal if their indices (`iᵈ`, `iᵛ`) match.
"""
Base.isequal(p::Vehicle, q::Vehicle) = isequal(p.iᵈ, q.iᵈ) && isequal(p.iᵛ, q.iᵛ)
"""
    isequal(p::Node, q::Node)

Return true if node `p` equals node `q`.
Two node are equal if their indices (`iⁿ`) match.
"""
Base.isequal(p::Node, q::Node)       = isequal(p.iⁿ, q.iⁿ)

"""
    isidentical(p::Vehicle, q::Vehicle)

Return true if vehicle `p` is identical to vehicle `q`.
Two vehicles are identical if they are of the same type,
i.e. if their type indices (`jᵛ`) match.
"""
isidentical(p::Vehicle, q::Vehicle) = isequal(p.jᵛ, q.jᵛ)

"""
    isdepot(n::Node)

Returns true if node `n` is a depot.
"""
isdepot(n::Node)    = isequal(typeof(n), DepotNode)
"""
    iscustomer(n::Node)
    
Returns true if node `n` is a customer.
"""
iscustomer(n::Node) = isequal(typeof(n), CustomerNode)



"""
    Route(v::Vehicle, d::DepotNode)

Returns a non-operational route traversed by vehicle `v` from depot node `d`.
"""
function Route(v::Vehicle, d::DepotNode)
    iʳ = length(v.R) + 1
    iᵛ = v.iᵛ
    iᵈ = d.iⁿ
    x  = 0.
    y  = 0. 
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
    q  = 0.
    l  = 0.
    r  = Route(iʳ, iᵛ, iᵈ, x, y, iˢ, iᵉ, θⁱ, θˢ, θᵉ, tⁱ, tˢ, tᵉ, τ, n, q, l)
    return r
end            
const NullRoute = Route(0, 0, 0, 0., 0., 0, 0, 0., 0., 0., Inf, Inf, Inf, 0., 0, 0, Inf)



"""
    Route(v::Vehicle, d::DepotNode)

    Returns a non-operational vehicle cloning vehicle `v` at depot node `d`.
"""
function Vehicle(v::Vehicle, d::DepotNode)
    iᵛ = length(d.V) + 1
    jᵛ = v.jᵛ
    iᵈ = v.iᵈ
    qᵛ = v.q
    lᵛ = v.l
    s  = v.s
    τᶠ = v.τᶠ
    τᵈ = v.τᵈ
    τᶜ = v.τᶜ
    τʷ = v.τʷ
    r̅  = v.r̅
    tˢ = d.tˢ
    tᵉ = d.tˢ
    πᵈ = v.πᵈ
    πᵗ = v.πᵗ
    πᶠ = v.πᶠ
    R  = Route[]
    v  = Vehicle(iᵛ, jᵛ, iᵈ, qᵛ, lᵛ, s, τᶠ, τᵈ, τᶜ, τʷ, r̅, tˢ, tᵉ, πᵈ, πᵗ, πᶠ, R)
    return v
end



"""
    vectorize(s::Solution)

Returns solution as a sequence of nodes in the order of visits.
"""
function vectorize(s::Solution)
    D = s.D
    C = s.C
    Z = [Int64[] for _ ∈ D]
    for d ∈ D
        iⁿ = d.iⁿ
        if !isopt(d) continue end
        for v ∈ d.V
            if !isopt(v) continue end
            for r ∈ v.R
                if !isopt(r) continue end
                cˢ, cᵉ = C[r.iˢ], C[r.iᵉ] 
                push!(Z[iⁿ], d.iⁿ)
                c = cˢ
                while true
                    push!(Z[iⁿ], c.iⁿ)
                    if isequal(c, cᵉ) break end
                    c = C[c.iʰ]
                end
            end
        end
        push!(Z[iⁿ], d.iⁿ)
    end
    return Z
end
"""
    hash(s::Solution)

Returns hash on vectorized solution.
"""
Base.hash(s::Solution) = hash(vectorize(s))



"""
    f(s::Solution; fixed=true, operational=true, penalty=true)

Returns objective function evaluation for solution `s`. Include `fixed`, 
`operational`, and `penalty` cost for constriant violation if `true`.
"""
function f(s::Solution; fixed=true, operational=true, penalty=true)
    πᶠ, πᵒ, πᵖ = 0., 0., 0.
    φᶠ, φᵒ, φᵖ = fixed, operational, penalty
    for d ∈ s.D
        πᶠ += isopt(d) * d.πᶠ
        qᵈ = 0
        nᵈ = 0
        for v ∈ d.V
            πᶠ += isopt(v) * v.πᶠ
            for r ∈ v.R 
                if !isopt(r) continue end
                qᵛ = r.q
                lᵛ = r.l
                qᵈ += qᵛ
                nᵈ += r.n
                πᵒ += r.l * v.πᵈ
                πᵖ += (qᵛ > v.q) * (qᵛ - v.q)                               # Vehicle capacity constraint
                πᵖ += (lᵛ > v.l) * (lᵛ - v.l)                               # Vehicle range constraint
            end
            πᵒ += (v.tᵉ - v.tˢ) * v.πᵗ
            πᵖ += (d.tˢ > v.tˢ) * (d.tˢ - v.tˢ)                             # Working-hours constraint (start time)
            πᵖ += (v.tᵉ > d.tᵉ) * (v.tᵉ - d.tᵉ)                             # Working-hours constraint (end time)
            πᵖ += (v.tᵉ - v.tˢ > v.τʷ) * (v.tᵉ - v.tˢ - v.τʷ)               # Working-hours constraint (duration)
        end
        pᵈ  = nᵈ/length(s.C)
        πᵒ += qᵈ * d.πᵒ
        πᵖ += (isone(d.φ) && !isopt(d)) * d.πᶠ                              # Depot use constraint
        πᵖ += (qᵈ > d.q) * (qᵈ - d.q)                                       # Depot capacity constraint
        πᵖ += (pᵈ < d.pˡ) * (d.pˡ - pᵈ)                                     # Depot customer share constraint
        πᵖ += (pᵈ > d.pᵘ) * (pᵈ - d.pᵘ)                                     # Depot customer share constraint
    end
    for c ∈ s.C πᵖ += isopen(c) ? 0. : (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) end    # Time-window constraint
    z = φᶠ * πᶠ + φᵒ * πᵒ + φᵖ * πᵖ * 10^(ceil(log10(πᶠ + πᵒ)))
    return z
end



"""
    isfeasible(s::Solution)

Returns true if node service, node flow, and sub-tour elimination
constraints; depot and vehicle capacity constriants; vehicle range 
and working-hours constraints; and time-window constraints are not 
violated.
"""
function isfeasible(s::Solution)
    X = zeros(Int64, eachindex(s.C))
    for d ∈ s.D
        qᵈ = 0
        nᵈ = 0.
        for v ∈ d.V
            for r ∈ v.R
                if !isopt(r) continue end
                qᵛ = r.q
                lᵛ = r.l
                if qᵛ > v.q return false end                                # Vehicle capacity constraint
                if lᵛ > v.l return false end                                # Vehicle range constraint
                qᵈ += r.q
                nᵈ += r.n
                cˢ = s.C[r.iˢ]
                cᵉ = s.C[r.iᵉ]
                cᵒ = cˢ
                while true
                    if cᵒ.tᵃ > cᵒ.tˡ return false end                       # Time-window constraint
                    X[cᵒ.iⁿ] += 1
                    if isequal(cᵒ, cᵉ) break end
                    cᵒ = s.C[cᵒ.iʰ]
                end
            end
            if d.tˢ > v.tˢ return false end                                 # Working-hours constraint (start time)
            if v.tᵉ > d.tᵉ return false end                                 # Working-hours constraint (end time)
            if v.tᵉ - v.tˢ > v.τʷ return false end                          # Working-hours constraint (duration)
        end
        pᵈ = nᵈ/length(s.C)
        if (isone(d.φ) && !isopt(d)) return false end                       # Depot use constraint
        if qᵈ > d.q return false end                                        # Depot capacity constraint
        if !(d.pˡ ≤ pᵈ ≤ d.pᵘ) return false end                             # Depot customer share constraint
    end
    if any(!isone, X) return false end                                      # Node service, customer flow, and sub-tour elimination constrinat
    return true
end



"""
    relatedness(c::CustomerNode, d::DepotNode, s::Solution)

Returns a measure of similarity between customer nodes `c` and depot node `d`.
"""
function relatedness(c::CustomerNode, d::DepotNode, s::Solution)
    l = s.A[(c.iⁿ,d.iⁿ)].l
    t = 0
    q = 0
    r = c.r
    φ = 1 + isequal(r.iᵈ, d.iⁿ)
    z = (q + φ)/(l + t)
    return z
end
"""
    relatedness(d::DepotNode, c::CustomerNode, s::Solution)

Returns a measure of similarity between depot node `d` and customer nodes `c`.
"""
relatedness(d::DepotNode, c::CustomerNode, s::Solution) = relatedness(c, d, s)
"""
    relatedness(c¹::CustomerNode, c²::CustomerNode, s::Solution)

Returns a measure of similarity between customer nodes `c¹` and `c²`.
"""
function relatedness(c¹::CustomerNode, c²::CustomerNode, s::Solution)
    if isequal(c¹, c²) return Inf end
    r¹ = c¹.r
    r² = c².r
    d¹ = s.D[r¹.iᵈ]
    d² = s.D[r².iᵈ]
    v¹ = d¹.V[r¹.iᵛ]
    v² = d².V[r².iᵛ]
    φᵈ = isequal(d¹, d²)
    φᵛ = isequal(v¹, v²)
    φʳ = isequal(r¹, r²)
    φ  = 1 + φᵈ + φᵛ + φʳ  
    l  = s.A[(c¹.iⁿ,c².iⁿ)].l
    t  = abs(c¹.tᵉ - c².tᵉ) + abs(c¹.tˡ - c².tˡ)
    q  = abs(c¹.q - c².q)
    z  = (q + φ)/(l + t)
    return z
end
"""
    relatedness(r¹::Route, r²::Route, s::Solution)

Returns a measure of similarity between routes `r¹` and `r²`.
"""
function relatedness(r¹::Route, r²::Route, s::Solution)
    if !isopt(r¹) || !isopt(r²) return -Inf end
    if isequal(r¹, r²) return Inf end
    d¹ = s.D[r¹.iᵈ]
    d² = s.D[r².iᵈ]
    v¹ = d¹.V[r¹.iᵛ]
    v² = d².V[r².iᵛ]
    φᵈ = isequal(d¹, d²)
    φᵛ = isequal(v¹, v²)
    φ  = 1 + φᵈ + φᵛ 
    l  = sqrt((r¹.x - r².x)^2 + (r¹.y - r².y)^2)
    t  = abs(r¹.tˢ - r².tˢ) + abs(r¹.tᵉ - r².tᵉ)
    q  = abs(r¹.q - r².q)
    z  = (q + φ)/(l + t)
    return z
end
"""
    relatedness(v¹::Vehicle, v²::Vehicle, s::Solution)

Returns a measure of similarity between vehicles `v¹` and `v²`.
"""
function relatedness(v¹::Vehicle, v²::Vehicle, s::Solution)
    if !isopt(v¹) || !isopt(v²) return -Inf end
    if isequal(v¹, v²) return Inf end
    d¹ = s.D[v¹.iᵈ]
    d² = s.D[v².iᵈ]
    φᵈ = isequal(d¹, d²)
    a  = 0.
    b  = 0.
    c  = 0
    for r ∈ v¹.R 
        a += r.n * r.x
        b += r.n * r.y
        c += r.n
    end
    x¹ = a/c
    y¹ = b/c 
    a  = 0.
    b  = 0.
    c  = 0
    for r ∈ v².R 
        a += r.n * r.x
        b += r.n * r.y
        c += r.n
    end
    x² = a/c
    y² = b/c 
    q¹ = 0
    q² = 0
    for r ∈ v¹.R q¹ += r.q end
    for r ∈ v².R q² += r.q end
    φ  = 1 + φᵈ
    l  = sqrt((x¹ - x²)^2 + (y¹ - y²)^2)
    t  = abs(v¹.tˢ - v².tˢ) + abs(v¹.tᵉ - v².tᵉ)
    q  = abs(q¹ - q²)
    z  = (q + φ)/(l + t)
    return z
end
