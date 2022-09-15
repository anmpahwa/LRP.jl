# Insert node nᵒ between tail node nᵗ and head node nʰ in route rᵒ in solution s.
function insertnode!(nᵒ::Node, nᵗ::Node, nʰ::Node, rᵒ::Route, s::Solution)
    dᵒ =  s.D[rᵒ.iᵈ]
    vᵒ = dᵒ.V[rᵒ.iᵛ]
    tᵉ = rᵒ.tˢ - vᵒ.τᶠ * (rᵒ.l/vᵒ.l) - vᵒ.τᵈ * rᵒ.q
    # update tail node and head node indices
    isdepot(nᵗ) ? rᵒ.iˢ = nᵒ.iⁿ : nᵗ.iʰ = nᵒ.iⁿ
    isdepot(nʰ) ? rᵒ.iᵉ = nᵒ.iⁿ : nʰ.iᵗ = nᵒ.iⁿ
    isdepot(nᵒ) ? (rᵒ.iˢ, rᵒ.iᵉ) = (nʰ.iⁿ, nᵗ.iⁿ) : (nᵒ.iʰ, nᵒ.iᵗ) = (nʰ.iⁿ, nᵗ.iⁿ)
    # update route
    if iscustomer(nᵒ)
        nᵒ.r  = rᵒ
        rᵒ.n += 1
        rᵒ.q += nᵒ.q
    end
    rᵒ.l += s.A[(nᵗ.iⁿ, nᵒ.iⁿ)].l + s.A[(nᵒ.iⁿ, nʰ.iⁿ)].l - s.A[(nᵗ.iⁿ, nʰ.iⁿ)].l
    # update arrival and departure time
    for r ∈ vᵒ.R
        if !isopt(r) continue end
        if r.tˢ < rᵒ.tˢ continue end
        r.tˢ = tᵉ + vᵒ.τᶠ * (r.l/vᵒ.l) + vᵒ.τᵈ * r.q
        cˢ = s.C[r.iˢ]
        cᵉ = s.C[r.iᵉ]   
        tᵈ = r.tˢ
        cᵒ = cˢ
        while true
            cᵒ.tᵃ = tᵈ + s.A[(cᵒ.iᵗ, cᵒ.iⁿ)].l/vᵒ.s
            cᵒ.tᵈ = cᵒ.tᵃ + max(0., cᵒ.tᵉ - cᵒ.tᵃ) + vᵒ.τᶜ
            if isequal(cᵒ, cᵉ) break end
            tᵈ = cᵒ.tᵈ
            cᵒ = s.C[cᵒ.iʰ]
        end
        r.tᵉ = cᵉ.tᵈ + s.A[(cᵉ.iⁿ, dᵒ.iⁿ)].l/vᵒ.s
        tᵉ = r.tᵉ 
    end
    if isempty(vᵒ.R) vᵒ.tˢ, vᵒ.tᵉ = 0., 0.
    else vᵒ.tˢ, vᵒ.tᵉ = vᵒ.R[1].tˢ, vᵒ.R[length(vᵒ.R)].tᵉ
    end
    return s
end

# Remove node nᵒ from its position between tail node nᵗ and head node nʰ in route rᵒ in solution s.
function removenode!(nᵒ::Node, nᵗ::Node, nʰ::Node, rᵒ::Route, s::Solution)
    dᵒ =  s.D[rᵒ.iᵈ]
    vᵒ = dᵒ.V[rᵒ.iᵛ]
    tᵉ = rᵒ.tˢ - vᵒ.τᶠ * (rᵒ.l/vᵒ.l) - vᵒ.τᵈ * rᵒ.q
    # update tail node and head node indices
    isdepot(nᵗ) ? rᵒ.iˢ = nʰ.iⁿ : nᵗ.iʰ = nʰ.iⁿ
    isdepot(nʰ) ? rᵒ.iᵉ = nᵗ.iⁿ : nʰ.iᵗ = nᵗ.iⁿ
    isdepot(nᵒ) ? false : (nᵒ.iʰ, nᵒ.iᵗ) = (0, 0)
    # update route
    if iscustomer(nᵒ)
        nᵒ.r  = NullRoute
        rᵒ.n -= 1
        rᵒ.q -= nᵒ.q
    end
    rᵒ.l -= s.A[(nᵗ.iⁿ, nᵒ.iⁿ)].l + s.A[(nᵒ.iⁿ, nʰ.iⁿ)].l - s.A[(nᵗ.iⁿ, nʰ.iⁿ)].l
    # update arrival and departure time
    for r ∈ vᵒ.R
        if !isopt(r) continue end
        if r.tˢ < rᵒ.tˢ continue end
        r.tˢ = tᵉ + vᵒ.τᶠ * (r.l/vᵒ.l) + vᵒ.τᵈ * r.q
        cˢ = s.C[r.iˢ]
        cᵉ = s.C[r.iᵉ]
        tᵈ = r.tˢ
        cᵒ = cˢ
        while true
            cᵒ.tᵃ = tᵈ + s.A[(cᵒ.iᵗ, cᵒ.iⁿ)].l/vᵒ.s
            cᵒ.tᵈ = cᵒ.tᵃ + max(0., cᵒ.tᵉ - cᵒ.tᵃ) + vᵒ.τᶜ
            if isequal(cᵒ, cᵉ) break end
            tᵈ = cᵒ.tᵈ
            cᵒ = s.C[cᵒ.iʰ]
        end
        r.tᵉ = cᵉ.tᵈ + s.A[(cᵉ.iⁿ, dᵒ.iⁿ)].l/vᵒ.s
        tᵉ = r.tᵉ
    end
    if iscustomer(nᵒ) nᵒ.tᵃ, nᵒ.tᵈ = Inf, Inf end
    if isempty(vᵒ.R) vᵒ.tˢ, vᵒ.tᵉ = 0., 0.
    else vᵒ.tˢ, vᵒ.tᵉ = vᵒ.R[1].tˢ, vᵒ.R[length(vᵒ.R)].tᵉ
    end
    return s
end

# Return true if vehicle vᵒ needs another route (adds conservatively)
function addroute(vᵒ::Vehicle, s::Solution)
    D  = s.D
    dᵒ = D[vᵒ.iᵈ]
    # condtions when route mustn't be added
    if any(!isopt, vᵒ.R) return false end
    for v ∈ dᵒ.V if v.tᵉ > v.w return false end end
    qᵈ = 0
    for v ∈ dᵒ.V for r ∈ v.R qᵈ += r.q end end
    if qᵈ ≥ dᵒ.q return false end
    # condition when route could be added
    if isempty(vᵒ.R) return true end
    for v ∈ dᵒ.V for r ∈ v.R if r.q > v.q return true end end end
    for d ∈ D
        qᵈ = 0
        if isequal(dᵒ, d) continue end
        for v ∈ d.V for r ∈ v.R qᵈ += r.q end end
        if qᵈ > d.q return true end
    end
    return false
end

# Return true if route rᵒ can be deleted (deletes liberally)
function deleteroute(rᵒ::Route, s::Solution)
    # condtions when route mustn't be deleted
    if isopt(rᵒ) return false end
    # condition when route could be deleted
    return true
end

# Return true if depot dᵒ needs another vehicle
function addvehicle(dᵒ::DepotNode, s::Solution)
    D = s.D
    C = s.C
    # condtions when vehicle mustn't be added
    if any(!isopt, dᵒ.V) return false end
    # condition when vehicle could be added
    for c ∈ C
        r = c.r
        if isopen(c) continue end
        d = D[r.iᵈ]
        if !isequal(dᵒ, d) continue end
        if c.tᵃ > c.tˡ return true end 
    end
    for v ∈ dᵒ.V if v.tᵉ > v.w return true end end
    return false
end

# Return false if vehicle vᵒ can be deleted
function deletevehicle(vᵒ::Vehicle, s::Solution)
    D  = s.D
    dᵒ = D[vᵒ.iᵈ]
    # condtions when vehicle mustn't be deleted
    if isopt(vᵒ) return false end
    # condition when vehicle could be deleted
    for v ∈ dᵒ.V
        if isequal(vᵒ, v) continue end
        if isidentical(vᵒ, v) return true end 
    end
    return false
end