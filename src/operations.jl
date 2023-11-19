"""
    insertnode!(nᵒ::Node, nᵗ::Node, nʰ::Node, rᵒ::Route, s::Solution)

Returns solution `s` after inserting node `nᵒ` between tail node `nᵗ` 
and head node `nʰ` in route `rᵒ` in solution `s`.
"""
function insertnode!(nᵒ::Node, nᵗ::Node, nʰ::Node, rᵒ::Route, s::Solution)
    # Note: if nᵒ is a depot node to be inserted between tail node nᵗ and
    # head node nʰ in route rᵒ in solution s, then the route must belong
    # to the depot node nᵒ in the first place
    if isdepot(nᵒ) && !isequal(nᵒ.iⁿ, rᵒ.iᵈ) return s end

    dᵒ =  s.D[rᵒ.iᵈ]
    vᵒ = dᵒ.V[rᵒ.iᵛ]
    tⁱ = rᵒ.tⁱ
    θⁱ = rᵒ.θⁱ
    aᵒ = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
    aᵗ = s.A[(nᵗ.iⁿ, nᵒ.iⁿ)]
    aʰ = s.A[(nᵒ.iⁿ, nʰ.iⁿ)]

    # update associated customer nodes, route, vehicle, and depot node
    if iscustomer(nᵒ)
        nᵒ.iʳ = rᵒ.iʳ
        nᵒ.iᵛ = rᵒ.iᵛ
        nᵒ.iᵈ = rᵒ.iᵈ
        isdepot(nᵗ) ? rᵒ.iˢ = nᵒ.iⁿ : nᵗ.iʰ = nᵒ.iⁿ
        isdepot(nʰ) ? rᵒ.iᵉ = nᵒ.iⁿ : nʰ.iᵗ = nᵒ.iⁿ
        nᵒ.iʰ = nʰ.iⁿ
        nᵒ.iᵗ = nᵗ.iⁿ
        nᵒ.r = rᵒ

        rᵒ.x = (rᵒ.n * rᵒ.x + nᵒ.x)/(rᵒ.n + 1)
        rᵒ.y = (rᵒ.n * rᵒ.y + nᵒ.y)/(rᵒ.n + 1)
        rᵒ.n += 1
        rᵒ.q += nᵒ.q
        rᵒ.l += aᵗ.l + aʰ.l - aᵒ.l

        vᵒ.n += 1
        vᵒ.q += nᵒ.q
        vᵒ.l += aᵗ.l + aʰ.l - aᵒ.l

        dᵒ.n += 1
        dᵒ.q += nᵒ.q
        dᵒ.l += aᵗ.l + aʰ.l - aᵒ.l
    end
    # Note: In the context of the above note, on removal of a depot node
    # the associated, routes', vehicles', and customers' hierarchical indices 
    # are not changed, and therefore on insertion of a depot node, the
    # corresponding changes are not required.
    if isdepot(nᵒ)
        nᵗ.iʰ = nᵒ.iⁿ
        nʰ.iᵗ = nᵒ.iⁿ

        rᵒ.iˢ = nʰ.iⁿ
        rᵒ.iᵉ = nᵗ.iⁿ
        rᵒ.l += aᵗ.l + aʰ.l - aᵒ.l

        vᵒ.l += aᵗ.l + aʰ.l - aᵒ.l

        dᵒ.l += aᵗ.l + aʰ.l - aᵒ.l
    end

    # update arrival and departure time
    if isequal(φᵀ::Bool, true)
        for r ∈ vᵒ.R
            if r.tⁱ < rᵒ.tⁱ continue end
            if isopt(r)
                r.θⁱ = θⁱ
                r.θˢ = θⁱ + max(0., (r.l/vᵒ.lᵛ - r.θⁱ))
                r.tⁱ = tⁱ
                r.tˢ = r.tⁱ + vᵒ.τᶠ * (r.θˢ - r.θⁱ) + vᵒ.τᵈ * r.q
                cˢ = s.C[r.iˢ]
                cᵉ = s.C[r.iᵉ]
                tᵈ = r.tˢ
                cᵒ = cˢ
                while true
                    cᵒ.tᵃ = tᵈ + s.A[(cᵒ.iᵗ, cᵒ.iⁿ)].l/vᵒ.sᵛ
                    cᵒ.tᵈ = cᵒ.tᵃ + vᵒ.τᶜ + max(0., cᵒ.tᵉ - cᵒ.tᵃ - vᵒ.τᶜ) + cᵒ.τᶜ
                    if isequal(cᵒ, cᵉ) break end
                    tᵈ = cᵒ.tᵈ
                    cᵒ = s.C[cᵒ.iʰ]
                end
                r.θᵉ = r.θˢ - r.l/vᵒ.lᵛ
                r.tᵉ = cᵉ.tᵈ + s.A[(cᵉ.iⁿ, dᵒ.iⁿ)].l/vᵒ.sᵛ
            else
                r.θⁱ = θⁱ
                r.θˢ = θⁱ
                r.θᵉ = θⁱ
                r.tⁱ = tⁱ
                r.tˢ = tⁱ
                r.tᵉ = tⁱ
            end
            tⁱ = r.tᵉ
            θⁱ = r.θᵉ
        end

        # update start and end time
        (vᵒ.tˢ, vᵒ.tᵉ) = isempty(vᵒ.R) ? (dᵒ.tˢ, dᵒ.tˢ) : (vᵒ.R[firstindex(vᵒ.R)].tⁱ, vᵒ.R[lastindex(vᵒ.R)].tᵉ)

        # update slack
        τ = dᵒ.tᵉ - vᵒ.tᵉ
        for r ∈ reverse(vᵒ.R)
            if !isopt(r) continue end
            cˢ = s.C[r.iˢ]
            cᵉ = s.C[r.iᵉ]
            cᵒ = cˢ
            while true
                τ  = min(τ, cᵒ.tˡ - cᵒ.tᵃ - vᵒ.τᶜ)
                if isequal(cᵒ, cᵉ) break end
                cᵒ = s.C[cᵒ.iʰ]
            end
            r.τ = τ
        end
        vᵒ.τ = min(τ, vᵒ.τ)
        dᵒ.τ = min(τ, dᵒ.τ)
    end

    return s
end
"""
    removenode!(nᵒ::Node, nᵗ::Node, nʰ::Node, rᵒ::Route, s::Solution)

Returns solution `s` after removing node `nᵒ` from its position between 
tail node `nᵗ` and head node `nʰ` in route `rᵒ` in solution `s`.
"""
function removenode!(nᵒ::Node, nᵗ::Node, nʰ::Node, rᵒ::Route, s::Solution)
    # Note: if nᵒ is a depot node to be removed from its position between
    # tail node nᵗ and head node nʰ in route rᵒ in solution s, then the 
    # route must belong to the depot node nᵒ in the first place
    if isdepot(nᵒ) && !isequal(nᵒ.iⁿ, rᵒ.iᵈ) return s end

    dᵒ =  s.D[rᵒ.iᵈ]
    vᵒ = dᵒ.V[rᵒ.iᵛ]
    tⁱ = rᵒ.tⁱ
    θⁱ = rᵒ.θⁱ
    aᵒ = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
    aᵗ = s.A[(nᵗ.iⁿ, nᵒ.iⁿ)]
    aʰ = s.A[(nᵒ.iⁿ, nʰ.iⁿ)]

    # update associated customer nodes, route, vehicle, and depot node
    if iscustomer(nᵒ)
        nᵒ.iʳ = 0
        nᵒ.iᵛ = 0
        nᵒ.iᵈ = 0
        isdepot(nᵗ) ? rᵒ.iˢ = nʰ.iⁿ : nᵗ.iʰ = nʰ.iⁿ
        isdepot(nʰ) ? rᵒ.iᵉ = nᵗ.iⁿ : nʰ.iᵗ = nᵗ.iⁿ
        nᵒ.iʰ = 0
        nᵒ.iᵗ = 0
        nᵒ.r = NullRoute

        rᵒ.x = isone(rᵒ.n) ? 0. : (rᵒ.n * rᵒ.x - nᵒ.x)/(rᵒ.n - 1)
        rᵒ.y = isone(rᵒ.n) ? 0. : (rᵒ.n * rᵒ.y - nᵒ.y)/(rᵒ.n - 1)
        rᵒ.n -= 1
        rᵒ.q -= nᵒ.q
        rᵒ.l -= aᵗ.l + aʰ.l - aᵒ.l

        vᵒ.n -= 1
        vᵒ.q -= nᵒ.q
        vᵒ.l -= aᵗ.l + aʰ.l - aᵒ.l

        dᵒ.n -= 1
        dᵒ.q -= nᵒ.q
        dᵒ.l -= aᵗ.l + aʰ.l - aᵒ.l
    end
    # Note: In the context of the above note, on removal of a depot node
    # the associated, routes', vehicles', and customers' hierarchical indices 
    # are not changed
    if isdepot(nᵒ)
        nᵗ.iʰ = nʰ.iⁿ
        nʰ.iᵗ = nᵗ.iⁿ

        rᵒ.iˢ = nʰ.iⁿ
        rᵒ.iᵉ = nᵗ.iⁿ
        rᵒ.l -= aᵗ.l + aʰ.l - aᵒ.l

        vᵒ.l -= aᵗ.l + aʰ.l - aᵒ.l

        dᵒ.l -= aᵗ.l + aʰ.l - aᵒ.l
    end

    # update arrival and departure time
    if isequal(φᵀ::Bool, true)
        if iscustomer(nᵒ) nᵒ.tᵃ, nᵒ.tᵈ = 0., 0. end
        for r ∈ vᵒ.R
            if r.tⁱ < rᵒ.tⁱ continue end
            if isopt(r)
                r.θⁱ = θⁱ
                r.θˢ = θⁱ + max(0., (r.l/vᵒ.lᵛ - r.θⁱ))
                r.tⁱ = tⁱ
                r.tˢ = r.tⁱ + vᵒ.τᶠ * (r.θˢ - r.θⁱ) + vᵒ.τᵈ * r.q
                cˢ = s.C[r.iˢ]
                cᵉ = s.C[r.iᵉ]
                tᵈ = r.tˢ
                cᵒ = cˢ
                while true
                    cᵒ.tᵃ = tᵈ + s.A[(cᵒ.iᵗ, cᵒ.iⁿ)].l/vᵒ.sᵛ
                    cᵒ.tᵈ = cᵒ.tᵃ + vᵒ.τᶜ + max(0., cᵒ.tᵉ - cᵒ.tᵃ - vᵒ.τᶜ) + cᵒ.τᶜ
                    if isequal(cᵒ, cᵉ) break end
                    tᵈ = cᵒ.tᵈ
                    cᵒ = s.C[cᵒ.iʰ]
                end
                r.θᵉ = r.θˢ - r.l/vᵒ.lᵛ
                r.tᵉ = cᵉ.tᵈ + s.A[(cᵉ.iⁿ, dᵒ.iⁿ)].l/vᵒ.sᵛ
            else
                r.θⁱ = θⁱ
                r.θˢ = θⁱ
                r.θᵉ = θⁱ
                r.tⁱ = tⁱ
                r.tˢ = tⁱ
                r.tᵉ = tⁱ
            end
            tⁱ = r.tᵉ
            θⁱ = r.θᵉ
        end

        # update start and end time
        (vᵒ.tˢ, vᵒ.tᵉ) = isempty(vᵒ.R) ? (dᵒ.tˢ, dᵒ.tˢ) : (vᵒ.R[firstindex(vᵒ.R)].tⁱ, vᵒ.R[lastindex(vᵒ.R)].tᵉ)

        # update slack
        τ = dᵒ.tᵉ - vᵒ.tᵉ
        for r ∈ reverse(vᵒ.R)
            if !isopt(r) continue end
            cˢ = s.C[r.iˢ]
            cᵉ = s.C[r.iᵉ]
            cᵒ = cˢ
            while true
                τ  = min(τ, cᵒ.tˡ - cᵒ.tᵃ - vᵒ.τᶜ)
                if isequal(cᵒ, cᵉ) break end
                cᵒ = s.C[cᵒ.iʰ]
            end
            r.τ = τ
        end
        vᵒ.τ = min(τ, vᵒ.τ)
        dᵒ.τ = min(τ, dᵒ.τ)
    end

    return s
end



"""
    addroute(rᵒ::Route, s::Solution)

Returns `true` if a route `rᵒ` can be added into solution `s` (conservative).
"""
function addroute(rᵒ::Route, s::Solution)
    dᵒ =  s.D[rᵒ.iᵈ]
    vᵒ = dᵒ.V[rᵒ.iᵛ]

    # condtions when route mustn't be added
    if isequal(length(vᵒ.R), vᵒ.r̅) return false end
    if any(!isopt, vᵒ.R) return false end
    if vᵒ.tᵉ > dᵒ.tᵉ return false end
    if (vᵒ.tᵉ - vᵒ.tˢ) > vᵒ.τʷ return false end
    if dᵒ.q ≥ dᵒ.qᵈ return false end

    # condition when route could be added
    if isempty(vᵒ.R) return true end
    for v ∈ dᵒ.V for r ∈ v.R if r.q > v.qᵛ return true end end end
    for d ∈ s.D if dᵒ.q > d.qᵈ return true end end

    return false
end
"""
    deleteroute(rᵒ::Route, s::Solution)

Returns `true` if route `rᵒ` can be deleted from solution `s` (liberal).
"""
function deleteroute(rᵒ::Route, s::Solution)
    # condtions when route mustn't be deleted
    if isopt(rᵒ) return false end

    # condition when route could be deleted


    return true
end



"""
    addvehicle(vᵒ::Vehicle, s::Solution)

Returns `true` if vehicle `vᵒ` can be added into solution `s` (conservative).
"""
function addvehicle(vᵒ::Vehicle, s::Solution)
    dᵒ = s.D[vᵒ.iᵈ]

    # condtions when vehicle mustn't be added
    if any(!isopt, filter(v -> isidentical(vᵒ, v), dᵒ.V)) return false end
    if dᵒ.q ≥ dᵒ.qᵈ return false end

    # condition when vehicle could be added
    if dᵒ.q < dᵒ.qᵈ return true end
    for v ∈ dᵒ.V
        if v.tᵉ > dᵒ.tᵉ return true end
        if (v.tᵉ - v.tˢ) > v.τʷ return true end
        for r ∈ v.R
            if !isopt(r) continue end
            cˢ = s.C[r.iˢ]
            cᵉ = s.C[r.iᵉ]
            cᵒ = cˢ
            while true
                if cᵒ.tᵃ > cᵒ.tˡ return true end
                if isequal(cᵒ, cᵉ) break end
                cᵒ = s.C[cᵒ.iʰ]
            end
        end
    end

    return false
end
"""
    deletevehicle(vᵒ::Vehicle, s::Solution)

Returns `true` if vehicle `vᵒ` can be deleted from solution `s` (liberal).
"""
function deletevehicle(vᵒ::Vehicle, s::Solution)
    dᵒ = s.D[vᵒ.iᵈ]

    # condtions when vehicle mustn't be deleted
    if isopt(vᵒ) return false end
    nʲ = 0
    for v ∈ dᵒ.V if isidentical(vᵒ, v) nʲ += 1 end end
    if isone(nʲ) return false end
    
    # condition when vehicle could be deleted


    return true
end



"""
    preinitialize!(s::Solution)

Pre-intialization procedures.
Returns solution `s` after adding new routes into the solution.
"""
function preinitialize!(s::Solution)
    for d ∈ s.D
        for v ∈ d.V
            rᵒ = Route(v, d)
            if addroute(rᵒ, s) push!(v.R, rᵒ) end
        end
    end
    return s
end
"""
    postnitialize!(s::Solution)

Post-intialization procedures.
Returns solution `s`.
"""
function postinitialize!(s::Solution)
    return s
end



"""
    preremove!(s::Solution)

Pre-removal procedures. 
Returns solution `s`.
"""
function preremove!(s::Solution)
    return s
end
"""
    postremove!(s::Solution)

Post-removal procedures.
Returns solution `s`.
"""
function postremove!(s::Solution)
    return s
end



"""
    preinsert!(s::Solution)

Pre-insertion procedures.
Returns solution `s` after adding new vehicles and routes into the solution.
"""
function preinsert!(s::Solution)
    for d ∈ s.D
        for v ∈ d.V
            rᵒ = Route(v, d)
            if addroute(rᵒ, s) push!(v.R, rᵒ) end
            vᵒ = Vehicle(v, d)
            if addvehicle(vᵒ, s) push!(d.V, vᵒ) end
        end
    end
    return s
end
"""
    postinsert!(s::Solution)

Post-insertion procedures. 
Retruns solution `s` after deleting routes and vehicles, and subsequently correcting vehicle and route indices.
"""
function postinsert!(s::Solution)
    for d ∈ s.D
        k = 1
        while true
            v = d.V[k]
            if deletevehicle(v, s) 
                deleteat!(d.V, k)
            else
                v.iᵛ = k
                for r ∈ v.R r.iᵛ = k end
                k += 1
            end
            if k > length(d.V) break end
        end
        for v ∈ d.V
            if isempty(v.R) continue end
            k = 1
            while true
                r = v.R[k]
                if deleteroute(r, s) 
                    deleteat!(v.R, k)
                else
                    r.iʳ = k
                    k += 1
                end
                if k > length(v.R) break end
            end
        end
    end
    for c ∈ s.C 
        c.iʳ = c.r.iʳ
        c.iᵛ = c.r.iᵛ
        c.iᵈ = c.r.iᵈ
    end
    return s
end



"""
    prelocalsearch!(s::Solution)

Pre-localsearch procedures. 
Returns solution `s`.
"""
function prelocalsearch!(s::Solution)
    return s
end
"""
    postlocalsearch!(s::Solution)

Post-localsearch procedures. 
Retruns solution `s`.
"""
function postlocalsearch!(s::Solution)
    return s
end