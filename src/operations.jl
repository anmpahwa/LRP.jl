"""
    insertnode!(n::Node, nᵗ::Node, nʰ::Node, r::Route, s::Solution)

Returns solution `s` after inserting node `n` between tail node `nᵗ` 
and head node `nʰ` in route `r` in solution `s`.
"""
function insertnode!(n::Node, nᵗ::Node, nʰ::Node, r::Route, s::Solution)
    if isclose(n) return s end
    # Note
    # If n is a depot node to be inserted between tail node nᵗ and
    # head node nʰ in route r in solution s, then the route must belong
    # to the depot node n in the first place
    if isdepot(n) && !isequal(n.iⁿ, r.iᵈ) return s end
    d  = s.D[r.iᵈ]
    v  = d.V[r.iᵛ]
    tⁱ = r.tⁱ
    θⁱ = r.θⁱ
    aᵒ = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
    aᵗ = s.A[(nᵗ.iⁿ, n.iⁿ)]
    aʰ = s.A[(n.iⁿ, nʰ.iⁿ)]
    # update associated customer nodes, route, vehicle, and depot node
    if iscustomer(n)
        n.iʳ = r.iʳ
        n.iᵛ = r.iᵛ
        n.iᵈ = r.iᵈ
        isdepot(nᵗ) ? r.iˢ = n.iⁿ : nᵗ.iʰ = n.iⁿ
        isdepot(nʰ) ? r.iᵉ = n.iⁿ : nʰ.iᵗ = n.iⁿ
        n.iʰ = nʰ.iⁿ
        n.iᵗ = nᵗ.iⁿ
        n.r  = r

        r.x  = (r.n * r.x + n.x)/(r.n + 1)
        r.y  = (r.n * r.y + n.y)/(r.n + 1)
        r.n += 1
        r.q += n.q
        r.l += aᵗ.l + aʰ.l - aᵒ.l

        v.n += 1
        v.q += n.q
        v.l += aᵗ.l + aʰ.l - aᵒ.l

        d.n += 1
        d.q += n.q
        d.l += aᵗ.l + aʰ.l - aᵒ.l
    end
    # Note
    # In the context of the above note, on removal of a depot node
    # the associated, routes', vehicles', and customers' hierarchical indices 
    # are not changed, and therefore on insertion of a depot node, the
    # corresponding changes are not required.
    if isdepot(n)
        nᵗ.iʰ = n.iⁿ
        nʰ.iᵗ = n.iⁿ

        r.iˢ = nʰ.iⁿ
        r.iᵉ = nᵗ.iⁿ
        r.l += aᵗ.l + aʰ.l - aᵒ.l

        v.l += aᵗ.l + aʰ.l - aᵒ.l

        d.l += aᵗ.l + aʰ.l - aᵒ.l
    end
    # update arrival and departure time
    rᵒ = r
    if isequal(φᵀ::Bool, true)
        for r ∈ v.R
            if r.tⁱ < rᵒ.tⁱ continue end
            if isopt(r)
                r.θⁱ = θⁱ
                r.θˢ = θⁱ + max(0., (r.l/v.lᵛ - r.θⁱ))
                r.tⁱ = tⁱ
                r.tˢ = r.tⁱ + v.τᶠ * (r.θˢ - r.θⁱ) + v.τᵈ * r.q
                cˢ = s.C[r.iˢ]
                cᵉ = s.C[r.iᵉ]
                tᵈ = r.tˢ
                cᵒ = cˢ
                while true
                    cᵒ.tᵃ = tᵈ + s.A[(cᵒ.iᵗ, cᵒ.iⁿ)].l/v.sᵛ
                    cᵒ.tᵈ = cᵒ.tᵃ + v.τᶜ + max(0., cᵒ.tᵉ - cᵒ.tᵃ - v.τᶜ) + cᵒ.τᶜ
                    if isequal(cᵒ, cᵉ) break end
                    tᵈ = cᵒ.tᵈ
                    cᵒ = s.C[cᵒ.iʰ]
                end
                r.θᵉ = r.θˢ - r.l/v.lᵛ
                r.tᵉ = cᵉ.tᵈ + s.A[(cᵉ.iⁿ, d.iⁿ)].l/v.sᵛ
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
        (v.tˢ, v.tᵉ) = isempty(v.R) ? (d.tˢ, d.tˢ) : (v.R[firstindex(v.R)].tⁱ, v.R[lastindex(v.R)].tᵉ)
        # update slack
        τ = d.tᵉ - v.tᵉ
        for r ∈ reverse(v.R)
            if !isopt(r) continue end
            cˢ = s.C[r.iˢ]
            cᵉ = s.C[r.iᵉ]
            cᵒ = cˢ
            while true
                τ  = min(τ, cᵒ.tˡ - cᵒ.tᵃ - v.τᶜ)
                if isequal(cᵒ, cᵉ) break end
                cᵒ = s.C[cᵒ.iʰ]
            end
            r.τ = τ
        end
        v.τ = min(τ, v.τ)
        d.τ = min(τ, d.τ)
    end
    return s
end
"""
    removenode!(n::Node, nᵗ::Node, nʰ::Node, r::Route, s::Solution)

Returns solution `s` after removing node `n` from its position between 
tail node `nᵗ` and head node `nʰ` in route `r` in solution `s`.
"""
function removenode!(n::Node, nᵗ::Node, nʰ::Node, r::Route, s::Solution)
    if isopen(n) return s end
    # Note
    # If n is a depot node to be removed from its position between
    # tail node nᵗ and head node nʰ in route r in solution s, then the 
    # route must belong to the depot node n in the first place
    if isdepot(n) && !isequal(n.iⁿ, r.iᵈ) return s end
    d = s.D[r.iᵈ]
    v = d.V[r.iᵛ]
    tⁱ = r.tⁱ
    θⁱ = r.θⁱ
    aᵒ = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
    aᵗ = s.A[(nᵗ.iⁿ, n.iⁿ)]
    aʰ = s.A[(n.iⁿ, nʰ.iⁿ)]
    # update associated customer nodes, route, vehicle, and depot node
    if iscustomer(n)
        n.iʳ = 0
        n.iᵛ = 0
        n.iᵈ = 0
        isdepot(nᵗ) ? r.iˢ = nʰ.iⁿ : nᵗ.iʰ = nʰ.iⁿ
        isdepot(nʰ) ? r.iᵉ = nᵗ.iⁿ : nʰ.iᵗ = nᵗ.iⁿ
        n.iʰ = 0
        n.iᵗ = 0
        n.r  = NullRoute

        r.x  = isone(r.n) ? 0. : (r.n * r.x - n.x)/(r.n - 1)
        r.y  = isone(r.n) ? 0. : (r.n * r.y - n.y)/(r.n - 1)
        r.n -= 1
        r.q -= n.q
        r.l -= aᵗ.l + aʰ.l - aᵒ.l

        v.n -= 1
        v.q -= n.q
        v.l -= aᵗ.l + aʰ.l - aᵒ.l

        d.n -= 1
        d.q -= n.q
        d.l -= aᵗ.l + aʰ.l - aᵒ.l
    end
    # Note
    # In the context of the above note, on removal of a depot node
    # the associated, routes', vehicles', and customers' hierarchical indices 
    # are not changed
    if isdepot(n)
        nᵗ.iʰ = nʰ.iⁿ
        nʰ.iᵗ = nᵗ.iⁿ

        r.iˢ = nʰ.iⁿ
        r.iᵉ = nᵗ.iⁿ
        r.l -= aᵗ.l + aʰ.l - aᵒ.l

        v.l -= aᵗ.l + aʰ.l - aᵒ.l

        d.l -= aᵗ.l + aʰ.l - aᵒ.l
    end
    # update arrival and departure time
    rᵒ = r
    if isequal(φᵀ::Bool, true)
        if iscustomer(n) n.tᵃ, n.tᵈ = 0., 0. end
        for r ∈ v.R
            if r.tⁱ < rᵒ.tⁱ continue end
            if isopt(r)
                r.θⁱ = θⁱ
                r.θˢ = θⁱ + max(0., (r.l/v.lᵛ - r.θⁱ))
                r.tⁱ = tⁱ
                r.tˢ = r.tⁱ + v.τᶠ * (r.θˢ - r.θⁱ) + v.τᵈ * r.q
                cˢ = s.C[r.iˢ]
                cᵉ = s.C[r.iᵉ]
                tᵈ = r.tˢ
                cᵒ = cˢ
                while true
                    cᵒ.tᵃ = tᵈ + s.A[(cᵒ.iᵗ, cᵒ.iⁿ)].l/v.sᵛ
                    cᵒ.tᵈ = cᵒ.tᵃ + v.τᶜ + max(0., cᵒ.tᵉ - cᵒ.tᵃ - v.τᶜ) + cᵒ.τᶜ
                    if isequal(cᵒ, cᵉ) break end
                    tᵈ = cᵒ.tᵈ
                    cᵒ = s.C[cᵒ.iʰ]
                end
                r.θᵉ = r.θˢ - r.l/v.lᵛ
                r.tᵉ = cᵉ.tᵈ + s.A[(cᵉ.iⁿ, d.iⁿ)].l/v.sᵛ
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
        (v.tˢ, v.tᵉ) = isempty(v.R) ? (d.tˢ, d.tˢ) : (v.R[firstindex(v.R)].tⁱ, v.R[lastindex(v.R)].tᵉ)
        # update slack
        τ = d.tᵉ - v.tᵉ
        for r ∈ reverse(v.R)
            if !isopt(r) continue end
            cˢ = s.C[r.iˢ]
            cᵉ = s.C[r.iᵉ]
            cᵒ = cˢ
            while true
                τ  = min(τ, cᵒ.tˡ - cᵒ.tᵃ - v.τᶜ)
                if isequal(cᵒ, cᵉ) break end
                cᵒ = s.C[cᵒ.iʰ]
            end
            r.τ = τ
        end
        v.τ = min(τ, v.τ)
        d.τ = min(τ, d.τ)
    end
    return s
end



"""
    addroute(r::Route, s::Solution)

Returns `true` if a route `r` can be added into solution `s`.
"""
function addroute(r::Route, s::Solution)
    d = s.D[r.iᵈ]
    v = d.V[r.iᵛ]
    if any(!isopt, v.R) return false end
    if isequal(length(v.R), v.r̅) return false end
    return true
end
"""
    deleteroute(r::Route, s::Solution)

Returns `true` if route `r` can be deleted from solution `s`.
"""
function deleteroute(r::Route, s::Solution)
    if isopt(r) return false end
    return true
end



"""
    addvehicle(v::Vehicle, s::Solution)

Returns `true` if vehicle `v` can be added into solution `s`.
"""
function addvehicle(v::Vehicle, s::Solution)
    d = s.D[v.iᵈ]
    if any(!isopt, filter(v′ -> isidentical(v′, v), d.V)) return false end
    return true 
end
"""
    deletevehicle(v::Vehicle, s::Solution)

Returns `true` if vehicle `v` can be deleted from solution `s`.
"""
function deletevehicle(v::Vehicle, s::Solution)
    d = s.D[v.iᵈ]
    if isopt(v) return false end
    if isone(count(v′ -> isidentical(v′, v), d.V)) return false end
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
            r = Route(v, d)
            if addroute(r, s) push!(v.R, r) end
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
            r = Route(v, d)
            if addroute(r, s) push!(v.R, r) end
            v = Vehicle(v, d)
            if addvehicle(v, s) push!(d.V, v) end
        end
    end
    return s
end
"""
    postinsert!(s::Solution)

Post-insertion procedures. 
Returns solution `s` after deleting routes and vehicles, and subsequently correcting vehicle and route indices.
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
Returns solution `s` correcting vehicle and route indices.
"""
function postlocalsearch!(s::Solution)
    return s
end