"""
    insertnode!(c::CustomerNode, nᵗ::Node, nʰ::Node, r::Route, s::Solution)

Returns solution `s` after inserting customer node `c` between tail node `nᵗ` 
and head node `nʰ` in route `r`.
"""
function insertnode!(c::CustomerNode, nᵗ::Node, nʰ::Node, r::Route, s::Solution)
    d  = s.D[r.iᵈ]
    v  = d.V[r.iᵛ]
    aᵒ = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
    aᵗ = s.A[(nᵗ.iⁿ, c.iⁿ)]
    aʰ = s.A[(c.iⁿ, nʰ.iⁿ)]
    # update associated customer nodes
    s.πᶠ -= 0.
    s.πᵒ -= 0.
    s.πᵖ -= 1.
    c.iʳ  = r.iʳ
    c.iᵛ  = r.iᵛ
    c.iᵈ  = r.iᵈ
    isdepot(nᵗ) ? r.iˢ = c.iⁿ : nᵗ.iʰ = c.iⁿ
    isdepot(nʰ) ? r.iᵉ = c.iⁿ : nʰ.iᵗ = c.iⁿ
    c.iʰ  = nʰ.iⁿ
    c.iᵗ  = nᵗ.iⁿ
    c.r   = r
    s.πᶠ += 0.
    s.πᵒ += 0.
    # update associated route
    s.πᶠ -= 0.
    s.πᵒ -= r.l * v.πᵈ
    s.πᵖ -= (r.q > v.qᵛ) * (r.q / v.qᵛ)
    s.πᵖ -= (r.l > v.lᵛ) * (r.l / v.lᵛ)
    r.x   = (r.n * r.x + c.x)/(r.n + 1)
    r.y   = (r.n * r.y + c.y)/(r.n + 1)
    r.n  += 1
    r.q  += c.qᶜ
    r.l  += aᵗ.l + aʰ.l - aᵒ.l
    s.πᶠ += 0.
    s.πᵒ += r.l * v.πᵈ
    s.πᵖ += (r.q > v.qᵛ) * (r.q / v.qᵛ)
    s.πᵖ += (r.l > v.lᵛ) * (r.l / v.lᵛ)
    # update associated vehicle
    s.πᶠ -= isopt(v) ? v.πᶠ : 0.
    s.πᵖ -= (length(v.R) > v.r̅) * (length(v.R) / v.r̅)
    v.n  += 1
    v.q  += c.qᶜ
    v.l  += aᵗ.l + aʰ.l - aᵒ.l
    s.πᶠ += isopt(v) ? v.πᶠ : 0.
    s.πᵖ += (length(v.R) > v.r̅) * (length(v.R) / v.r̅)
    # update associated depot
    s.πᶠ -= isopt(d) ? d.πᶠ : 0.
    s.πᵒ -= d.q * d.πᵒ
    s.πᵖ -= (d.q > d.qᵈ) * (d.q / d.qᵈ)
    d.n  += 1
    d.q  += c.qᶜ
    d.l  += aᵗ.l + aʰ.l - aᵒ.l
    s.πᶠ += isopt(d) ? d.πᶠ : 0.
    s.πᵒ += d.q * d.πᵒ
    s.πᵖ += (d.q > d.qᵈ) * (d.q / d.qᵈ)
    # update en-route parameters
    if isequal(φᵉ::Bool, true)
        s.πᵒ -= (v.tᵉ - v.tˢ) * v.πᵗ
        s.πᵖ -= (d.tˢ > v.tˢ) * (d.tˢ / v.tˢ)
        s.πᵖ -= (v.tᵉ > d.tᵉ) * (v.tᵉ / d.tᵉ)
        s.πᵖ -= (v.tᵉ - v.tˢ > v.τʷ) * ((v.tᵉ - v.tˢ) / v.τʷ)
        tᵒ = r.tⁱ
        tⁱ = r.tⁱ
        θⁱ = r.θⁱ
        for r ∈ v.R
            if r.tⁱ < tᵒ continue end
            if isopt(r)
                r.θⁱ = θⁱ
                r.θˢ = θⁱ + max(0., (r.l/v.lᵛ - r.θⁱ))
                r.tⁱ = tⁱ
                r.tˢ = r.tⁱ + v.τᶠ * (r.θˢ - r.θⁱ) + v.τᵈ * r.q
                cˢ = s.C[r.iˢ]
                cᵉ = s.C[r.iᵉ]
                tᵈ = r.tˢ
                c  = cˢ
                while true
                    s.πᵖ -= (c.tᵃ > c.tˡ) * (c.tᵃ / c.tˡ)
                    c.tᵃ  = tᵈ + s.A[(c.iᵗ, c.iⁿ)].l/v.sᵛ
                    c.tᵈ  = c.tᵃ + v.τᶜ + max(0., c.tᵉ - c.tᵃ - v.τᶜ) + c.τᶜ
                    s.πᵖ += (c.tᵃ > c.tˡ) * (c.tᵃ / c.tˡ)
                    if isequal(c, cᵉ) break end
                    tᵈ = c.tᵈ
                    c  = s.C[c.iʰ]
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
        (v.tˢ, v.tᵉ) = isempty(v.R) ? (d.tˢ, d.tˢ) : (v.R[firstindex(v.R)].tⁱ, v.R[lastindex(v.R)].tᵉ)
        s.πᵒ += (v.tᵉ - v.tˢ) * v.πᵗ
        s.πᵖ += (d.tˢ > v.tˢ) * (d.tˢ / v.tˢ)
        s.πᵖ += (v.tᵉ > d.tᵉ) * (v.tᵉ / d.tᵉ)
        s.πᵖ += (v.tᵉ - v.tˢ > v.τʷ) * ((v.tᵉ - v.tˢ) / v.τʷ)
    end
    return s
end
"""
    removenode!(c::CustomerNode, nᵗ::Node, nʰ::Node, r::Route, s::Solution)

Returns solution `s` after removing customer node `c` between tail node `nᵗ` 
and head node `nʰ` in route `r`.
"""
function removenode!(c::CustomerNode, nᵗ::Node, nʰ::Node, r::Route, s::Solution)
    d  = s.D[r.iᵈ]
    v  = d.V[r.iᵛ]
    aᵒ = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
    aᵗ = s.A[(nᵗ.iⁿ, c.iⁿ)]
    aʰ = s.A[(c.iⁿ, nʰ.iⁿ)]
    # update associated customer nodes
    s.πᶠ -= 0.
    s.πᵒ -= 0.
    c.iʳ  = 0
    c.iᵛ  = 0
    c.iᵈ  = 0
    isdepot(nᵗ) ? r.iˢ = nʰ.iⁿ : nᵗ.iʰ = nʰ.iⁿ
    isdepot(nʰ) ? r.iᵉ = nᵗ.iⁿ : nʰ.iᵗ = nᵗ.iⁿ
    c.iʰ  = 0
    c.iᵗ  = 0
    c.r   = NullRoute
    s.πᶠ += 0.
    s.πᵒ += 0.
    s.πᵖ += 1.
    # update associated route
    s.πᶠ -= 0.
    s.πᵒ -= r.l * v.πᵈ
    s.πᵖ -= (r.q > v.qᵛ) * (r.q / v.qᵛ)
    s.πᵖ -= (r.l > v.lᵛ) * (r.l / v.lᵛ)
    r.x   = isone(r.n) ? 0. : (r.n * r.x - c.x)/(r.n - 1)
    r.y   = isone(r.n) ? 0. : (r.n * r.y - c.y)/(r.n - 1)
    r.n  -= 1
    r.q  -= c.qᶜ
    r.l  -= aᵗ.l + aʰ.l - aᵒ.l
    s.πᶠ += 0.
    s.πᵒ += r.l * v.πᵈ
    s.πᵖ += (r.q > v.qᵛ) * (r.q / v.qᵛ)
    s.πᵖ += (r.l > v.lᵛ) * (r.l / v.lᵛ)
    # update associated vehicle
    s.πᶠ -= isopt(v) ? v.πᶠ : 0.
    s.πᵖ -= (length(v.R) > v.r̅) * (length(v.R) / v.r̅)
    v.n  -= 1
    v.q  -= c.qᶜ
    v.l  -= aᵗ.l + aʰ.l - aᵒ.l
    s.πᶠ += isopt(v) ? v.πᶠ : 0.
    s.πᵖ += (length(v.R) > v.r̅) * (length(v.R) / v.r̅)
    # update associated depot
    s.πᶠ -= isopt(d) ? d.πᶠ : 0.
    s.πᵒ -= d.q * d.πᵒ
    s.πᵖ -= (d.q > d.qᵈ) * (d.q / d.qᵈ)
    d.n  -= 1
    d.q  -= c.qᶜ
    d.l  -= aᵗ.l + aʰ.l - aᵒ.l
    s.πᶠ += isopt(d) ? d.πᶠ : 0.
    s.πᵒ += d.q * d.πᵒ
    s.πᵖ += (d.q > d.qᵈ) * (d.q / d.qᵈ)
    # update en-route variables
    if isequal(φᵉ::Bool, true)
        s.πᵖ -= (c.tᵃ > c.tˡ) * (c.tᵃ / c.tˡ)
        c.tᵃ = 0.
        c.tᵈ = 0.
        s.πᵖ += (c.tᵃ > c.tˡ) * (c.tᵃ / c.tˡ)
        s.πᵒ -= (v.tᵉ - v.tˢ) * v.πᵗ
        s.πᵖ -= (d.tˢ > v.tˢ) * (d.tˢ / v.tˢ)
        s.πᵖ -= (v.tᵉ > d.tᵉ) * (v.tᵉ / d.tᵉ)
        s.πᵖ -= (v.tᵉ - v.tˢ > v.τʷ) * ((v.tᵉ - v.tˢ) / v.τʷ)
        tᵒ = r.tⁱ
        tⁱ = r.tⁱ
        θⁱ = r.θⁱ
        for r ∈ v.R
            if r.tⁱ < tᵒ continue end
            if isopt(r)
                r.θⁱ = θⁱ
                r.θˢ = θⁱ + max(0., (r.l/v.lᵛ - r.θⁱ))
                r.tⁱ = tⁱ
                r.tˢ = r.tⁱ + v.τᶠ * (r.θˢ - r.θⁱ) + v.τᵈ * r.q
                cˢ = s.C[r.iˢ]
                cᵉ = s.C[r.iᵉ]
                tᵈ = r.tˢ
                c  = cˢ
                while true
                    s.πᵖ -= (c.tᵃ > c.tˡ) * (c.tᵃ / c.tˡ)
                    c.tᵃ  = tᵈ + s.A[(c.iᵗ, c.iⁿ)].l/v.sᵛ
                    c.tᵈ  = c.tᵃ + v.τᶜ + max(0., c.tᵉ - c.tᵃ - v.τᶜ) + c.τᶜ
                    s.πᵖ += (c.tᵃ > c.tˡ) * (c.tᵃ / c.tˡ)
                    if isequal(c, cᵉ) break end
                    tᵈ = c.tᵈ
                    c  = s.C[c.iʰ]
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
        (v.tˢ, v.tᵉ) = isempty(v.R) ? (d.tˢ, d.tˢ) : (v.R[firstindex(v.R)].tⁱ, v.R[lastindex(v.R)].tᵉ)
        s.πᵒ += (v.tᵉ - v.tˢ) * v.πᵗ
        s.πᵖ += (d.tˢ > v.tˢ) * (d.tˢ / v.tˢ)
        s.πᵖ += (v.tᵉ > d.tᵉ) * (v.tᵉ / d.tᵉ)
        s.πᵖ += (v.tᵉ - v.tˢ > v.τʷ) * ((v.tᵉ - v.tˢ) / v.τʷ)
    end
    return s
end



"""
    movevehicle!(v::Vehicle, d₁::DepotNode, d₂::DepotNode, s::Solution)

Returns solution `s` after moving vehicle `v` from fleet of `d₁` into fleet 
of depot node `d₂`.
"""
function movevehicle!(v::Vehicle, d₁::DepotNode, d₂::DepotNode, s::Solution)
    # remove vehicle from d₁
    d = d₁
    deleteat!(d.V, findfirst(isequal(v), d.V))
    v.iᵈ = 0
    for r ∈ v.R
        if isopt(r)
            nᵗ = s.C[r.iᵉ]
            nʰ = s.C[r.iˢ]
            aᵒ = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
            aᵗ = s.A[(nᵗ.iⁿ, d.iⁿ)]
            aʰ = s.A[(d.iⁿ, nʰ.iⁿ)]
            # update associated customer nodes
            nᵗ.iʰ = nʰ.iⁿ
            nʰ.iᵗ = nᵗ.iⁿ
            cˢ = nʰ
            cᵉ = nᵗ
            c  = cˢ
            while true
                c.iᵈ = 0
                if isequal(c, cᵉ) break end
                c = s.C[c.iʰ]
            end
            # update associated route
            s.πᶠ -= 0.
            s.πᵒ -= r.l * v.πᵈ
            s.πᵖ -= (r.q > v.qᵛ) * (r.q / v.qᵛ)
            s.πᵖ -= (r.l > v.lᵛ) * (r.l / v.lᵛ)
            r.iᵈ  = 0
            r.iˢ  = nʰ.iⁿ
            r.iᵉ  = nᵗ.iⁿ
            r.l  -= aᵗ.l + aʰ.l - aᵒ.l
            s.πᶠ += 0.
            s.πᵒ += r.l * v.πᵈ
            s.πᵖ += (r.q > v.qᵛ) * (r.q / v.qᵛ)
            s.πᵖ += (r.l > v.lᵛ) * (r.l / v.lᵛ)
            # update associated vehicle
            s.πᶠ -= isopt(v) ? v.πᶠ : 0.
            s.πᵖ -= (length(v.R) > v.r̅) * (length(v.R) / v.r̅)
            v.l  -= aᵗ.l + aʰ.l - aᵒ.l
            s.πᶠ += isopt(v) ? v.πᶠ : 0.
            s.πᵖ += (length(v.R) > v.r̅) * (length(v.R) / v.r̅)
            # update associated depot
            s.πᶠ -= isopt(d) ? d.πᶠ : 0.
            s.πᵒ -= d.q * d.πᵒ
            s.πᵖ -= (d.q > d.qᵈ) * (d.q / d.qᵈ)
            d.n  -= r.n
            d.q  -= r.q
            d.l  -= aᵗ.l + aʰ.l - aᵒ.l
            s.πᶠ += isopt(d) ? d.πᶠ : 0.
            s.πᵒ += d.q * d.πᵒ
            s.πᵖ += (d.q > d.qᵈ) * (d.q / d.qᵈ)
        else
            # update associated route
            s.πᶠ -= 0.
            s.πᵒ -= r.l * v.πᵈ
            s.πᵖ -= (r.q > v.qᵛ) * (r.q / v.qᵛ)
            s.πᵖ -= (r.l > v.lᵛ) * (r.l / v.lᵛ)
            r.iᵈ  = 0
            r.iˢ  = 0
            r.iᵉ  = 0
            s.πᶠ += 0.
            s.πᵒ += r.l * v.πᵈ
            s.πᵖ += (r.q > v.qᵛ) * (r.q / v.qᵛ)
            s.πᵖ += (r.l > v.lᵛ) * (r.l / v.lᵛ)
        end
    end
    # add vehicle to d₂
    d = d₂
    push!(d.V, v)
    v.iᵈ = d.iⁿ
    for r ∈ v.R
        if isopt(r)
            nᵗ = s.C[r.iᵉ]
            nʰ = s.C[r.iˢ]
            aᵒ = s.A[(nᵗ.iⁿ, nʰ.iⁿ)]
            aᵗ = s.A[(nᵗ.iⁿ, d.iⁿ)]
            aʰ = s.A[(d.iⁿ, nʰ.iⁿ)]
            # update associated customer nodes
            nᵗ.iʰ = d.iⁿ
            nʰ.iᵗ = d.iⁿ
            cˢ = nʰ
            cᵉ = nᵗ
            c  = cˢ
            while true
                c.iᵈ = d.iⁿ
                if isequal(c, cᵉ) break end
                c = s.C[c.iʰ]
            end
            # update associated route
            s.πᶠ -= 0.
            s.πᵒ -= r.l * v.πᵈ
            s.πᵖ -= (r.q > v.qᵛ) * (r.q / v.qᵛ)
            s.πᵖ -= (r.l > v.lᵛ) * (r.l / v.lᵛ)
            r.iᵈ  = d.iⁿ
            r.iˢ  = nʰ.iⁿ
            r.iᵉ  = nᵗ.iⁿ
            r.l  += aᵗ.l + aʰ.l - aᵒ.l
            s.πᶠ += 0.
            s.πᵒ += r.l * v.πᵈ
            s.πᵖ += (r.q > v.qᵛ) * (r.q / v.qᵛ)
            s.πᵖ += (r.l > v.lᵛ) * (r.l / v.lᵛ)
            # update associated vehicle
            s.πᶠ -= isopt(v) ? v.πᶠ : 0.
            s.πᵖ -= (length(v.R) > v.r̅) * (length(v.R) / v.r̅)
            v.l  += aᵗ.l + aʰ.l - aᵒ.l
            s.πᶠ += isopt(v) ? v.πᶠ : 0.
            s.πᵖ += (length(v.R) > v.r̅) * (length(v.R) / v.r̅)
            # update associated depot
            s.πᶠ -= isopt(d) ? d.πᶠ : 0.
            s.πᵒ -= d.q * d.πᵒ
            s.πᵖ -= (d.q > d.qᵈ) * (d.q / d.qᵈ)
            d.n  += r.n
            d.q  += r.q
            d.l  += aᵗ.l + aʰ.l - aᵒ.l
            s.πᶠ += isopt(d) ? d.πᶠ : 0.
            s.πᵒ += d.q * d.πᵒ
            s.πᵖ += (d.q > d.qᵈ) * (d.q / d.qᵈ)
        else
            # update associated route
            s.πᶠ -= 0.
            s.πᵒ -= r.l * v.πᵈ
            s.πᵖ -= (r.q > v.qᵛ) * (r.q / v.qᵛ)
            s.πᵖ -= (r.l > v.lᵛ) * (r.l / v.lᵛ)
            r.iᵈ  = d.iⁿ
            r.iˢ  = d.iⁿ
            r.iᵉ  = d.iⁿ
            s.πᶠ += 0.
            s.πᵒ += r.l * v.πᵈ
            s.πᵖ += (r.q > v.qᵛ) * (r.q / v.qᵛ)
            s.πᵖ += (r.l > v.lᵛ) * (r.l / v.lᵛ)
        end
    end
    # update en-route variables
    if isequal(φᵉ::Bool, true)
        s.πᵒ -= (v.tᵉ - v.tˢ) * v.πᵗ
        s.πᵖ -= (d.tˢ > v.tˢ) * (d.tˢ / v.tˢ)
        s.πᵖ -= (v.tᵉ > d.tᵉ) * (v.tᵉ / d.tᵉ)
        s.πᵖ -= (v.tᵉ - v.tˢ > v.τʷ) * ((v.tᵉ - v.tˢ) / v.τʷ)
        tⁱ = d.tˢ
        θⁱ = 1.
        for r ∈ v.R
            if isopt(r)
                r.θⁱ = θⁱ
                r.θˢ = θⁱ + max(0., (r.l/v.lᵛ - r.θⁱ))
                r.tⁱ = tⁱ
                r.tˢ = r.tⁱ + v.τᶠ * (r.θˢ - r.θⁱ) + v.τᵈ * r.q
                cˢ = s.C[r.iˢ]
                cᵉ = s.C[r.iᵉ]
                tᵈ = r.tˢ
                c  = cˢ
                while true
                    s.πᵖ -= (c.tᵃ > c.tˡ) * (c.tᵃ / c.tˡ)
                    c.tᵃ  = tᵈ + s.A[(c.iᵗ, c.iⁿ)].l/v.sᵛ
                    c.tᵈ  = c.tᵃ + v.τᶜ + max(0., c.tᵉ - c.tᵃ - v.τᶜ) + c.τᶜ
                    s.πᵖ += (c.tᵃ > c.tˡ) * (c.tᵃ / c.tˡ)
                    if isequal(c, cᵉ) break end
                    tᵈ = c.tᵈ
                    c  = s.C[c.iʰ]
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
        (v.tˢ, v.tᵉ) = isempty(v.R) ? (d.tˢ, d.tˢ) : (v.R[firstindex(v.R)].tⁱ, v.R[lastindex(v.R)].tᵉ)
        s.πᵒ += (v.tᵉ - v.tˢ) * v.πᵗ
        s.πᵖ += (d.tˢ > v.tˢ) * (d.tˢ / v.tˢ)
        s.πᵖ += (v.tᵉ > d.tᵉ) * (v.tᵉ / d.tᵉ)
        s.πᵖ += (v.tᵉ - v.tˢ > v.τʷ) * ((v.tᵉ - v.tˢ) / v.τʷ)
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
    if any(!isopt, filter(v′ -> isequal(v′.jᵛ, v.jᵛ), d.V)) return false end
    return true 
end
"""
    deletevehicle(v::Vehicle, s::Solution)

Returns `true` if vehicle `v` can be deleted from solution `s`.
"""
function deletevehicle(v::Vehicle, s::Solution)
    d = s.D[v.iᵈ]
    if isopt(v) return false end
    if isone(count(v′ -> isequal(v′.jᵛ, v.jᵛ), d.V)) return false end
    return true
end



"""
    preinitialize!(s::Solution)

Returns solution `s` after performing pre-intialization procedures.
Adds new routes and vehicles into the solution.
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

Returns solution `s` after performing post-intialization procedures. 
Deletes routes and vehicles if possible, and subsequently updates indices.
Additionally, updates route, vehicle, and depot slack time.
"""
function postinitialize!(s::Solution)
    # update indices
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
    for c ∈ s.C c.iᵛ, c.iʳ = c.r.iᵛ, c.r.iʳ end
    # update slack
    if isequal(φᵉ::Bool, false) return s end
    for d ∈ s.D
        τ = Inf
        for v ∈ d.V
            τ = d.tᵉ - v.tᵉ
            for r ∈ reverse(v.R)
                if !isopt(r) continue end
                cˢ = s.C[r.iˢ]
                cᵉ = s.C[r.iᵉ]
                c  = cˢ
                while true
                    τ = min(τ, c.tˡ - c.tᵃ - v.τᶜ)
                    if isequal(c, cᵉ) break end
                    c = s.C[c.iʰ]
                end
                r.τ = τ
            end
            v.τ = τ
        end
        d.τ = τ
    end
    return s
end



"""
    preremove!(s::Solution)

Returns solution `s` after performing pre-removal procedures. 
"""
function preremove!(s::Solution)
    return s
end
"""
    postremove!(s::Solution)

Returns solution `s` after performing post-removal procedures.
Updates route, vehicle, and depot slack time.
"""
function postremove!(s::Solution)
    if isequal(φᵉ::Bool, false) return s end
    for d ∈ s.D
        τ = Inf
        for v ∈ d.V
            τ = d.tᵉ - v.tᵉ
            for r ∈ reverse(v.R)
                if !isopt(r) continue end
                cˢ = s.C[r.iˢ]
                cᵉ = s.C[r.iᵉ]
                c  = cˢ
                while true
                    τ = min(τ, c.tˡ - c.tᵃ - v.τᶜ)
                    if isequal(c, cᵉ) break end
                    c = s.C[c.iʰ]
                end
                r.τ = τ
            end
            v.τ = τ
        end
        d.τ = τ
    end
    return s
end



"""
    preinsert!(s::Solution)

Returns solution `s` after performing pre-insertion procedures.
Adds new routes and vehicles into the solution.
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

Returns solution `s` after performing post-insertion procedures. 
Deletes routes and vehicles if possible, and subsequently updates indices.
Additionally, updates route, vehicle, and depot slack time.
"""
function postinsert!(s::Solution)
    # update indices
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
    for c ∈ s.C c.iᵛ, c.iʳ = c.r.iᵛ, c.r.iʳ end
    # update slack
    if isequal(φᵉ::Bool, false) return s end
    for d ∈ s.D
        τ = Inf
        for v ∈ d.V
            τ = d.tᵉ - v.tᵉ
            for r ∈ reverse(v.R)
                if !isopt(r) continue end
                cˢ = s.C[r.iˢ]
                cᵉ = s.C[r.iᵉ]
                c  = cˢ
                while true
                    τ = min(τ, c.tˡ - c.tᵃ - v.τᶜ)
                    if isequal(c, cᵉ) break end
                    c = s.C[c.iʰ]
                end
                r.τ = τ
            end
            v.τ = τ
        end
        d.τ = τ
    end
    return s
end



"""
    prelocalsearch!(s::Solution)

Returns solution `s` after performing pre-localsearch procedures. 
"""
function prelocalsearch!(s::Solution)
    return s
end
"""
    postlocalsearch!(s::Solution)

Returns solution `s` after performing post-localsearch procedures.
Updates route, vehicle, and depot slack time.
"""
function postlocalsearch!(s::Solution)
    if isequal(φᵉ::Bool, false) return s end
    for d ∈ s.D
        τ = Inf
        for v ∈ d.V
            τ = d.tᵉ - v.tᵉ
            for r ∈ reverse(v.R)
                if !isopt(r) continue end
                cˢ = s.C[r.iˢ]
                cᵉ = s.C[r.iᵉ]
                c  = cˢ
                while true
                    τ = min(τ, c.tˡ - c.tᵃ - v.τᶜ)
                    if isequal(c, cᵉ) break end
                    c = s.C[c.iʰ]
                end
                r.τ = τ
            end
            v.τ = τ
        end
        d.τ = τ
    end
    return s
end