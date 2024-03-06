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
    s.πᵖ -= c.qᶜ
    if iscustomer(nᵗ) nᵗ.iʰ = c.iⁿ end
    if iscustomer(nʰ) nʰ.iᵗ = c.iⁿ end
    c.iᵗ  = nᵗ.iⁿ
    c.iʰ  = nʰ.iⁿ
    c.r   = r
    s.πᶠ += 0.
    s.πᵒ += 0.
    s.πᵖ += 0.
    # update associated route
    s.πᶠ -= 0.
    s.πᵒ -= r.l * v.πᵈ
    s.πᵖ -= (r.q > v.qᵛ) ? (r.q - v.qᵛ) : 0.
    s.πᵖ -= (r.l > v.lᵛ) ? (r.l - v.lᵛ) : 0.
    r.x   = (r.n * r.x + c.x)/(r.n + 1)
    r.y   = (r.n * r.y + c.y)/(r.n + 1)
    if isdepot(nᵗ) r.iˢ = c.iⁿ end
    if isdepot(nʰ) r.iᵉ = c.iⁿ end
    r.n  += 1
    r.q  += c.qᶜ
    r.l  += aᵗ.l + aʰ.l - aᵒ.l
    s.πᶠ += 0.
    s.πᵒ += r.l * v.πᵈ
    s.πᵖ += (r.q > v.qᵛ) ? (r.q - v.qᵛ) : 0.
    s.πᵖ += (r.l > v.lᵛ) ? (r.l - v.lᵛ) : 0.
    # update associated vehicle
    s.πᶠ -= isopt(v) * v.πᶠ
    s.πᵒ -= 0.
    s.πᵖ -= (length(v.R) > v.r̅) ? (length(v.R) - v.r̅) : 0.
    v.n  += 1
    v.q  += c.qᶜ
    v.l  += aᵗ.l + aʰ.l - aᵒ.l
    s.πᶠ += isopt(v) * v.πᶠ
    s.πᵒ += 0.
    s.πᵖ += (length(v.R) > v.r̅) ? (length(v.R) - v.r̅) : 0.
    # update associated depot
    s.πᶠ -= isopt(d) * d.πᶠ
    s.πᵒ -= d.q * d.πᵒ
    s.πᵖ -= (d.q > d.qᵈ) ? (d.q - d.qᵈ) : 0.
    d.n  += 1
    d.q  += c.qᶜ
    d.l  += aᵗ.l + aʰ.l - aᵒ.l
    s.πᶠ += isopt(d) * d.πᶠ
    s.πᵒ += d.q * d.πᵒ
    s.πᵖ += (d.q > d.qᵈ) ? (d.q - d.qᵈ) : 0.
    # update en-route parameters
    if isequal(s.φ, false) return s end
    s.πᶠ -= 0.
    s.πᵒ -= (v.tᵉ - v.tˢ) * v.πᵗ
    s.πᵖ -= (d.tˢ > v.tˢ) ? (d.tˢ - v.tˢ) : 0.
    s.πᵖ -= (v.tᵉ > d.tᵉ) ? (v.tᵉ - d.tᵉ) : 0.
    s.πᵖ -= ((v.tᵉ - v.tˢ) > v.τʷ) ? ((v.tᵉ - v.tˢ) - v.τʷ) : 0.
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
                s.πᶠ -= 0.
                s.πᵒ -= 0.
                s.πᵖ -= (c.tᵃ > c.tˡ) ? (c.tᵃ - c.tˡ) : 0.
                c.tᵃ  = tᵈ + s.A[(c.iᵗ, c.iⁿ)].l/v.sᵛ
                c.tᵈ  = c.tᵃ + v.τᶜ + max(0., c.tᵉ - c.tᵃ - v.τᶜ) + c.τᶜ
                s.πᶠ += 0.
                s.πᵒ += 0.
                s.πᵖ += (c.tᵃ > c.tˡ) ? (c.tᵃ - c.tˡ) : 0.
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
    s.πᶠ += 0.
    s.πᵒ += (v.tᵉ - v.tˢ) * v.πᵗ
    s.πᵖ += (d.tˢ > v.tˢ) ? (d.tˢ - v.tˢ) : 0.
    s.πᵖ += (v.tᵉ > d.tᵉ) ? (v.tᵉ - d.tᵉ) : 0.
    s.πᵖ += ((v.tᵉ - v.tˢ) > v.τʷ) ? ((v.tᵉ - v.tˢ) - v.τʷ) : 0.
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
    s.πᵖ -= (c.tᵃ > c.tˡ) ? (c.tᵃ - c.tˡ) : 0.
    if iscustomer(nᵗ) nᵗ.iʰ = nʰ.iⁿ end
    if iscustomer(nʰ) nʰ.iᵗ = nᵗ.iⁿ end
    c.iᵗ  = 0
    c.iʰ  = 0
    c.tᵃ  = 0.
    c.tᵈ  = 0.
    c.r   = NullRoute
    s.πᶠ += 0.
    s.πᵒ += 0.
    s.πᵖ += (c.tᵃ > c.tˡ) ? (c.tᵃ - c.tˡ) : 0.
    s.πᵖ += c.qᶜ
    # update associated route
    s.πᶠ -= 0.
    s.πᵒ -= r.l * v.πᵈ
    s.πᵖ -= (r.q > v.qᵛ) ? (r.q - v.qᵛ) : 0.
    s.πᵖ -= (r.l > v.lᵛ) ? (r.l - v.lᵛ) : 0.
    r.x   = isone(r.n) ? 0. : (r.n * r.x - c.x)/(r.n - 1)
    r.y   = isone(r.n) ? 0. : (r.n * r.y - c.y)/(r.n - 1)
    if isdepot(nᵗ) r.iˢ = nʰ.iⁿ end
    if isdepot(nʰ) r.iᵉ = nᵗ.iⁿ end
    r.n  -= 1
    r.q  -= c.qᶜ
    r.l  -= aᵗ.l + aʰ.l - aᵒ.l
    s.πᶠ += 0.
    s.πᵒ += r.l * v.πᵈ
    s.πᵖ += (r.q > v.qᵛ) ? (r.q - v.qᵛ) : 0.
    s.πᵖ += (r.l > v.lᵛ) ? (r.l - v.lᵛ) : 0.
    # update associated vehicle
    s.πᶠ -= isopt(v) * v.πᶠ
    s.πᵒ -= 0.
    s.πᵖ -= (length(v.R) > v.r̅) ? (length(v.R) - v.r̅) : 0.
    v.n  -= 1
    v.q  -= c.qᶜ
    v.l  -= aᵗ.l + aʰ.l - aᵒ.l
    s.πᶠ += isopt(v) * v.πᶠ
    s.πᵒ += 0.
    s.πᵖ += (length(v.R) > v.r̅) ? (length(v.R) - v.r̅) : 0.
    # update associated depot
    s.πᶠ -= isopt(d) * d.πᶠ
    s.πᵒ -= d.q * d.πᵒ
    s.πᵖ -= (d.q > d.qᵈ) ? (d.q - d.qᵈ) : 0.
    d.n  -= 1
    d.q  -= c.qᶜ
    d.l  -= aᵗ.l + aʰ.l - aᵒ.l
    s.πᶠ += isopt(d) * d.πᶠ
    s.πᵒ += d.q * d.πᵒ
    s.πᵖ += (d.q > d.qᵈ) ? (d.q - d.qᵈ) : 0.
    # update en-route variables
    if isequal(s.φ, false) return s end
    s.πᶠ -= 0.
    s.πᵒ -= (v.tᵉ - v.tˢ) * v.πᵗ
    s.πᵖ -= (d.tˢ > v.tˢ) ? (d.tˢ - v.tˢ) : 0.
    s.πᵖ -= (v.tᵉ > d.tᵉ) ? (v.tᵉ - d.tᵉ) : 0.
    s.πᵖ -= ((v.tᵉ - v.tˢ) > v.τʷ) ? ((v.tᵉ - v.tˢ) - v.τʷ) : 0.
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
                s.πᶠ -= 0.
                s.πᵒ -= 0.
                s.πᵖ -= (c.tᵃ > c.tˡ) ? (c.tᵃ - c.tˡ) : 0.
                c.tᵃ  = tᵈ + s.A[(c.iᵗ, c.iⁿ)].l/v.sᵛ
                c.tᵈ  = c.tᵃ + v.τᶜ + max(0., c.tᵉ - c.tᵃ - v.τᶜ) + c.τᶜ
                s.πᶠ += 0.
                s.πᵒ += 0.
                s.πᵖ += (c.tᵃ > c.tˡ) ? (c.tᵃ - c.tˡ) : 0.
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
    s.πᶠ += 0.
    s.πᵒ += (v.tᵉ - v.tˢ) * v.πᵗ
    s.πᵖ += (d.tˢ > v.tˢ) ? (d.tˢ - v.tˢ) : 0.
    s.πᵖ += (v.tᵉ > d.tᵉ) ? (v.tᵉ - d.tᵉ) : 0.
    s.πᵖ += ((v.tᵉ - v.tˢ) > v.τʷ) ? ((v.tᵉ - v.tˢ) - v.τʷ) : 0.
    return s
end



"""
    movevehicle!(v::Vehicle, d₁::DepotNode, d₂::DepotNode, s::Solution)

Returns solution `s` after moving vehicle `v` from fleet of `d₁` into fleet 
of depot node `d₂`.
"""
function movevehicle!(v::Vehicle, d₁::DepotNode, d₂::DepotNode, s::Solution)
    deleteat!(d₁.V, findfirst(isequal(v), d₁.V))
    push!(d₂.V, v)
    v.iᵈ = d₂.iⁿ
    for r ∈ v.R
        if isopt(r)
            nᵗ  = s.C[r.iᵉ]
            nʰ  = s.C[r.iˢ]
            aᵗ₁ = s.A[(nᵗ.iⁿ, d₁.iⁿ)]
            aᵗ₂ = s.A[(nᵗ.iⁿ, d₂.iⁿ)]
            aʰ₁ = s.A[(d₁.iⁿ, nʰ.iⁿ)]
            aʰ₂ = s.A[(d₂.iⁿ, nʰ.iⁿ)]
            # update depot d₁
            s.πᶠ -= isopt(d₁) * d₁.πᶠ
            s.πᵒ -= d₁.q * d₁.πᵒ
            s.πᵖ -= (d₁.q > d₁.qᵈ) ? (d₁.q - d₁.qᵈ) : 0.
            d₁.n -= r.n
            d₁.q -= r.q
            d₁.l -= r.l
            s.πᶠ += isopt(d₁) * d₁.πᶠ
            s.πᵒ += d₁.q * d₁.πᵒ
            s.πᵖ += (d₁.q > d₁.qᵈ) ? (d₁.q - d₁.qᵈ) : 0.
            # update associated customer nodes
            nᵗ.iʰ = d₂.iⁿ
            nʰ.iᵗ = d₂.iⁿ
            # update associated route
            s.πᶠ -= 0.
            s.πᵒ -= r.l * v.πᵈ
            s.πᵖ -= (r.q > v.qᵛ) ? (r.q - v.qᵛ) : 0.
            s.πᵖ -= (r.l > v.lᵛ) ? (r.l - v.lᵛ) : 0.
            r.iᵈ  = d₂.iⁿ
            r.iˢ  = nʰ.iⁿ
            r.iᵉ  = nᵗ.iⁿ
            r.l  -= aᵗ₁.l + aʰ₁.l
            r.l  += aᵗ₂.l + aʰ₂.l
            s.πᶠ += 0.
            s.πᵒ += r.l * v.πᵈ
            s.πᵖ += (r.q > v.qᵛ) ? (r.q - v.qᵛ) : 0.
            s.πᵖ += (r.l > v.lᵛ) ? (r.l - v.lᵛ) : 0.
            # update associated vehicle
            s.πᶠ -= isopt(v) * v.πᶠ
            s.πᵒ -= 0.
            s.πᵖ -= (length(v.R) > v.r̅) ? (length(v.R) - v.r̅) : 0.
            v.l  -= aᵗ₁.l + aʰ₁.l
            v.l  += aᵗ₂.l + aʰ₂.l
            s.πᶠ += isopt(v) * v.πᶠ
            s.πᵒ += 0.
            s.πᵖ += (length(v.R) > v.r̅) ? (length(v.R) - v.r̅) : 0.
            # update depot d₂
            s.πᶠ -= isopt(d₂) * d₂.πᶠ
            s.πᵒ -= d₂.q * d₂.πᵒ
            s.πᵖ -= (d₂.q > d₂.qᵈ) ? (d₂.q - d₂.qᵈ) : 0.
            d₂.n += r.n
            d₂.q += r.q
            d₂.l += r.l
            s.πᶠ += isopt(d₂) * d₂.πᶠ
            s.πᵒ += d₂.q * d₂.πᵒ
            s.πᵖ += (d₂.q > d₂.qᵈ) ? (d₂.q - d₂.qᵈ) : 0.
        else
            r.iᵈ = d₂.iⁿ
            r.iˢ = d₂.iⁿ
            r.iᵉ = d₂.iⁿ
        end
    end
    # update en-route variables
    if isequal(s.φ, false) return s end
    s.πᶠ -= 0.
    s.πᵒ -= (v.tᵉ - v.tˢ) * v.πᵗ
    s.πᵖ -= (d₁.tˢ > v.tˢ) ? (d₁.tˢ - v.tˢ) : 0.
    s.πᵖ -= (v.tᵉ > d₁.tᵉ) ? (v.tᵉ - d₁.tᵉ) : 0.
    s.πᵖ -= ((v.tᵉ - v.tˢ) > v.τʷ) ? ((v.tᵉ - v.tˢ) - v.τʷ) : 0.
    tⁱ = d₂.tˢ
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
                s.πᶠ -= 0.
                s.πᵒ -= 0.
                s.πᵖ -= (c.tᵃ > c.tˡ) ? (c.tᵃ - c.tˡ) : 0.
                c.tᵃ  = tᵈ + s.A[(c.iᵗ, c.iⁿ)].l/v.sᵛ
                c.tᵈ  = c.tᵃ + v.τᶜ + max(0., c.tᵉ - c.tᵃ - v.τᶜ) + c.τᶜ
                s.πᶠ += 0.
                s.πᵒ += 0.
                s.πᵖ += (c.tᵃ > c.tˡ) ? (c.tᵃ - c.tˡ) : 0.
                if isequal(c, cᵉ) break end
                tᵈ = c.tᵈ
                c  = s.C[c.iʰ]
            end
            r.θᵉ = r.θˢ - r.l/v.lᵛ
            r.tᵉ = cᵉ.tᵈ + s.A[(cᵉ.iⁿ, d₂.iⁿ)].l/v.sᵛ
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
    (v.tˢ, v.tᵉ) = isempty(v.R) ? (d₂.tˢ, d₂.tˢ) : (v.R[firstindex(v.R)].tⁱ, v.R[lastindex(v.R)].tᵉ)
    s.πᶠ += 0.
    s.πᵒ += (v.tᵉ - v.tˢ) * v.πᵗ
    s.πᵖ += (d₂.tˢ > v.tˢ) ? (d₂.tˢ - v.tˢ) : 0.
    s.πᵖ += (v.tᵉ > d₂.tᵉ) ? (v.tᵉ - d₂.tᵉ) : 0.
    s.πᵖ += ((v.tᵉ - v.tˢ) > v.τʷ) ? ((v.tᵉ - v.tˢ) - v.τʷ) : 0.
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
    # update slack
    if isequal(s.φ, false) return s end
    for d ∈ s.D
        τ = Inf
        for v ∈ d.V
            d.tᵉ - v.tᵉ
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
    if isequal(s.φ, false) return s end
    for d ∈ s.D
        τ = Inf
        for v ∈ d.V
            d.tᵉ - v.tᵉ
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
    # update slack
    if isequal(s.φ, false) return s end
    for d ∈ s.D
        τ = Inf
        for v ∈ d.V
            d.tᵉ - v.tᵉ
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
    if isequal(s.φ, false) return s end
    for d ∈ s.D
        τ = Inf
        for v ∈ d.V
            d.tᵉ - v.tᵉ
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