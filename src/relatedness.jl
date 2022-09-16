# Relatedness
function relatedness(d¹::DepotNode, d²::DepotNode, a::Arc)
    l = a.l
    t = 0
    q = abs(d¹.q - d².q)
    ϕ = false
    z = 1/(l + t - q - ϕ)
    return z
end
function relatedness(c¹::CustomerNode, c²::CustomerNode, a::Arc)
    l = a.l
    t = abs(c¹.tᵉ - c².tᵉ) + abs(c¹.tˡ - c².tˡ)
    q = abs(c¹.q - c².q)
    r¹ = c¹.r
    r² = c².r
    ϕʳ = isequal(r¹, r²)
    ϕᵛ = isequal(r¹.iᵛ, r².iᵛ)
    ϕᵈ = isequal(r¹.iᵈ, r².iᵈ)
    ϕ  = isequal(c¹, c²) ? 0 : (ϕʳ + ϕᵛ + ϕᵈ)
    z = 1/(l + t - q - ϕ)
    return z
end
function relatedness(c::CustomerNode, d::DepotNode, a::Arc)
    l = a.l
    t = 0
    q = 0
    ϕ = false
    r = c.r
    ϕ = isequal(r.iᵈ, d.iⁿ)
    z = 1/(l + t - q - ϕ)
    return z
end
function relatedness(d::DepotNode, c::CustomerNode, a::Arc)
    l = a.l
    t = 0
    q = 0
    ϕ = false
    r = c.r
    ϕ = isequal(r.iᵈ, d.iⁿ)
    z = 1/(l + t - q - ϕ)
    return z
end
function relatedness(r¹::Route, r²::Route)
    if !isopt(r¹) || !isopt(r²) return -Inf end
    if isequal(r¹, r²) return Inf end
    l = abs(r¹.l - r².l)
    t = abs(r¹.tˢ - r².tˢ) + abs(r¹.tᵉ - r².tᵉ)
    q = abs(r¹.q - r².q)
    ϕ = isequal(r¹.iᵛ, r².iᵛ)
    z = 1/(l + t - q - ϕ)
    return z
end
function relatedness(v¹::Vehicle, v²::Vehicle)
    if !isopt(v¹) || !isopt(v²) return -Inf end
    if isequal(v¹, v²) return Inf end
    l¹ = 0.
    t¹ = 0.
    q¹ = 0.
    for r ∈ v¹.R 
        l¹ += r.l
        t¹ += r.l/v¹.s
        q¹ += r.q
    end
    l² = 0.
    t² = 0.
    q² = 0.
    for r ∈ v².R
        l² += r.l
        t² += r.l/v².s
        q² += r.q
    end
    l = abs(l¹ - l²)
    t = abs(t¹ - t²)
    q = abs(q¹ - q²)
    ϕ = isequal(v¹.iᵈ, v².iᵈ)
    z = 1/(l + t - q - ϕ)
    return z
end