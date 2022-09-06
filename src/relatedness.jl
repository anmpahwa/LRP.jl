# Relatedness
function relatedness(d¹::DepotNode, d²::DepotNode, a::Arc)
    l = a.l
    q = abs(d¹.q - d².q)
    ϕ = false
    z = 1/(l - q - ϕ)
    return z
end
function relatedness(c¹::CustomerNode, c²::CustomerNode, a::Arc)
    l = a.l
    q = abs(c¹.q - c².q)
    ϕ = !isequal(c¹, c²) & isequal(c¹.r, c².r)
    z = 1/(l - q - ϕ)
    return z
end
function relatedness(c::CustomerNode, d::DepotNode, a::Arc)
    l = a.l
    q = 0
    ϕ = false
    r = c.r
    for v ∈ d.V 
        ϕ = isequal(r.iᵛ, v.iᵛ) 
        ϕ ? break : continue
    end
    z = 1/(l - q - ϕ)
    return z
end
function relatedness(d::DepotNode, c::CustomerNode, a::Arc)
    l = a.l
    q = 0
    ϕ = false
    r = c.r
    for v ∈ d.V 
        ϕ = isequal(r.iᵛ, v.iᵛ) 
        ϕ ? break : continue
    end
    z = 1/(l - q - ϕ)
    return z
end
function relatedness(r¹::Route, r²::Route)
    if !isopt(r¹) || !isopt(r²) return -Inf end
    if isequal(r¹, r²) return Inf end
    l = abs(r¹.l - r².l)
    q = abs(r¹.q - r².q)
    ϕ = isequal(r¹.iᵛ, r².iᵛ)
    z = 1/(l - q - ϕ)
    return z
end
function relatedness(v¹::Vehicle, v²::Vehicle)
    if !isopt(v¹) || !isopt(v²) return -Inf end
    if isequal(v¹, v²) return Inf end
    l¹ = 0.
    q¹ = 0.
    for r ∈ v¹.R 
        l¹ += r.l
        q¹ += r.q
    end
    l² = 0.
    q² = 0.
    for r ∈ v².R
        l² += r.l
        q² += r.q
    end
    l = abs(l¹ - l²)
    q = abs(q¹ - q²)
    ϕ = isequal(v¹.iᵈ, v².iᵈ)
    z = 1/(l - q - ϕ)
    return z
end