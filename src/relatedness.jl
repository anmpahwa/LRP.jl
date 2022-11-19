# Relatedness
function relatedness(d¹::DepotNode, d²::DepotNode, a::Arc)
    l = a.l
    t = 0
    q = abs(d¹.q - d².q)
    φ = false
    z = (q + φ)/(l + t)
    return z
end
function relatedness(c¹::CustomerNode, c²::CustomerNode, a::Arc)
    l = a.l
    t = abs(c¹.tᵉ - c².tᵉ) + abs(c¹.tˡ - c².tˡ)
    q = abs(c¹.q - c².q)
    r¹ = c¹.r
    r² = c².r
    φʳ = isequal(r¹, r²)
    φᵛ = isequal(r¹.iᵛ, r².iᵛ)
    φᵈ = isequal(r¹.iᵈ, r².iᵈ)
    φ  = isequal(c¹, c²) ? 0 : (φʳ + φᵛ + φᵈ)
    z = (q + φ)/(l + t)
    return z
end
function relatedness(c::CustomerNode, d::DepotNode, a::Arc)
    l = a.l
    t = 0
    q = 0
    r = c.r
    φ = isequal(r.iᵈ, d.iⁿ)
    z = (q + φ)/(l + t)
    return z
end
function relatedness(d::DepotNode, c::CustomerNode, a::Arc)
    l = a.l
    t = 0
    q = 0
    r = c.r
    φ = isequal(r.iᵈ, d.iⁿ)
    z = (q + φ)/(l + t)
    return z
end
function relatedness(r¹::Route, r²::Route)
    if !isopt(r¹) || !isopt(r²) return -Inf end
    if isequal(r¹, r²) return Inf end
    l = abs(r¹.l - r².l)
    t = abs(r¹.tˢ - r².tˢ) + abs(r¹.tᵉ - r².tᵉ)
    q = abs(r¹.q - r².q)
    φ = isequal(r¹.iᵛ, r².iᵛ)
    z = (q + φ)/(l + t)
    return z
end
function relatedness(v¹::Vehicle, v²::Vehicle)
    if !isopt(v¹) || !isopt(v²) return -Inf end
    if isequal(v¹, v²) return Inf end
    l¹ = 0.
    t¹ = v¹.tᵉ - v¹.tˢ
    q¹ = 0.
    for r ∈ v¹.R 
        l¹ += r.l
        q¹ += r.q
    end
    l² = 0.
    t² = v².tᵉ - v².tˢ
    q² = 0.
    for r ∈ v².R
        l² += r.l
        q² += r.q
    end
    l = abs(l¹ - l²)
    t = abs(t¹ - t²)
    q = abs(q¹ - q²)
    φ = isequal(v¹.iᵈ, v².iᵈ)
    z = (q + φ)/(l + t)
    return z
end