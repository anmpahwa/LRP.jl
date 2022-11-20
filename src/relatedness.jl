# Relatedness
function relatedness(c::CustomerNode, d::DepotNode, s::Solution)
    l = s.A[(c.iⁿ,d.iⁿ)].l
    t = 0
    q = 0
    r = c.r
    φ = 1 + isequal(r.iᵈ, d.iⁿ)
    z = (q + φ)/(l + t)
    return z
end
relatedness(d::DepotNode, c::CustomerNode, s::Solution) = relatedness(c, d, s)
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