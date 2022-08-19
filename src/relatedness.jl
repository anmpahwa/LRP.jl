# Relatedness
function relatedness(d₁::DepotNode, d₂::DepotNode, a::Arc)
    l = a.l
    t = a.t
    f = a.f
    q = abs(d₁.q - d₂.q)
    ϕ = false
    z = 1/(0.482 * l + 35.0 * t + 3.826 * f - 1.0 * q - 1.0 * ϕ)
    return z
end
function relatedness(c₁::CustomerNode, c₂::CustomerNode, a::Arc)
    l = a.l
    t = a.t
    f = a.f
    q = abs(c₁.q - c₂.q)
    ϕ = !isequal(c₁, c₂) & isequal(c₁.r, c₂.r)
    z = 1/(0.482 * l + 35.0 * t + 3.826 * f - 1.0 * q - 1.0 * ϕ)
    return z
end
function relatedness(c::CustomerNode, d::DepotNode, a::Arc)
    l = a.l
    t = a.t
    f = a.f
    q = 0
    ϕ = false
    r = c.r
    for v ∈ d.V 
        ϕ = isequal(r.o, v.i) 
        ϕ ? break : continue
    end
    z = 1/(0.482 * l + 35.0 * t + 3.826 * f - 1.0 * q - 1.0 * ϕ)
    return z
end
function relatedness(d::DepotNode, c::CustomerNode, a::Arc)
    l = a.l
    t = a.t
    f = a.f
    q = 0
    ϕ = false
    r = c.r
    for v ∈ d.V 
        ϕ = isequal(r.o, v.i) 
        ϕ ? break : continue
    end
    z = 1/(0.482 * l + 35.0 * t + 3.826 * f - 1.0 * q - 1.0 * ϕ)
    return z
end
function relatedness(r₁::Route, r₂::Route)
    if isclose(r₁) || isclose(r₂) return -Inf end
    if isequal(r₁, r₂) return Inf end
    l = abs(r₁.l - r₂.l)
    t = abs(r₁.t - r₂.t)
    f = abs(r₁.f - r₂.f)
    q = abs(r₁.q - r₂.q)
    ϕ = isequal(r₁.o, r₂.o)
    z = 1/(0.482 * l + 35.0 * t + 3.826 * f - 1.0 * q - 1.0 * ϕ)
    return z
end
function relatedness(v₁::Vehicle, v₂::Vehicle)
    if isclose(v₁) || isclose(v₂) return -Inf end
    if isequal(v₁, v₂) return Inf end
    z = 0.
    for r₁ ∈ v₁.R 
        if isclose(r₁) continue end
        for r₂ ∈ v₂.R
            if isclose(r₂) continue end
            z += relatedness(r₁, r₂) 
        end 
    end
    return z
end