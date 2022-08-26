# Objective function evaluation
"""
    f(s::Solution; fixed=true, operational=true, constraint=true)

Objective function evaluation for solution `s`.
Include `fixed`, `operational`, and `constraint` violation cost if `true`.
"""
function f(s::Solution; fixed=true, operational=true, constraint=true)
    zᶠ, zᵒ, zᶜ = 0., 0., 0.
    ϕᶠ, ϕᵒ, ϕᶜ = fixed, operational, constraint
    for d ∈ s.D
        if !isopt(d) continue end 
        qᵈ = 0
        zᶠ += d.πᶠ
        for v ∈ d.V 
            if !isopt(v) continue end
            zᶠ += v.πᶠ
            for r ∈ v.R 
                if !isopt(r) continue end
                qʳ  = r.q
                qᵈ += qʳ
                zᵒ += r.l * v.πᵒ
                zᶜ += (qʳ > v.q) * (qʳ - v.q)
            end
        end
        zᵒ += qᵈ * d.πᵒ
        zᶜ += (qᵈ > d.q) * (qᵈ - d.q)
    end
    z = ϕᶠ * zᶠ + ϕᵒ * zᵒ + ϕᶜ * zᶜ * (zᶠ + zᵒ)
    return z
end