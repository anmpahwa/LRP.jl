# Objective function evaluation
"""
    f(s::Solution; fixed=true, operational=true, penalty=true)

Objective function evaluation for solution `s`.
Include `fixed`, `operational`, and `penalty` cost for constriant 
violation if `true`.
"""
function f(s::Solution; fixed=true, operational=true, constraint=true)
    πᶠ, πᵒ, πᵖ = 0., 0., 0.
    ϕᶠ, ϕᵒ, ϕᵖ = fixed, operational, constraint
    for d ∈ s.D
        if !isopt(d) continue end 
        qᵈ = 0
        πᶠ += d.πᶠ
        for v ∈ d.V 
            if !isopt(v) continue end
            πᶠ += v.πᶠ
            for r ∈ v.R 
                if !isopt(r) continue end
                qᵛ  = r.q
                qᵈ += qᵛ
                πᵒ += r.l * v.πᵒ
                πᵖ += (qᵛ > v.q) * (qᵛ - v.q)
            end
        end
        πᵒ += qᵈ * d.πᵒ
        πᵖ += (qᵈ > d.q) * (qᵈ - d.q)
    end
    for c ∈ s.C πᵖ += isopen(c) ? 0. : (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) end
    z = ϕᶠ * πᶠ + ϕᵒ * πᵒ + ϕᵖ * πᵖ * (πᶠ + πᵒ)
    return z
end