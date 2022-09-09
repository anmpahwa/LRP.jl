# Objective function evaluation
"""
    f(s::Solution; fixed=true, operational=true, penalty=true)

Objective function evaluation for solution `s`. Include `fixed`, 
`operational`, and `penalty` cost for constriant violation if `true`.
"""
function f(s::Solution; fixed=true, operational=true, penalty=true)
    πᶠ, πᵒ, πᵖ = 0., 0., 0.
    ϕᶠ, ϕᵒ, ϕᵖ = fixed, operational, penalty
    for d ∈ s.D
        if !isopt(d) continue end 
        πᶠ += d.πᶠ
        qᵈ = 0
        for v ∈ d.V
            if !isopt(v) continue end 
            πᶠ += v.πᶠ
            tˢ = 0.
            tᵉ = 0.
            for r ∈ v.R 
                if !isopt(r) continue end
                qᵛ = r.q
                lᵛ = r.l
                tᵉ = r.tᵉ
                qᵈ += qᵛ
                πᵒ += r.l * v.πᵒ
                πᵖ += (qᵛ > v.q) * (qᵛ - v.q)
                πᵖ += (lᵛ > v.l) * (lᵛ - v.l)
            end
            tᵛ = tᵉ - tˢ
            πᵖ += (tᵛ > v.w) * (tᵛ - v.w)
        end
        πᵒ += qᵈ * d.πᵒ
        πᵖ += (qᵈ > d.q) * (qᵈ - d.q)
    end
    for c ∈ s.C πᵖ += isopen(c) ? 0. : (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) end
    z = ϕᶠ * πᶠ + ϕᵒ * πᵒ + ϕᵖ * πᵖ * (πᶠ + πᵒ)
    return z
end