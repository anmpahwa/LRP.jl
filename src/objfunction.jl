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
        πᶠ += isopt(d) * d.πᶠ
        qᵈ = 0
        nᵈ = 0
        for v ∈ d.V
            πᶠ += isopt(v) * v.πᶠ
            tˢ = 0.
            tᵉ = 0.
            for r ∈ v.R 
                if !isopt(r) continue end
                qᵛ = r.q
                lᵛ = r.l
                tᵉ = r.tᵉ
                qᵈ += qᵛ
                nᵈ += r.n
                πᵒ += r.l * v.πᵒ
                πᵖ += (qᵛ > v.q) * (qᵛ - v.q)                               # Vehicle capacity constraint
                πᵖ += (lᵛ > v.l) * (lᵛ - v.l)                               # Vehicle range constraint
            end
            tᵛ = tᵉ - tˢ
            πᵖ += (tᵛ > v.w) * (tᵛ - v.w)                                   # Working-hours constraint 
        end
        pᵈ  = nᵈ/length(s.C)
        πᵒ += qᵈ * d.πᵒ
        πᵖ += (isone(s.ϕᴱ) && isone(d.jⁿ) && !isopt(d)) * d.πᶠ              # Depot use constraint
        πᵖ += (qᵈ > d.q) * (qᵈ - d.q)                                       # Depot capacity constraint
        πᵖ += (pᵈ < d.pˡ) * (d.pˡ - pᵈ)                                     # Depot customer share constraint
        πᵖ += (pᵈ > d.pᵘ) * (pᵈ - d.pᵘ)                                     # Depot customer share constraint
    end
    for c ∈ s.C πᵖ += isopen(c) ? 0. : (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) end    # Time-window constraint
    z = ϕᶠ * πᶠ + ϕᵒ * πᵒ + ϕᵖ * πᵖ * 10^(ceil(log10(πᶠ + πᵒ)))
    return z
end