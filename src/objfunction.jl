# Objective function evaluation
"""
    f(s::Solution; fixed=true, operational=true, penalty=true)

Objective function evaluation for solution `s`. Include `fixed`, 
`operational`, and `penalty` cost for constriant violation if `true`.
"""
function f(s::Solution; fixed=true, operational=true, penalty=true)
    πᶠ, πᵒ, πᵖ = 0., 0., 0.
    φᶠ, φᵒ, φᵖ = fixed, operational, penalty
    for d ∈ s.D
        πᶠ += isopt(d) * d.πᶠ
        qᵈ = 0
        nᵈ = 0
        for v ∈ d.V
            πᶠ += isopt(v) * v.πᶠ
            for r ∈ v.R 
                if !isopt(r) continue end
                qᵛ = r.q
                lᵛ = r.l
                qᵈ += qᵛ
                nᵈ += r.n
                πᵒ += r.l * v.πᵈ
                πᵖ += (qᵛ > v.q) * (qᵛ - v.q)                               # Vehicle capacity constraint
                πᵖ += (lᵛ > v.l) * (lᵛ - v.l)                               # Vehicle range constraint
            end
            πᵒ += (v.tᵉ - v.tˢ) * v.πᵗ
            πᵖ += (d.tˢ > v.tˢ) * (d.tˢ - v.tˢ)                             # Working-hours constraint (start time)
            πᵖ += (v.tᵉ > d.tᵉ) * (v.tᵉ - d.tᵉ)                             # Working-hours constraint (end time)
            πᵖ += (v.tᵉ - v.tˢ > v.τʷ) * (v.tᵉ - v.tˢ - v.τʷ)               # Working-hours constraint (duration)
        end
        pᵈ  = nᵈ/length(s.C)
        πᵒ += qᵈ * d.πᵒ
        πᵖ += (isone(s.φᴱ) && isone(d.jⁿ) && !isopt(d)) * d.πᶠ              # Depot use constraint
        πᵖ += (qᵈ > d.q) * (qᵈ - d.q)                                       # Depot capacity constraint
        πᵖ += (pᵈ < d.pˡ) * (d.pˡ - pᵈ)                                     # Depot customer share constraint
        πᵖ += (pᵈ > d.pᵘ) * (pᵈ - d.pᵘ)                                     # Depot customer share constraint
    end
    for c ∈ s.C πᵖ += isopen(c) ? 0. : (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) end    # Time-window constraint
    z = φᶠ * πᶠ + φᵒ * πᵒ + φᵖ * πᵖ * 10^(ceil(log10(πᶠ + πᵒ)))
    return z
end