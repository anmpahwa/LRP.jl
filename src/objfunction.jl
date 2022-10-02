# Objective function evaluation
"""
    f(s::Solution; fixed=true, operational=true, penalty=true)

Objective function evaluation for solution `s`. Include `fixed`, 
`operational`, and `penalty` cost for constriant violation if `true`.
"""
function f(s::Solution; fixed=true, operational=true, penalty=true)
    πᶠ, πᵒ, πᵖ = 0., 0., 0.
    ϕᶠ, ϕᵒ, ϕᵖ = fixed, operational, penalty
    n = length(s.C)
    for d ∈ s.D
        πᶠ += isopt(d) * d.πᶠ
        qᵈ = 0
        pᵈ = 0.
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
                pᵈ += r.n/n
                πᵒ += r.l * v.πᵒ
                πᵖ += (qᵛ > v.q) * (qᵛ - v.q)
                πᵖ += (lᵛ > v.l) * (lᵛ - v.l)
            end
            tᵛ = tᵉ - tˢ
            πᵖ += (tᵛ > v.w) * (tᵛ - v.w)
        end
        πᵒ += qᵈ * d.πᵒ
        πᵖ += (isone(s.ϕᴱ) && isone(d.jⁿ) && !isopt(d)) * d.πᶠ
        πᵖ += (qᵈ > d.q) * (qᵈ - d.q)
        πᵖ += (pᵈ < d.pˡ) * (d.pˡ - pᵈ)
        πᵖ += (pᵈ > d.pᵘ) * (pᵈ - d.pᵘ)
    end
    for c ∈ s.C πᵖ += isopen(c) ? 0. : (c.tᵃ > c.tˡ) * (c.tᵃ - c.tˡ) end
    z = ϕᶠ * πᶠ + ϕᵒ * πᵒ + ϕᵖ * πᵖ * 10^(ceil(log10(πᶠ + πᵒ)))
    return z
end