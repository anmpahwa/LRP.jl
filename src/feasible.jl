"""
    isfeasible(s::Solution)

Returns true if 
node service, node flow, and sub-tour elimination constraints; 
depot and vehicle capacity constriants; 
vehicle range and working-hours constraints; and
time-window constraints
are not violated.
"""
function isfeasible(s::Solution)
    x = zeros(Int64, eachindex(s.C))
    n = length(s.C)
    for d ∈ s.D
        qᵈ = 0
        pᵈ = 0.
        for v ∈ d.V
            tˢ = 0.
            tᵉ = 0.
            for r ∈ v.R
                if !isopt(r) continue end
                qᵛ = r.q
                lᵛ = r.l
                if qᵛ > v.q return false end                                # Vehicle capacity constraint
                if lᵛ > v.l return false end                                # Vehicle range constraint
                qᵈ += r.q
                pᵈ += r.n/n
                tᵉ = r.tᵉ
                cˢ = s.C[r.iˢ]
                cᵉ = s.C[r.iᵉ]
                cᵒ = cˢ
                while true
                    if cᵒ.tᵃ > cᵒ.tˡ return false end                       # Time-window constraint
                    x[cᵒ.iⁿ] += 1
                    if isequal(cᵒ, cᵉ) break end
                    cᵒ = s.C[cᵒ.iʰ]
                end
            end
            tᵛ = tᵉ - tˢ
            if tᵛ > v.w return false end                                    # Working-hours constraint 
        end
        if (isone(s.ϕᴱ) && isone(d.jⁿ) && !isopt(d)) return false end       # Depot use constraint
        if qᵈ > d.q return false end                                        # Depot capacity constraint
        if !(d.pˡ ≤ pᵈ ≤ d.pᵘ) return false end                             # Depot customer share constraint
    end
    if any(!isone, x) return false end                                      # Node service, customer flow, and sub-tour elimination constrinat
    return true
end