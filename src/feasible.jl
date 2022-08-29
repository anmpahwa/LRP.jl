"""
    isfeasible(s::Solution)

Is the solution feasible?
Returns true if node service constraint, node flow constraint,
sub-tour elimination, and capacity constraints are not violated.
"""
function isfeasible(s::Solution)
    C = s.C
    # Customer node service and flow constraints
    x = zeros(Int64, eachindex(C))
    for d ∈ s.D
        for v ∈ d.V
            for r ∈ v.R
                if !isopt(r) continue end
                cₛ = C[r.s]
                cₑ = C[r.e]
                c  = cₛ
                while true
                    k = c.i
                    x[k] += 1
                    if isequal(c, cₑ) break end
                    c = C[c.h]
                end
            end
        end
    end
    if any(!isone, x) return false end
    # Capacity constraints
    for d ∈ s.D
        qᵈ = 0
        for v ∈ d.V
            for r ∈ v.R 
                if !isopt(r) continue end
                qᵛ = r.q
                qᵈ += qᵛ
                if qᵛ > v.q return false end
            end
        end
        if qᵈ > d.q return false end
    end
    return true
end