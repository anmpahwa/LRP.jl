# Insert node nₒ between tail node nₜ and head node nₕ in route rₒ in solution s.
function insertnode!(nₒ::Node, nₜ::Node, nₕ::Node, rₒ::Route, s::Solution)
    isdepot(nₜ) ? rₒ.iₛ = nₒ.i : nₜ.iₕ = nₒ.i
    isdepot(nₕ) ? rₒ.iₑ = nₒ.i : nₕ.iₜ = nₒ.i
    isdepot(nₒ) ? (rₒ.iₛ, rₒ.iₑ) = (nₕ.i, nₜ.i) : (nₒ.iₕ, nₒ.iₜ) = (nₕ.i, nₜ.i)
    rₒ.n += iscustomer(nₒ) 
    rₒ.q += iscustomer(nₒ) ? nₒ.q : 0
    rₒ.l += s.A[(nₜ.i, nₒ.i)].l + s.A[(nₒ.i, nₕ.i)].l - s.A[(nₜ.i, nₕ.i)].l
    if iscustomer(nₒ) nₒ.r = rₒ end
    return s
end

# Remove node nₒ from its position between tail node nₜ and head node nₕ in route rₒ in solution s.
function removenode!(nₒ::Node, nₜ::Node, nₕ::Node, rₒ::Route, s::Solution)
    isdepot(nₜ) ? rₒ.iₛ = nₕ.i : nₜ.iₕ = nₕ.i
    isdepot(nₕ) ? rₒ.iₑ = nₜ.i : nₕ.iₜ = nₜ.i
    isdepot(nₒ) ? (rₒ.iₛ, rₒ.iₑ) = (0, 0) : (nₒ.iₕ, nₒ.iₜ) = (0, 0)
    rₒ.n -= iscustomer(nₒ) 
    rₒ.q -= iscustomer(nₒ) ? nₒ.q : 0
    rₒ.l -= s.A[(nₜ.i, nₒ.i)].l + s.A[(nₒ.i, nₕ.i)].l - s.A[(nₜ.i, nₕ.i)].l
    if iscustomer(nₒ) nₒ.r = NullRoute end
    return s
end

# Return true if vehicle v needs another route (adds conservatively)
function addroute(v::Vehicle, s::Solution)
    D = s.D
    d = D[v.o]
    # condtions when route mustn't be added
    if any(!isopt, v.R) return false end
    qᵈ = 0
    for v ∈ d.V for r ∈ v.R qᵈ += r.q end end
    if qᵈ ≥ d.q return false end
    # condition when route could be added
    if isempty(v.R) return true end
    for v ∈ d.V for r ∈ v.R if r.q > v.q return true end end end
    for d ∈ D 
        qᵈ = 0
        if isequal(v.o, d.i) continue end
        for v ∈ d.V for r ∈ v.R qᵈ += r.q end end
        if qᵈ > d.q return true end
    end
    return false
end

# Return true if route r can be deleted (deletes liberally)
function deleteroute(r::Route)
    if isopt(r) return false end
    return true
end