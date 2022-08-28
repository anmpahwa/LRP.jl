# Insert node nₒ between tail node nₜ and head node nₕ in solution s.
function insertnode!(nₒ::Node, nₜ::Node, nₕ::Node, r::Route, s::Solution)
    aₒ = s.A[(nₜ.i, nₕ.i)]
    aₜ = s.A[(nₜ.i, nₒ.i)]
    aₕ = s.A[(nₒ.i, nₕ.i)]
    isdepot(nₜ) ? r.s = nₒ.i : nₜ.h = nₒ.i
    isdepot(nₕ) ? r.e = nₒ.i : nₕ.t = nₒ.i
    isdepot(nₒ) ? (r.s, r.e) = (nₕ.i, nₜ.i) : (nₒ.h, nₒ.t) = (nₕ.i, nₜ.i)
    r.n += iscustomer(nₒ) 
    r.q += iscustomer(nₒ) * nₒ.q
    r.l += aₜ.l + aₕ.l - aₒ.l
    if iscustomer(nₒ) nₒ.r = r end
    return s
end

# Remove node nₒ from its position between tail node nₜ and head node nₕ in solution s.
function removenode!(nₒ::Node, nₜ::Node, nₕ::Node, r::Route, s::Solution)
    aₒ = s.A[(nₜ.i, nₕ.i)]
    aₜ = s.A[(nₜ.i, nₒ.i)]
    aₕ = s.A[(nₒ.i, nₕ.i)]
    isdepot(nₜ) ? r.s = nₕ.i : nₜ.h = nₕ.i
    isdepot(nₕ) ? r.e = nₜ.i : nₕ.t = nₜ.i
    isdepot(nₒ) ? (r.s, r.e) = (0, 0) : (nₒ.h, nₒ.t) = (0, 0)
    r.n -= iscustomer(nₒ) 
    r.q -= iscustomer(nₒ) * nₒ.q
    r.l -= aₜ.l + aₕ.l - aₒ.l
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