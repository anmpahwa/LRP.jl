# Insert node nᵒ between tail node nᵗ and head node nʰ in route rᵒ in solution s.
function insertnode!(nᵒ::Node, nᵗ::Node, nʰ::Node, rᵒ::Route, s::Solution)
    isdepot(nᵗ) ? rᵒ.iˢ = nᵒ.iⁿ : nᵗ.iʰ = nᵒ.iⁿ
    isdepot(nʰ) ? rᵒ.iᵉ = nᵒ.iⁿ : nʰ.iᵗ = nᵒ.iⁿ
    isdepot(nᵒ) ? (rᵒ.iˢ, rᵒ.iᵉ) = (nʰ.iⁿ, nᵗ.iⁿ) : (nᵒ.iʰ, nᵒ.iᵗ) = (nʰ.iⁿ, nᵗ.iⁿ)
    rᵒ.n += iscustomer(nᵒ) 
    rᵒ.q += iscustomer(nᵒ) ? nᵒ.q : 0
    rᵒ.l += s.A[(nᵗ.iⁿ, nᵒ.iⁿ)].l + s.A[(nᵒ.iⁿ, nʰ.iⁿ)].l - s.A[(nᵗ.iⁿ, nʰ.iⁿ)].l
    if iscustomer(nᵒ) nᵒ.r = rᵒ end
    return s
end

# Remove node nᵒ from its position between tail node nᵗ and head node nʰ in route rᵒ in solution s.
function removenode!(nᵒ::Node, nᵗ::Node, nʰ::Node, rᵒ::Route, s::Solution)
    isdepot(nᵗ) ? rᵒ.iˢ = nʰ.iⁿ : nᵗ.iʰ = nʰ.iⁿ
    isdepot(nʰ) ? rᵒ.iᵉ = nᵗ.iⁿ : nʰ.iᵗ = nᵗ.iⁿ
    isdepot(nᵒ) ? (rᵒ.iˢ, rᵒ.iᵉ) = (0, 0) : (nᵒ.iʰ, nᵒ.iᵗ) = (0, 0)
    rᵒ.n -= iscustomer(nᵒ) 
    rᵒ.q -= iscustomer(nᵒ) ? nᵒ.q : 0
    rᵒ.l -= s.A[(nᵗ.iⁿ, nᵒ.iⁿ)].l + s.A[(nᵒ.iⁿ, nʰ.iⁿ)].l - s.A[(nᵗ.iⁿ, nʰ.iⁿ)].l
    if iscustomer(nᵒ) nᵒ.r = NullRoute end
    return s
end

# Return true if vehicle v needs another route (adds conservatively)
function addroute(v::Vehicle, s::Solution)
    D = s.D
    d = D[v.iᵈ]
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
        if isequal(v.iᵈ, d.iⁿ) continue end
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