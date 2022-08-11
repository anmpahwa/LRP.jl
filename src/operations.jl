# Insert operation
# insert customer node c between tail node nₜ and head node nₕ in route r
function insertnode!(c::CustomerNode, nₜ::Node, nₕ::Node, r::Route, s::Solution)
    A  = s.A
    V  = s.V
    aₒ = A[(nₜ.i, nₕ.i)]
    aₜ = A[(nₜ.i, c.i)]
    aₕ = A[(c.i, nₕ.i)]
    v  = V[r.o]
    Δˡ = aₜ.l + aₕ.l - aₒ.l
    Δᵗ = aₜ.t + aₕ.t - aₒ.t
    Δᶠ = aₜ.f + aₕ.f - aₒ.f
    Δᶜ = v.πᵐ * Δˡ + v.πʷ * Δᵗ + v.πᶠ * Δᶠ
    iscustomer(nₜ) ? nₜ.h = c.i : r.s = c.i
    iscustomer(nₕ) ? nₕ.t = c.i : r.e = c.i
    c.t = nₜ.i
    c.h = nₕ.i
    c.r = r
    r.n += 1
    r.q += c.q
    r.l += Δˡ
    r.t += Δᵗ
    r.f += Δᶠ
    r.c += Δᶜ
    return
end

# Remove operation
# remove customer node c between tail node nₜ and head node nₕ in route r
function removenode!(c::CustomerNode, nₜ::Node, nₕ::Node, r::Route, s::Solution)
    A  = s.A
    V  = s.V
    aₒ = A[(nₜ.i, nₕ.i)]
    aₜ = A[(nₜ.i, c.i)]
    aₕ = A[(c.i, nₕ.i)]
    v  = V[r.o]
    Δˡ = aₜ.l + aₕ.l - aₒ.l
    Δᵗ = aₜ.t + aₕ.t - aₒ.t
    Δᶠ = aₜ.f + aₕ.f - aₒ.f
    Δᶜ = v.πᵐ * Δˡ + v.πʷ * Δᵗ + v.πᶠ * Δᶠ
    iscustomer(nₜ) ? nₜ.h = c.h : r.s = c.h
    iscustomer(nₕ) ? nₕ.t = c.t : r.e = c.t
    c.t = 0
    c.h = 0
    # c.r = Route() # Note: It is not necessary to assign customer to a null route. Avoid for speedup.
    r.n -= 1
    r.q -= c.q
    r.l -= Δˡ
    r.t -= Δᵗ
    r.f -= Δᶠ
    r.c -= Δᶜ
    return
end