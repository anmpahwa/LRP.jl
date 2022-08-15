# INSERT OPERATION
# Insert customer node nₒ between depot tail node nₜ and depot head node nₕ in route r
function insertnode!(nₒ::CustomerNode, nₜ::DepotNode, nₕ::DepotNode, r::Route, s::Solution)
    if isclose(nₒ) throw(ArgumentError("Customer node nₒ ($(nₒ.i)) is closed")) end
    if s.V[r.o] ∉ nₜ.V throw(ArgumentError("Route r $(r.i) does not belong to a vehicle origination from the depot tail node nₜ ($(nₜ.i))")) end
    if s.V[r.o] ∉ nₕ.V throw(ArgumentError("Route r $(r.i) does not belong to a vehicle origination from the depot head node nₕ ($(nₕ.i))")) end
    if !isequal(r.e, nₜ.i) throw(ArgumentError("Tail index of depot head node nₕ ($(nₕ.i)) does not match the index of depot tail node nₜ ($(nₜ.i)) ")) end
    if !isequal(r.s, nₕ.i) throw(ArgumentError("Head index of depot tail node nₜ ($(nₜ.i)) does not match the index of depot head node nₕ ($(nₕ.i)) ")) end
    A = s.A
    V = s.V
    v = V[r.o]
    aₒ = A[(nₜ.i, nₕ.i)]
    aₜ = A[(nₜ.i, nₒ.i)]
    aₕ = A[(nₒ.i, nₕ.i)]
    Δˡ = aₜ.l + aₕ.l - aₒ.l
    Δᵗ = aₜ.t + aₕ.t - aₒ.t
    Δᶠ = aₜ.f + aₕ.f - aₒ.f
    Δᶜ = v.πᵐ * Δˡ + v.πʷ * Δᵗ + v.πᶠ * Δᶠ
    r.s  = nₒ.i
    r.e  = nₒ.i
    nₒ.t = nₜ.i
    nₒ.h = nₕ.i
    nₒ.r = r
    r.n += 1
    r.q += nₒ.q
    r.l += Δˡ
    r.t += Δᵗ
    r.f += Δᶠ
    r.c += Δᶜ
    return s
end
# Insert customer node nₒ between customer tail node nₜ and customer head node nₕ in route r
function insertnode!(nₒ::CustomerNode, nₜ::CustomerNode, nₕ::CustomerNode, r::Route, s::Solution)
    if isclose(nₒ) throw(ArgumentError("Customer node nₒ ($(nₒ.i)) is closed")) end
    if !isequal(nₜ.r, r) throw(ArgumentError("Customer tail node nₜ ($(nₜ.i)) does not belong to route r ($(r.i)")) end
    if !isequal(nₕ.r, r) throw(ArgumentError("Customer head node nₕ ($(nₕ.i)) does not belong to route r ($(r.i)")) end
    if !isequal(nₜ.h, nₕ.i) throw(ArgumentError("Head index of customer tail node nₜ ($(nₜ.i)) does not match the index of customer head node nₕ ($(nₕ.i)) ")) end
    if !isequal(nₕ.t, nₜ.i) throw(ArgumentError("Tail index of customer head node nₕ ($(nₕ.i)) does not match the index of customer tail node nₜ ($(nₜ.i)) ")) end
    A = s.A
    V = s.V
    v = V[r.o]
    aₒ = A[(nₜ.i, nₕ.i)]
    aₜ = A[(nₜ.i, nₒ.i)]
    aₕ = A[(nₒ.i, nₕ.i)]
    Δˡ = aₜ.l + aₕ.l - aₒ.l
    Δᵗ = aₜ.t + aₕ.t - aₒ.t
    Δᶠ = aₜ.f + aₕ.f - aₒ.f
    Δᶜ = v.πᵐ * Δˡ + v.πʷ * Δᵗ + v.πᶠ * Δᶠ
    nₜ.h = nₒ.i
    nₕ.t = nₒ.i
    nₒ.t = nₜ.i
    nₒ.h = nₕ.i
    nₒ.r = r
    r.n += 1
    r.q += nₒ.q
    r.l += Δˡ
    r.t += Δᵗ
    r.f += Δᶠ
    r.c += Δᶜ
    return s
end
# Insert customer node nₒ between depot node nₜ and customer node nₕ in route r
function insertnode!(nₒ::CustomerNode, nₜ::DepotNode, nₕ::CustomerNode, r::Route, s::Solution)
    if isclose(nₒ) throw(ArgumentError("Customer node nₒ ($(nₒ.i)) is closed")) end
    if s.V[r.o] ∉ nₜ.V throw(ArgumentError("Route r $(r.i) does not belong to a vehicle origination from the depot tail node nₜ ($(nₜ.i))")) end
    if !isequal(nₕ.r, r) throw(ArgumentError("Customer head node nₕ ($(nₕ.i)) does not belong to route r ($(r.i)")) end
    if !isequal(nₕ.t, nₜ.i) throw(ArgumentError("Tail index of customer head node nₕ ($(nₕ.i)) does not match the index of depot tail node nₜ ($(nₜ.i)) ")) end
    if !isequal(r.s, nₕ.i) throw(ArgumentError("Head index of depot tail node nₜ ($(nₜ.i)) does not match the index of customer head node nₕ ($(nₕ.i)) ")) end
    A = s.A
    V = s.V
    v = V[r.o]
    aₒ = A[(nₜ.i, nₕ.i)]
    aₜ = A[(nₜ.i, nₒ.i)]
    aₕ = A[(nₒ.i, nₕ.i)]
    Δˡ = aₜ.l + aₕ.l - aₒ.l
    Δᵗ = aₜ.t + aₕ.t - aₒ.t
    Δᶠ = aₜ.f + aₕ.f - aₒ.f
    Δᶜ = v.πᵐ * Δˡ + v.πʷ * Δᵗ + v.πᶠ * Δᶠ
    r.s  = nₒ.i
    nₕ.t = nₒ.i
    nₒ.t = nₜ.i
    nₒ.h = nₕ.i
    nₒ.r = r
    r.n += 1
    r.q += nₒ.q
    r.l += Δˡ
    r.t += Δᵗ
    r.f += Δᶠ
    r.c += Δᶜ
    return s
end
# Insert customer node nₒ between tail customer node nₜ and head customer node nₕ in route r
function insertnode!(nₒ::CustomerNode, nₜ::CustomerNode, nₕ::DepotNode, r::Route, s::Solution)
    if isclose(nₒ) throw(ArgumentError("Customer node nₒ ($(nₒ.i)) is closed")) end
    if !isequal(nₜ.r, r) throw(ArgumentError("Customer tail node nₜ ($(nₜ.i)) does not belong to route r ($(r.i)")) end
    if s.V[r.o] ∉ nₕ.V throw(ArgumentError("Route r $(r.i) does not belong to a vehicle origination from the depot head node nₕ ($(nₕ.i))")) end
    if !isequal(r.e, nₜ.i) throw(ArgumentError("Tail index of depot head node nₕ ($(nₕ.i)) does not match the index of customer tail node nₜ ($(nₜ.i)) ")) end
    if !isequal(nₜ.h, nₕ.i) throw(ArgumentError("Head index of customer tail node nₜ ($(nₜ.i)) does not match the index of depot head node nₕ ($(nₕ.i)) ")) end
    A = s.A
    V = s.V
    v = V[r.o]
    aₒ = A[(nₜ.i, nₕ.i)]
    aₜ = A[(nₜ.i, nₒ.i)]
    aₕ = A[(nₒ.i, nₕ.i)]
    Δˡ = aₜ.l + aₕ.l - aₒ.l
    Δᵗ = aₜ.t + aₕ.t - aₒ.t
    Δᶠ = aₜ.f + aₕ.f - aₒ.f
    Δᶜ = v.πᵐ * Δˡ + v.πʷ * Δᵗ + v.πᶠ * Δᶠ
    nₜ.h = nₒ.i
    r.e  = nₒ.i
    nₒ.t = nₜ.i
    nₒ.h = nₕ.i
    nₒ.r = r
    r.n += 1
    r.q += nₒ.q
    r.l += Δˡ
    r.t += Δᵗ
    r.f += Δᶠ
    r.c += Δᶜ
    return s
end
# Insert depot node nₒ between customer tail node nₜ and customer head node nₕ in route r
function insertnode!(nₒ::DepotNode, nₜ::CustomerNode, nₕ::CustomerNode, r::Route, s::Solution)
    if s.V[r.o] ∉ nₒ.V throw(ArgumentError("Route r ($(r.i)) does not belong to a vehicle origination from the depot node nₒ ($(nₒ.i))")) end
    if !isequal(nₜ.r, r) throw(ArgumentError("Customer tail node nₜ ($(nₜ.i)) does not belong to route r ($(r.i)")) end
    if !isequal(nₕ.r, r) throw(ArgumentError("Customer head node nₕ ($(nₕ.i)) does not belong to route r ($(r.i)")) end
    if !isequal(nₜ.h, nₕ.i) throw(ArgumentError("Head index of customer tail node nₜ ($(nₜ.i)) does not match the index of customer head node nₕ ($(nₕ.i)) ")) end
    if !isequal(nₕ.t, nₜ.i) throw(ArgumentError("Tail index of customer head node nₕ ($(nₕ.i)) does not match the index of customer tail node nₜ ($(nₜ.i)) ")) end
    A = s.A
    V = s.V
    v = V[r.o]
    aₒ = A[(nₜ.i, nₕ.i)]
    aₜ = A[(nₜ.i, nₒ.i)]
    aₕ = A[(nₒ.i, nₕ.i)]
    Δˡ = aₜ.l + aₕ.l - aₒ.l
    Δᵗ = aₜ.t + aₕ.t - aₒ.t
    Δᶠ = aₜ.f + aₕ.f - aₒ.f
    Δᶜ = v.πᵐ * Δˡ + v.πʷ * Δᵗ + v.πᶠ * Δᶠ
    r.s  = nₕ.i
    r.e  = nₜ.i
    nₜ.h = nₒ.i
    nₕ.t = nₒ.i
    r.l += Δˡ
    r.t += Δᵗ
    r.f += Δᶠ
    r.c += Δᶜ
    return s
end

# REMOVE OPERATION
# Remove customer node nₒ between depot tail node nₜ and depot head node nₕ in route r
function removenode!(nₒ::CustomerNode, nₜ::DepotNode, nₕ::DepotNode, r::Route, s::Solution)
    if isopen(nₒ) throw(ArgumentError("Customer node nₒ ($(nₒ.i)) is open")) end
    if !isequal(nₒ.r, r) throw(ArgumentError("Customer node nₒ ($(nₒ.i)) does not belong to route r ($(r.i))")) end
    if s.V[r.o] ∉ nₜ.V throw(ArgumentError("Route r ($(r.i)) does not belong to a vehicle origination from the depot tail node nₜ ($(nₜ.i))")) end
    if s.V[r.o] ∉ nₕ.V throw(ArgumentError("Route r ($(r.i)) does not belong to a vehicle origination from the depot head node nₕ ($(nₕ.i))")) end
    if !isequal(nₒ.t, nₜ.i) throw(ArgumentError("Tail index of customer node nₒ ($(nₒ.i)) does not match the index of depot tail node nₜ ($(nₜ.i))")) end
    if !isequal(nₒ.h, nₕ.i) throw(ArgumentError("Head index of customer node nₒ ($(nₒ.i)) does not match the index of depot head node nₕ ($(nₕ.i))")) end
    if !isequal(r.s, nₒ.i) throw(ArgumentError("Head index of depot tail node nₜ ($(nₜ.i)) does not match the index of the customer node nₒ ($(nₒ.i))")) end
    if !isequal(r.e, nₒ.i) throw(ArgumentError("Tail index of depot head node nₕ ($(nₕ.i)) does not match the index of the customer node nₒ ($(nₒ.i))")) end
    A = s.A
    V = s.V
    v = V[r.o]
    aₒ = A[(nₜ.i, nₕ.i)]
    aₜ = A[(nₜ.i, nₒ.i)]
    aₕ = A[(nₒ.i, nₕ.i)]
    Δˡ = aₜ.l + aₕ.l - aₒ.l
    Δᵗ = aₜ.t + aₕ.t - aₒ.t
    Δᶠ = aₜ.f + aₕ.f - aₒ.f
    Δᶜ = v.πᵐ * Δˡ + v.πʷ * Δᵗ + v.πᶠ * Δᶠ
    r.s = nₕ.i
    r.e = nₜ.i
    nₒ.t = 0
    nₒ.h = 0
    # nₒ.r = Route()
    r.n -= 1
    r.q -= nₒ.q
    r.l -= Δˡ
    r.t -= Δᵗ
    r.f -= Δᶠ
    r.c -= Δᶜ
    return s
end
# Remove customer node nₒ between customer tail node nₜ and customer head node nₕ in route r
function removenode!(nₒ::CustomerNode, nₜ::CustomerNode, nₕ::CustomerNode, r::Route, s::Solution)
    if isopen(nₒ) throw(ArgumentError("Customer node nₒ ($(nₒ.i)) is open")) end
    if !isequal(nₒ.r, r) throw(ArgumentError("Customer node nₒ ($(nₒ.i)) does not belong to route r ($(r.i))")) end
    if !isequal(nₜ.r, r) throw(ArgumentError("Customer tail node nₜ ($(nₜ.i)) does not belong to route r ($(r.i))")) end
    if !isequal(nₕ.r, r) throw(ArgumentError("Customer head node nₕ ($(nₕ.i)) does not belong to route r ($(r.i))")) end 
    if !isequal(nₒ.t, nₜ.i) throw(ArgumentError("Tail index of customer node nₒ ($(nₒ.i)) does not match the index of customer tail node nₜ ($(nₜ.i))")) end
    if !isequal(nₒ.h, nₕ.i) throw(ArgumentError("Head index of customer node nₒ ($(nₒ.i)) does not match the index of customer head node nₕ ($(nₕ.i))")) end
    if !isequal(nₜ.h, nₒ.i) throw(ArgumentError("Head index of customer tail node nₜ ($(nₜ.i)) does not match the index of the customer node nₒ ($(nₒ.i))")) end
    if !isequal(nₕ.t, nₒ.i) throw(ArgumentError("Tail index of customer head node nₕ ($(nₕ.i)) does not match the index of the customer node nₒ ($(nₒ.i))")) end
    A = s.A
    V = s.V
    v = V[r.o]
    aₒ = A[(nₜ.i, nₕ.i)]
    aₜ = A[(nₜ.i, nₒ.i)]
    aₕ = A[(nₒ.i, nₕ.i)]
    Δˡ = aₜ.l + aₕ.l - aₒ.l
    Δᵗ = aₜ.t + aₕ.t - aₒ.t
    Δᶠ = aₜ.f + aₕ.f - aₒ.f
    Δᶜ = v.πᵐ * Δˡ + v.πʷ * Δᵗ + v.πᶠ * Δᶠ
    nₜ.h = nₕ.i
    nₕ.t = nₜ.i
    nₒ.t = 0
    nₒ.h = 0
    # nₒ.r = Route()
    r.n -= 1
    r.q -= nₒ.q
    r.l -= Δˡ
    r.t -= Δᵗ
    r.f -= Δᶠ
    r.c -= Δᶜ
    return s
end
# Remove customer node nₒ between depot tail node nₜ and customer head node nₕ in route r
function removenode!(nₒ::CustomerNode, nₜ::DepotNode, nₕ::CustomerNode, r::Route, s::Solution)
    if isopen(nₒ) throw(ArgumentError("Customer node nₒ ($(nₒ.i)) is open")) end
    if !isequal(nₒ.r, r) throw(ArgumentError("Customer node nₒ ($(nₒ.i)) does not belong to route r ($(r.i))")) end
    if s.V[r.o] ∉ nₜ.V throw(ArgumentError("Route r ($(r.i)) does not belong to a vehicle origination from the depot tail node nₜ ($(nₜ.i))")) end
    if !isequal(nₕ.r, r) throw(ArgumentError("Customer head node nₕ ($(nₕ.i)) does not belong to route r ($(r.i))")) end 
    if !isequal(nₒ.t, nₜ.i) throw(ArgumentError("Tail index of customer node nₒ ($(nₒ.i)) does not match the index of depot tail node nₜ ($(nₜ.i))")) end
    if !isequal(nₒ.h, nₕ.i) throw(ArgumentError("Head index of customer node nₒ ($(nₒ.i)) does not match the index of customer head node nₕ ($(nₕ.i))")) end
    if !isequal(r.s, nₒ.i) throw(ArgumentError("Head index of depot tail node nₜ ($(nₜ.i)) does not match the index of the customer node nₒ ($(nₒ.i))")) end
    if !isequal(nₕ.t, nₒ.i) throw(ArgumentError("Tail index of customer head node nₕ ($(nₕ.i)) does not match the index of the customer node nₒ ($(nₒ.i))")) end
    A = s.A
    V = s.V
    v = V[r.o]
    aₒ = A[(nₜ.i, nₕ.i)]
    aₜ = A[(nₜ.i, nₒ.i)]
    aₕ = A[(nₒ.i, nₕ.i)]
    Δˡ = aₜ.l + aₕ.l - aₒ.l
    Δᵗ = aₜ.t + aₕ.t - aₒ.t
    Δᶠ = aₜ.f + aₕ.f - aₒ.f
    Δᶜ = v.πᵐ * Δˡ + v.πʷ * Δᵗ + v.πᶠ * Δᶠ
    r.s  = nₕ.i
    nₕ.t = nₜ.i
    nₒ.t = 0
    nₒ.h = 0
    # nₒ.r = Route()
    r.n -= 1
    r.q -= nₒ.q
    r.l -= Δˡ
    r.t -= Δᵗ
    r.f -= Δᶠ
    r.c -= Δᶜ
    return s
end
# Remove customer node nₒ between customer tail node nₜ and depot head node nₕ in route r
function removenode!(nₒ::CustomerNode, nₜ::CustomerNode, nₕ::DepotNode, r::Route, s::Solution)
    if isopen(nₒ) throw(ArgumentError("Customer node nₒ ($(nₒ.i)) is open")) end
    if !isequal(nₒ.r, r) throw(ArgumentError("Customer node nₒ ($(nₒ.i)) does not belong to route r ($(r.i))")) end
    if !isequal(nₜ.r, r) throw(ArgumentError("Customer tail node nₜ ($(nₜ.i)) does not belong to route r ($(r.i))")) end
    if s.V[r.o] ∉ nₕ.V throw(ArgumentError("Route r ($(r.i)) does not belong to a vehicle origination from the depot head node nₕ ($(nₕ.i))")) end
    if !isequal(nₒ.t, nₜ.i) throw(ArgumentError("Tail index of customer node nₒ ($(nₒ.i)) does not match the index of customer tail node nₜ ($(nₜ.i))")) end
    if !isequal(nₒ.h, nₕ.i) throw(ArgumentError("Head index of customer node nₒ ($(nₒ.i)) does not match the index of depot head node nₕ ($(nₕ.i))")) end
    if !isequal(nₜ.h, nₒ.i) throw(ArgumentError("Head index of customer tail node nₜ ($(nₜ.i)) does not match the index of the customer node nₒ $(nₒ.i)")) end
    if !isequal(r.e, nₒ.i) throw(ArgumentError("Tail index of depot head node nₕ ($(nₕ.i)) does not match the index of the customer node nₒ ($(nₒ.i))")) end
    A = s.A
    V = s.V
    v = V[r.o]
    aₒ = A[(nₜ.i, nₕ.i)]
    aₜ = A[(nₜ.i, nₒ.i)]
    aₕ = A[(nₒ.i, nₕ.i)]
    Δˡ = aₜ.l + aₕ.l - aₒ.l
    Δᵗ = aₜ.t + aₕ.t - aₒ.t
    Δᶠ = aₜ.f + aₕ.f - aₒ.f
    Δᶜ = v.πᵐ * Δˡ + v.πʷ * Δᵗ + v.πᶠ * Δᶠ
    nₜ.h = nₕ.i
    r.e  = nₜ.i
    nₒ.t = 0
    nₒ.h = 0
    # nₒ.r = Route()
    r.n -= 1
    r.q -= nₒ.q
    r.l -= Δˡ
    r.t -= Δᵗ
    r.f -= Δᶠ
    r.c -= Δᶜ
    return s
end
# Remove depot node nₒ between customer tail node nₜ and customer head node nₕ in route r
function removenode!(nₒ::DepotNode, nₜ::CustomerNode, nₕ::CustomerNode, r::Route, s::Solution)
    if s.V[r.o] ∉ nₒ.V throw(ArgumentError("Route r ($(r.i)) does not belong to a vehicle origination from the depot node nₒ ($(nₒ.i))")) end
    if !isequal(nₜ.r, r) throw(ArgumentError("Customer tail node nₜ ($(nₜ.i)) does not belong to route r ($(r.i))")) end
    if !isequal(nₕ.r, r) throw(ArgumentError("Customer head node nₕ ($(nₕ.i)) does not belong to route r ($(r.i))")) end 
    if !isequal(r.s, nₕ.i) throw(ArgumentError("Head index of depot node nₒ ($(nₒ.i)) does not match the index of customer head node nₕ ($(nₕ.i))")) end
    if !isequal(r.e, nₜ.i) throw(ArgumentError("Tail index of depot node nₒ ($(nₒ.i)) does not match the index of customer tail node nₜ ($(nₜ.i))")) end
    if !isequal(nₜ.h, nₒ.i) throw(ArgumentError("Head index of customer tail node nₜ ($(nₜ.i)) does not match the index of the depot node nₒ ($(nₒ.i))")) end
    if !isequal(nₕ.t, nₒ.i) throw(ArgumentError("Tail index of customer head node nₕ ($(nₕ.i)) does not match the index of the depot node nₒ ($(nₒ.i))")) end
    A = s.A
    V = s.V
    v = V[r.o]
    aₒ = A[(nₜ.i, nₕ.i)]
    aₜ = A[(nₜ.i, nₒ.i)]
    aₕ = A[(nₒ.i, nₕ.i)]
    Δˡ = aₜ.l + aₕ.l - aₒ.l
    Δᵗ = aₜ.t + aₕ.t - aₒ.t
    Δᶠ = aₜ.f + aₕ.f - aₒ.f
    Δᶜ = v.πᵐ * Δˡ + v.πʷ * Δᵗ + v.πᶠ * Δᶠ
    nₜ.h = nₕ.i
    nₕ.t = nₜ.i
    r.s  = 0
    r.e  = 0
    r.l -= Δˡ
    r.t -= Δᵗ
    r.f -= Δᶠ
    r.c -= Δᶜ
    return s
end