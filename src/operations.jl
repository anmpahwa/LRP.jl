"""
    insertnode!(nₒ::Node, nₜ::Node, nₕ::Node, s::Solution)

Insert node `nₒ` between tail node `nₜ` and head node `nₕ` in solution `s`.
"""
function insertnode!(nₒ::Node, nₜ::Node, nₕ::Node, s::Solution)
    A  = s.A 
    aₒ = A[(nₜ.i, nₕ.i)]
    aₜ = A[(nₜ.i, nₒ.i)]
    aₕ = A[(nₒ.i, nₕ.i)]
    nₜ.h = nₒ.i
    nₕ.t = nₒ.i
    nₒ.t = nₜ.i
    nₒ.h = nₕ.i
    s.c += aₜ.c + aₕ.c - aₒ.c
    return
end

"""
    removenode!(nₒ::Node, nₜ::Node, nₕ::Node, s::Solution)

Remove node `nₒ` from its position between tail node `nₜ` and head node `nₕ` in solution `s`.
"""
function removenode!(nₒ::Node, nₜ::Node, nₕ::Node, s::Solution)
    A  = s.A 
    aₒ = A[(nₜ.i, nₕ.i)]
    aₜ = A[(nₜ.i, nₒ.i)]
    aₕ = A[(nₒ.i, nₕ.i)]
    nₜ.h = nₒ.h
    nₕ.t = nₒ.t
    nₒ.t = nₒ.i
    nₒ.h = nₒ.i
    s.c -= aₜ.c + aₕ.c - aₒ.c
    return
end