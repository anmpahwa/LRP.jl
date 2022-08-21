# INSERT OPERATION
# Insert node nₒ between tail node nₜ and head node nₕ in route r
function insertnode!(nₒ::Node, nₜ::Node, nₕ::Node, r::Route, s::Solution)
    aₒ = s.A[(nₜ.i, nₕ.i)]
    aₜ = s.A[(nₜ.i, nₒ.i)]
    aₕ = s.A[(nₒ.i, nₕ.i)]
    if isdepot(nₜ) r.s = nₒ.i 
    else nₜ.h = nₒ.i
    end
    if isdepot(nₕ) r.e = nₒ.i
    else nₕ.t = nₒ.i
    end
    if isdepot(nₒ) 
        r.s = nₕ.i
        r.e = nₜ.i
    else
        nₒ.h = nₕ.i
        nₒ.t = nₜ.i
        nₒ.r = r
        r.n += 1
        r.q += nₒ.q
    end
    r.l += aₜ.l + aₕ.l - aₒ.l
    return s
end

# REMOVE OPERATION
# Remove node nₒ between tail node nₜ and head node nₕ in route r
function removenode!(nₒ::Node, nₜ::Node, nₕ::Node, r::Route, s::Solution)
    aₒ = s.A[(nₜ.i, nₕ.i)]
    aₜ = s.A[(nₜ.i, nₒ.i)]
    aₕ = s.A[(nₒ.i, nₕ.i)]
    if isdepot(nₜ) r.s = nₕ.i
    else nₜ.h = nₕ.i
    end
    if isdepot(nₕ) r.e = nₜ.i
    else nₕ.t = nₜ.i
    end
    if isdepot(nₒ) 
        r.s = 0
        r.e = 0
    else
        nₒ.h = 0
        nₒ.t = 0
        #nₒ.r = Route()
        r.n -= 1
        r.q -= nₒ.q
    end
    
    r.l -= aₜ.l + aₕ.l - aₒ.l
    return s
end