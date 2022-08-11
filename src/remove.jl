# NOTE: The first node, assumed to be the depot node, is prohibited from being removed.

"""
    remove!([rng], q, s::Solution, method::Symbol)

Return solution removing `q` nodes from solution `s` using the given `method`.

Available methods include,
- Random Removal    : `random_remove!`
- Worst Removal     : `worst_remove!`
- Shaw Removal      : `shaw_remove!`

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
remove!(rng::AbstractRNG, q, s::Solution, method::Symbol)::Solution = getfield(LRP, method)(rng, q, s)
remove!(q, s::Solution, method::Symbol) = remove!(Random.GLOBAL_RNG, q, s, method)

# Random Removal
# Randomly select q nodes to remove
function random_remove!(rng::AbstractRNG, q, s::Solution)
    N = s.N
    I = length(N)
    d = N[1]
    w = [if isequal(n, d) 0 else 1 end for n ∈ N]   # w[i]: selection weight of node N[i]
    # Step 1: Randomly select customer nodes to remove until q nodes have been removed
    for _ ∈ 1:q
        i  = sample(rng, 1:I, Weights(w))
        nₒ = N[i]
        nₜ = N[nₒ.t]
        nₕ = N[nₒ.h]
        removenode!(nₒ, nₜ, nₕ, s)
        w[i] = 0
    end
    # Step 2: Return solution
    return s
end

# Shaw Removal
# For a randomly selected customer node, remove q most related customer nodes
function shaw_remove!(rng::AbstractRNG, q, s::Solution)
    N = s.N
    A = s.A
    I = length(N)-1
    # Step 1: Randomly select a pivot customer node
    j = rand(rng, 2:I)
    # Step 2: For each customer node, evaluate relatedness to this pivot customer node
    x = fill(-Inf, I)                               # x[i]: relatedness of node N[i] with node N[j]
    for i ∈ 2:I
        a = A[(i,j)]
        x[i] = 1/a.c
    end
    # Step 3: Remove q most related customer nodes
    for _ ∈ 1:q
        i  = argmax(x)
        nₒ = N[i]
        nₜ = N[nₒ.t]
        nₕ = N[nₒ.h]
        removenode!(nₒ, nₜ, nₕ, s)
        x[i] = -Inf
    end
    # Step 4: Return solution
    return s
end

# Worst Removal
# Remove q nodes with highest removal cost
function worst_remove!(rng::AbstractRNG, q, s::Solution)
    N = s.N
    I = length(N)
    # Step 1: Iterate until q nodes have been removed
    x = fill(-Inf, I)                               # x[i]: removal cost of node N[i]
    for _ ∈ 1:q
        # Step 1.1: For every closed node evaluate removal cost
        z = f(s)
        for i ∈ 2:I
            nₒ = N[i]
            if isopen(nₒ) continue end
            # Step 1.1.1: Remove closed node nₒ between tail node nₜ and head node nₕ
            nₜ = N[nₒ.t]
            nₕ = N[nₒ.h]
            removenode!(nₒ, nₜ, nₕ, s)
            # Step 1.1.2: Evaluate the removal cost
            z⁻ = f(s)
            x[i] = z - z⁻
            # Step 1.1.3: Re-insert node nₒ between tail node nₜ and head node nₕ
            insertnode!(nₒ, nₜ, nₕ, s)
        end
        # Step 1.2: Remove the node with highest removal cost
        i  = argmax(x)
        nₒ = N[i]
        nₜ = N[nₒ.t]
        nₕ = N[nₒ.h]
        removenode!(nₒ, nₜ, nₕ, s)
        # Step 1.3: Update cost vector
        x[i] = -Inf
    end
    # Step 2: Return solution
    return s
end