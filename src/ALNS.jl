"""
    ALNS([rng::AbstractRNG], χ::ALNSparameters, sₒ::Solution)

Adaptive Large Neighborhood Search (ALNS)

Given ALNS optimization parameters `χ` and an initial solution `sₒ`, 
ALNS returns a vector of solutions with current solution from every 
iteration.

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
function ALNS(rng::AbstractRNG, χ::ALNSparameters, sₒ::Solution)
    # Step 0: Pre-initialize
    j, k = χ.j, χ.k
    n, m = χ.n, χ.m
    Ψᵣ, Ψᵢ, Ψₗ = χ.Ψᵣ, χ.Ψᵢ, χ.Ψₗ
    σ₁, σ₂, σ₃ = χ.σ₁, χ.σ₂, χ.σ₃
    μ̲, C̲ = χ.μ̲, χ.C̲
    μ̅, C̅ = χ.μ̅, χ.C̅
    ω̅, τ̅ = χ.ω̅, χ.τ̅
    ω̲, τ̲ = χ.ω̲, χ.τ̲
    φ, θ, ρ = χ.φ, χ.θ, χ.ρ   
    R = eachindex(Ψᵣ)
    I = eachindex(Ψᵢ)
    L = eachindex(Ψₗ)
    H = OffsetVector{UInt64}(undef, 0:j*n)
    S = OffsetVector{Solution}(undef, 0:j*n)
    # Step 1: Initialize
    s = deepcopy(sₒ)
    s⃰ = s
    z = f(sₒ)
    z⃰ = z
    h = hash(s)
    S[0] = s
    H[0] = h
    T = ω̅ * z⃰/log(1/τ̅)
    cᵣ, pᵣ, πᵣ, wᵣ = zeros(Int64, R), zeros(R), zeros(R), ones(R)
    cᵢ, pᵢ, πᵢ, wᵢ = zeros(Int64, I), zeros(I), zeros(I), ones(I)
    # Step 2: Loop over segments.
    p = Progress(n * j, desc="Computing...", color=:blue, showspeed=true)
    for u ∈ 1:j
        # Step 2.1: Reset count and score for every removal and insertion operator
        for r ∈ R cᵣ[r], πᵣ[r] = 0, 0. end
        for i ∈ I cᵢ[i], πᵢ[i] = 0, 0. end
        # Step 2.2: Update selection probability for every removal and insertion operator
        for r ∈ R pᵣ[r] = wᵣ[r]/sum(values(wᵣ)) end
        for i ∈ I pᵢ[i] = wᵢ[i]/sum(values(wᵢ)) end
        # Step 2.3: Loop over iterations within the segment
        for v ∈ 1:n
            # Step 2.3.1: Randomly select a removal and an insertion operator based on operator selection probabilities, and consequently update count for the selected operators.
            r = sample(rng, 1:length(Ψᵣ), Weights(pᵣ))
            i = sample(rng, 1:length(Ψᵢ), Weights(pᵢ))
            cᵣ[r] += 1
            cᵢ[i] += 1
            # Step 2.3.2: Using the selected removal and insertion operators destroy and repair the current solution to develop a new solution.
            η = rand(rng)
            q = Int64(floor(((1 - η) * min(C̲, μ̲ * length(s.C)) + η * min(C̅, μ̅ * length(s.C)))))
            s′= deepcopy(s)
            remove!(rng, q, s′, Ψᵣ[r])
            insert!(rng, s′, Ψᵢ[i])
            z′ = f(s′)
            h  = hash(s′)
            # Step 2.3.3: If this new solution is better than the best solution, then set the best solution and the current solution to the new solution, and accordingly update scores of the selected removal and insertion operators by σ₁.
            if z′ < z⃰
                s = s′
                s⃰ = s′
                z = z′
                z⃰ = z′
                πᵣ[r] += σ₁
                πᵢ[i] += σ₂
            # Step 2.3.4: Else if this new solution is only better than the current solution, then set the current solution to the new solution and accordingly update scores of the selected removal and insertion operators by σ₂.
            elseif z′ < z
                s = s′
                z = z′
                if h ∉ H
                    πᵣ[r] += σ₂
                    πᵢ[i] += σ₂
                end
            # Step 2.3.5: Else accept the new solution with simulated annealing acceptance criterion. Further, if the new solution is also newly found then update operator scores by σ₃.
            else
                η = rand(rng)
                pr = exp(-(z′ - z)/T)
                if η < pr
                    s = s′
                    z = z′
                    if h ∉ H
                        πᵣ[r] += σ₃
                        πᵢ[i] += σ₃
                    end
                end
            end
            S[(u - 1) * n + v] = s
            H[(u - 1) * n + v] = h
            T = max(T * θ, ω̲ * z⃰/log(1/τ̲))
            next!(p)
        end
        # Step 2.4: Update weights for every removal and insertion operator.
        for r ∈ R if !iszero(cᵣ[r]) wᵣ[r] = ρ * πᵣ[r] / cᵣ[r] + (1 - ρ) * wᵣ[r] end end
        for i ∈ I if !iszero(cᵢ[i]) wᵢ[i] = ρ * πᵢ[i] / cᵢ[i] + (1 - ρ) * wᵢ[i] end end
        # Step 2.5: Reset current solution.
        if iszero(u % k) s, z = deepcopy(s⃰), z⃰ end
        # Step 2.6: Perform local search.
        if z ≤ z⃰ * (1 + φ)
            s′ = deepcopy(s)
            for l ∈ L localsearch!(rng, m, s′, Ψₗ[l]) end
            z′ = f(s′)
            if z′ < z⃰
                s = s′
                s⃰ = s′
                z = z′
                z⃰ = z′
            elseif z′ < z
                s = s′
                z = z′
            end
        end
    end
    # Step 3: Return vector of solutions
    return S
end
ALNS(χ::ALNSparameters, s::Solution) = ALNS(Random.GLOBAL_RNG, χ, s)