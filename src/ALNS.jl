"""
    ALNS([rng::AbstractRNG], s::Solution, χ::ALNSParameters)

Adaptive Large Neighborhood Search (ALNS)

Given ALNS optimization parameters `χ` and an initial solution `s`, 
ALNS returns a vector of solutions with best found solution from every 
iteration.

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
function ALNS(rng::AbstractRNG, s::Solution, χ::ALNSParameters)
    # Step 0: Pre-initialize
    k̲, k̅, k̲ₛ, k̅ₛ = χ.k̲, χ.k̅, χ.k̲ₛ, χ.k̅ₛ 
    Ψᵣ, Ψᵢ, Ψₛ = χ.Ψᵣ, χ.Ψᵢ, χ.Ψₛ
    σ₁, σ₂, σ₃ = χ.σ₁, χ.σ₂, χ.σ₃
    ω, τ, 𝜃 = χ.ω, χ.τ, χ.𝜃
    C̲, C̅ = χ.C̲, χ.C̅
    μ̲, μ̅ = χ.μ̲, χ.μ̅
    ρ  = χ.ρ
    χₒ = χ.χₒ
    # Step 1: Initialize
    S  = Solution[]
    H  = UInt64[]
    s⃰  = deepcopy(s)
    h  = hash(s)
    j̅  = k̅÷k̲
    jₛ = k̲ₛ÷k̲
    T  = ω * f(s, χₒ)/log(ℯ, 1/τ)
    wᵣ = ones(length(Ψᵣ))
    wᵢ = ones(length(Ψᵢ))
    pᵣ = zeros(length(Ψᵣ))
    pᵢ = zeros(length(Ψᵢ))
    πᵣ = zeros(length(Ψᵣ))
    πᵢ = zeros(length(Ψᵢ))
    cᵣ = zeros(Int64, length(Ψᵣ))
    cᵢ = zeros(Int64, length(Ψᵢ))
    # Step 2: Loop over segments.
    push!(S, s⃰)
    push!(H, h)
    p = Progress(k̅, desc="Computing...", color=:blue, showspeed=true)
    for j ∈ 1:j̅
        # Step 2.1: Set operator scores and count, and update operator probabilities.
        for r ∈ 1:length(Ψᵣ) pᵣ[r], πᵣ[r], cᵣ[r] = wᵣ[r]/sum(values(wᵣ)), 0., 0 end
        for i ∈ 1:length(Ψᵢ) pᵢ[i], πᵢ[i], cᵢ[i] = wᵢ[i]/sum(values(wᵢ)), 0., 0 end
        # Step 2.2: Loop over iterations.
        for k ∈ 1:k̲
            # Step 2.2.1: Select removal and insertion operators using roulette wheel selection and update operator counts.
            r = sample(rng, 1:length(Ψᵣ), Weights(pᵣ))
            i = sample(rng, 1:length(Ψᵢ), Weights(pᵢ))
            R = Ψᵣ[r]
            I = Ψᵢ[i]
            cᵣ[r] += 1
            cᵢ[i] += 1
            # Step 2.2.2: Using the removal and insertion operators create new solution.
            η = rand(rng)
            q = Int64(floor(((1 - η) * min(C̲, μ̲ * length(s.C)) + η * min(C̅, μ̅ * length(s.C)))))
            s′= deepcopy(s)
            remove!(rng, q, s′, χₒ, R)
            insert!(rng, s′, χₒ, I)
            # Step 2.2.3: If the new solution is better than the best found then update the best and current solutions, and update the operator scores by σ₁.
            if f(s′, χₒ) < f(s⃰, χₒ)
                s = s′
                s⃰ = s
                h = hash(s)
                πᵣ[r] += σ₁
                πᵢ[i] += σ₂
                push!(H, h)
            # Step 2.2.4: Else if the new solution is better than current solution, update the current solution. If the new solution is also newly found then update the operator scores by σ₂.
            elseif f(s′, χₒ) < f(s, χₒ)
                s = s′
                h = hash(s)
                if h ∉ H
                    πᵣ[r] += σ₂
                    πᵢ[i] += σ₂
                end
                push!(H, h)
            # Step 2.2.5: Else accept the new solution with simulated annealing acceptance criterion. Further, if the new solution is also newly found then update operator scores by σ₃.
            else
                η = rand(rng)
                pr = exp(-(f(s′, χₒ) - f(s, χₒ))/T)
                if η > pr
                    s = s′
                    h = hash(s)
                    if h ∉ H
                        πᵣ[r] += σ₃
                        πᵢ[i] += σ₃
                    end
                    push!(H, h)
                end
            end
            # Step 2.2.6: Update annealing tempertature.
            T *= 𝜃
            # Step 2.2.7: Miscellaneous
            push!(S, s⃰)
            next!(p)
        end
        # Step 2.3: Update operator weights.
        for r ∈ 1:length(Ψᵣ) if !iszero(cᵣ[r]) wᵣ[r] = ρ * πᵣ[r] / cᵣ[r] + (1 - ρ) * wᵣ[r] end end
        for i ∈ 1:length(Ψᵢ) if !iszero(cᵢ[i]) wᵢ[i] = ρ * πᵢ[i] / cᵢ[i] + (1 - ρ) * wᵢ[i] end end
        # Step 2.4: Local search.
        if iszero(j%jₛ)
            for ls ∈ Ψₛ localsearch!(rng, k̅ₛ, s, χₒ, ls) end
            h = hash(s)
            if f(s, χₒ) < f(s⃰, χₒ)
                s⃰ = s
                push!(S, s⃰) 
            end
            push!(H, h)
        end
    end
    # Step 3: Return vector of solutions
    return S
end
ALNS(s::Solution, χ::ALNSParameters) = ALNS(Random.GLOBAL_RNG, s, χ)