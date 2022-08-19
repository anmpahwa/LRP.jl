"""
    ALNS([rng::AbstractRNG], χ::ALNSParameters, s::Solution)

Adaptive Large Neighborhood Search (ALNS)

Given ALNS optimization parameters `χ` and an initial solution `s`, 
ALNS returns a vector of solutions with best found solution from every 
iteration.

Optionally specify a random number generator `rng` as the first argument
(defaults to `Random.GLOBAL_RNG`).
"""
function ALNS(rng::AbstractRNG, χ::ALNSParameters, s::Solution)
    # Step 0: Pre-initialize
    k̲, l̲, l̅, k̅ = χ.k̲, χ.l̲, χ.l̅, χ.k̅
    Ψᵣ, Ψᵢ, Ψₗ = χ.Ψᵣ, χ.Ψᵢ, χ.Ψₗ
    σ₁, σ₂, σ₃ = χ.σ₁, χ.σ₂, χ.σ₃
    ω, τ, 𝜃 = χ.ω, χ.τ, χ.𝜃
    C̲, C̅ = χ.C̲, χ.C̅
    μ̲, μ̅ = χ.μ̲, χ.μ̅
    ρ  = χ.ρ
    # Step 1: Initialize
    S  = Solution[]
    H  = UInt64[]
    s⃰  = deepcopy(s)
    h  = hash(s)
    j̅  = k̅ ÷ k̲
    jₗ = l̲ ÷ k̲
    T  = ω * f(s)/log(ℯ, 1/τ)
    R = eachindex(Ψᵣ)
    I = eachindex(Ψᵢ)
    L = eachindex(Ψᵢ)
    wᵣ = ones(R)
    wᵢ = ones(I)
    pᵣ = zeros(R)
    pᵢ = zeros(I)
    πᵣ = zeros(R)
    πᵢ = zeros(I)
    cᵣ = zeros(Int64, R)
    cᵢ = zeros(Int64, I)
    # Step 2: Loop over segments.
    push!(S, s⃰)
    push!(H, h)
    p = Progress(k̅, desc="Computing...", color=:blue, showspeed=true)
    for j ∈ 1:j̅
        # Step 2.1: Set operator scores and count.
        for r ∈ R pᵣ[r], πᵣ[r], cᵣ[r] = wᵣ[r]/sum(values(wᵣ)), 0., 0 end
        for i ∈ I pᵢ[i], πᵢ[i], cᵢ[i] = wᵢ[i]/sum(values(wᵢ)), 0., 0 end
        # Step 2.2: Update operator probabilities
        for r ∈ R pᵣ[r] = wᵣ[r]/sum(values(wᵣ)) end
        for i ∈ I pᵢ[i] = wᵢ[i]/sum(values(wᵢ)) end
        # Step 2.2: Loop over iterations.
        for k ∈ 1:k̲
            # Step 2.2.1: Select removal and insertion operators using roulette wheel selection and update operator counts.
            r = sample(rng, 1:length(Ψᵣ), Weights(pᵣ))
            i = sample(rng, 1:length(Ψᵢ), Weights(pᵢ))
            cᵣ[r] += 1
            cᵢ[i] += 1
            # Step 2.3.2: Using the removal and insertion operators create new solution.
            η = rand(rng)
            q = Int64(floor(((1 - η) * min(C̲, μ̲ * length(s.C)) + η * min(C̅, μ̅ * length(s.C)))))
            s′= deepcopy(s)
            remove!(rng, q, s′, Ψᵣ[r])
            insert!(rng, s′, Ψᵢ[i])
            # Step 2.3.3: If the new solution is better than the best found then update the best and current solutions, and update the operator scores by σ₁.
            if f(s′) < f(s⃰)
                s = s′
                s⃰ = s
                h = hash(s)
                πᵣ[r] += σ₁
                πᵢ[i] += σ₂
                push!(H, h)
            # Step 2.3.4: Else if the new solution is better than current solution, update the current solution. If the new solution is also newly found then update the operator scores by σ₂.
            elseif f(s′) < f(s)
                s = s′
                h = hash(s)
                if h ∉ H
                    πᵣ[r] += σ₂
                    πᵢ[i] += σ₂
                end
                push!(H, h)
            # Step 2.3.5: Else accept the new solution with simulated annealing acceptance criterion. Further, if the new solution is also newly found then update operator scores by σ₃.
            else
                η = rand(rng)
                pr = exp(-(f(s′) - f(s))/T)
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
            # Step 2.3.6: Update annealing tempertature.
            T *= 𝜃
            # Step 2.3.7: Miscellaneous
            push!(S, s⃰)
            next!(p)
        end
        # Step 2.4: Update operator weights.
        for r ∈ R if !iszero(cᵣ[r]) wᵣ[r] = ρ * πᵣ[r] / cᵣ[r] + (1 - ρ) * wᵣ[r] end end
        for i ∈ I if !iszero(cᵢ[i]) wᵢ[i] = ρ * πᵢ[i] / cᵢ[i] + (1 - ρ) * wᵢ[i] end end
        # Step 2.5: Local search.
        if iszero(j%jₗ)
            for l ∈ L localsearch!(rng, l̅, s, Ψₗ[l]) end
            h = hash(s)
            if f(s) < f(s⃰)
                s⃰ = s
                push!(S, s⃰) 
            end
            push!(H, h)
        end
    end
    # Step 3: Return vector of solutions
    return S
end
ALNS(χ::ALNSParameters, s::Solution) = ALNS(Random.GLOBAL_RNG, χ, s)