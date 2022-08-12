"""
    vectorize(s::Solution)

Returns solution as a sequence of nodes in the order of visits.
"""
function vectorize(s::Solution)
    D = s.D
    C = s.C
    V = [Int64[] for _ ∈ D]
    for d ∈ D
        if isclose(d) continue end
        for v ∈ d.V
            for r ∈ v.R
                if isclose(r) continue end
                cₛ, cₑ = C[r.s], C[r.e]
                push!(V[d.i], d.i)
                c = cₛ
                while true
                    push!(V[d.i], c.i)
                    if isequal(c, cₑ) break end
                    c = C[c.h]
                end
            end
        end
        push!(V[d.i], d.i)
    end
    return V
end

"""
    visualize(s::Solution)

Plots solution depicting route and unvisited nodes.
"""
function visualize(s::Solution)
    D = s.D
    C = s.C
    fig = plot(legend=:none)
    # Open nodes
    V = vectorize(s)
    for d ∈ s.D
        Z = V[d.i]
        K = length(Z)
        X = zeros(Float64, K)
        Y = zeros(Float64, K)
        W = fill("color", K)
        for k ∈ 1:K
            i = Z[k]
            n = i ≤ length(D) ? D[i] : C[i]
            X[k] = n.x
            Y[k] = n.y
            if isdepot(n) W[k] = "DarkRed"
            else W[k] = "DarkBlue"
            end
        end
        scatter!(X, Y, markersize=5, markerstrokewidth=0, color=W)
        plot!(X, Y, color="SteelBlue")
    end
    # Closed nodes
    L  = [c.i for c ∈ C if isopen(c)]
    for d ∈ D if isclose(d) push!(L, d.i) end end
    Z′ = L
    K′ = length(Z′)
    X′ = zeros(Float64, K′)
    Y′ = zeros(Float64, K′)
    W′ = fill("color", K′)
    for k ∈ 1:K′
        i = Z′[k]
        n = i ≤ length(D) ? D[i] : C[i]
        X′[k] = n.x
        Y′[k] = n.y
        if isdepot(n) W′[k] = "IndianRed"
        else W′[k] = "LightBlue"
        end
    end
    scatter!(X′, Y′, markersize=5, markerstrokewidth=0, color=W′)
    return fig
end

"""
    animate(S::Vector{Solution})

Iteratively plots solutions in `S` to develop a gif.
"""
function animate(S::Vector{Solution})
    K = 0:(length(S)-1)
    figs = Vector{Plots.Plot{Plots.GRBackend}}(undef, length(S))
    for (k, s) ∈ enumerate(S)
        fig = visualize(s)
        plot!(title="iter:$(K[k])")
        figs[k] = fig
    end
    anim =  @animate for fig in figs
        plot(fig)
    end
    gif(anim, fps = 10, show_msg=false)
end

"""
    convergence(S::Vector{Solution}, χₒ::ObjectiveFunctionParameters)

Plots objective function values for solutions in `S` using objective function parameters `χₒ`.
"""
function convergence(S::Vector{Solution}, χₒ::ObjectiveFunctionParameters)
    Y = [f(s, χₒ) for s ∈ S]
    X = 0:(length(S)-1)
    plot(X,Y)
end