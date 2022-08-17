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
    visualize(s::Solution; backend=gr)

Plots solution depicting route and unvisited nodes (if any).
Uses given backend to plot (default backend gr).
"""
function visualize(s::Solution; backend=gr)
    backend()
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
    animate(S::Vector{Solution}, fps=10)

Iteratively plots solutions in `S` to develop a gif at given `fps`.
"""
function animate(S::Vector{Solution}, fps=10)
    K = 0:(length(S)-1)
    figs = Vector(undef, length(S))
    for (k, s) ∈ enumerate(S)
        fig = visualize(s, backend=gr)
        plot!(title="Iteration #$(K[k])", titlefontsize=11)
        figs[k] = fig
    end
    anim = @animate for fig in figs
        plot(fig)
    end
    gif(anim, fps=fps, show_msg=false)
end

"""
    convergence(S::Vector{Solution}; backend=gr)

Plots objective function values for solutions in `S`.
Uses given backend to plot (default backend gr).
"""
function convergence(S::Vector{Solution}; backend=gr)
    backend()
    Y = [f(s) for s ∈ S]
    X = 0:(length(S)-1)
    fig = plot(legend=:none)
    plot!(X,Y, xlabel="iterations", ylabel="objective function value")
    return fig
end