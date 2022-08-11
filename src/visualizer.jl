"""
    vectorize(s::Solution)

Returns solution as a sequence of nodes in the order of visits.
"""
function vectorize(s::Solution)
    N = s.N
    d = N[1]
    V = Int64[]
    if isopen(d) return V end
    nₜ = d
    nₕ = N[nₜ.h]
    push!(V, d.i)
    while true
        push!(V, nₕ.i)
        if isequal(nₕ, d) break end
        nₜ = nₕ
        nₕ = N[nₜ.h]
    end
    return V
end

"""
    visualize(s::Solution)

Plots solution depicting route and unvisited nodes (if any).
"""
function visualize(s::Solution)
    N = s.N
    fig = plot(legend=:none)
    # Open nodes
    V = vectorize(s)
    Z = V
    K = length(Z)
    X = zeros(Float64, K)
    Y = zeros(Float64, K)
    W = fill("color", K)
    for k ∈ 1:K
        i = Z[k]
        n = N[i]
        X[k] = n.x
        Y[k] = n.y
        W[k] = "DarkBlue"
    end
    scatter!(X, Y, markersize=5, markerstrokewidth=0, color=W)
    plot!(X, Y, color="SteelBlue")
    # Closed nodes
    L  = [n.i for n ∈ N if isopen(n)]
    Z′ = L
    K′ = length(Z′)
    X′ = zeros(Float64, K′)
    Y′ = zeros(Float64, K′)
    W′ = fill("color", K′)
    for k ∈ 1:K′
        i = Z′[k]
        n = N[i]
        X′[k] = n.x
        Y′[k] = n.y
        W′[k] = "LightBlue"
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
    figs = Vector{Plots.Plot{Plots.GRBackend}}(undef, length(S))
    for (k, s) ∈ enumerate(S)
        fig = visualize(s)
        plot!(title="iter:$(K[k])")
        figs[k] = fig
    end
    anim = @animate for fig in figs
        plot(fig)
    end
    gif(anim, fps=fps, show_msg=false)
end

"""
    convergence(S::Vector{Solution})

Plots objective function values for solutions in `S`.
"""
function convergence(S::Vector{Solution})
    X = [f(s) for s ∈ S]
    Y = 0:(length(S)-1)
    plot(X,Y)
end