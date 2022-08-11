@doc """
    Node(i::Int64, x::Float64, y::Float64, t=i, h=i)

A `Node` is a point on the graph at `(x,y)` with index `i`, 
tail node index `t`, and head node index `h`.
"""
mutable struct Node
    i::Int64                                                    # Node index
    x::Float64                                                  # Node location on the x-axis
    y::Float64                                                  # Node location on the y-axis
    t::Int64                                                    # Tail node index
    h::Int64                                                    # Head node index
    Node(i, x, y) = new(i, x, y, i, i)
end
Base.isequal(p::Node, q::Node) = isequal(p.i, q.i)
isopen(n::Node) = isequal(n.t, n.i) & isequal(n.i, n.h)
isclose(n::Node) = !isopen(n)

@doc """
    Arc(i::Int64, j::Int64, c::Float64)

An `Arc` is a connection between node `i` and `j` with 
traversal cost `c`.
"""
struct Arc
    i::Int64                                                    # Tail node index
    j::Int64                                                    # Head node index
    c::Float64                                                  # Cost
end

@doc """
    Solution(N::Vector{Node}, A::Dict{Tuple{Int64,Int64}, Arc}, c=0.)

A `Solution` is a graph with nodes `N`, arcs `A`, and TSP route 
cost `c`.
"""
mutable struct Solution
    N::Vector{Node}                                             # Vector of depot nodes
    A::Dict{Tuple{Int64,Int64}, Arc}                            # Set of arcs
    c::Float64                                                  # TSP route cost
    Solution(N, A) = new(N, A, 0.)
end
"""
    f(s::Solution)

Returns objective function value (solution cost)
"""
f(s::Solution) = s.c
"""
    infeasible(s::Solution)

Returns true if node service constraint, node flow constraint, and
sub-tour elimination constraint are not violated.
"""
function isfeasible(s::Solution) 
    V = vectorize(s)
    pop!(V)
    return allunique(V)
end