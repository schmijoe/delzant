using LinearAlgebra
using GLMakie

Base.rationalize(x::Rational) = x
Base.rationalize(x::Integer) = x//1

struct Edge
  λ::Vector{<:Integer}
  c::Rational{Int128}
  """
    `Edge(λ::Vector{<:Real}, c::<:Real)`

  Construct edge/halfspace satisfying equation `λ'*x + c >= 0`.
  """
  function Edge(λ::Vector{<:Real}, c::T) where T<:Real
    λ = rationalize.(λ); c = rationalize(c)
    d = gcd(λ...);
    new(λ .÷ d, c//d)
  end
end
"""
  `Edge(p1::Vector{<:Real}, p2::Vector{<:Real})`

Consturct edge/hafspace containing points `p1` and `p2`.
"""
function Edge(p1::Vector{<:Real}, p2::Vector{<:Real})
  λ = [0 -1; 1 0] * (p2-p1)
  Edge(λ, -λ'*p1)
end
Base.hash(e::Edge, h::UInt) = hash((e.λ,e.c),h)
Base.isequal(e1::Edge,e2::Edge) = (e1.λ == e2.λ && e1.c == e2.c)
Base.:(==)(e1::Edge,e2::Edge) = Base.isequal(e1,e2)

ev(e::Edge, p::Vector{<:Real}) = p'*e.λ + e.c


"""
  `intersect(e1::Edge, e2::Edge)`
Intersection of `e1`, `e2` with orientation.
"""
function intersect(e1::Edge, e2::Edge)
  A = [e1.λ'//1; e2.λ'//1]; detA = det(A);
  (detA == 0 ? nothing : -A^-1*[e1.c; e2.c], sign(detA))
end
e1::Edge ∩ e2::Edge = intersect(e1,e2)


struct Polygon
  edges::Vector{<:Edge}
  vertices::Vector{Vector{<:Rational}}
end

"""
Construct the polygon bounded by `edges`
"""
function Polygon(edges::Edge...)
  edges = unique(edges) # remove duplicate edges
  intersections = []
  # Collect all intersections
  for (n,e1) in enumerate(edges) # iterate over all unordered pairs
    for e2 in Iterators.drop(edges,n)
      i,orientation = e1 ∩ e2
      if !isnothing(i)
        # intersections are labelled by their edges, in counter-clockwise order
        push!(intersections,
              (i, orientation==1 ? [e1,e2] : [e2,e1])
             )
      end
    end
  end

  # discard intersections outside the polygon
  for e in edges
    filter!(i->(ev(e,i[1]) >= 0),intersections)
  end

  # filter out redundant edges
  uint = unique(first.(intersections))
  redundant_edges = [e for e in edges if count(int->(ev(e,int)==0), uint) < 2]
  for e in redundant_edges
    filter!(i->!(e in i[2]), intersections)
  end

  if isempty(intersections)
    @warn "Polygon was empty!" edges
    return nothing
  end

  # Order intersections by outgoing edge orientation
  sort!(intersections, by=i->(atan(i[2][2].λ[2],i[2][2].λ[1])))
  vertices=first.(intersections)
  edges=[int[2][2] for int in intersections]
  Polygon(edges, vertices)
end

"""
Construct polygon given by `vertices`. They must be in counter-clockwise order.

If `ensure_convex` is true, the resulting edges are passed to `Polygon(edges::Edge...)`
"""
function Polygon(vertices::Vector{<:Real}...; ensure_convex=false)
  vertices = [rationalize.(v) for v in vertices]
  edges = Edge[]
  for i in eachindex(vertices)
    push!(edges, Edge(vertices[i], vertices[i%end+1]))
  end
  ensure_convex ? Polygon(edges...) : Polygon(edges, vertices)
end

Base.:*(M::Matrix{<:Real}, Δ::Polygon) = Polygon( (det(M) > 0 ? [M * v for v in Δ.vertices] : [M * v for v in reverse(Δ.vertices)] )...)
Base.:+(u::Vector{<:Real}, Δ::Polygon) = Polygon([u + v for v in Δ.vertices]...)

"""
  `intersect(e::Edge, Δ::Polygon)`

Intersections of `e` with `Δ`
"""
function intersect(e::Edge, Δ::Polygon)
  intersections = []
  for e1 in Δ.edges
    int = e1 ∩ e
    if !isnothing(int[1])
      push!(intersections,int)
    end
  end
  for e1 in Δ.edges
    filter!(i->(ev(e1,i[1])>=0), intersections)
  end
  unique!(i->i[1],intersections)
  sort!(intersections, by=i->i[2])
  first.(intersections)
end
e::Edge ∩ Δ::Polygon = intersect(e,Δ)


# Enable Makie to plot Polygon
Makie.convert_arguments(P::PointBased, Δ::Polygon) = convert_arguments(P, Point2.(Δ.vertices))

"""
Get some possible probe direction from `e`
"""
get_default_probe(e::Edge) = collect(gcdx(e.λ...)[2:3])

"""
  `get_probe_intervall(e::Vector{Edge}, probe::Vector{<:Integer})`

Get range of `k` for which `probe + k*[0 -1; 1 0]*e[2].λ` is a reasonable probe in the tripple of edges `e[1]`, `e[2]` ,`e[3]`.

"""
function get_probe_intervall(e::Vector{Edge}, probe::Vector{<:Integer})
  v = [0 -1; 1 0] * e[2].λ
  a = - (probe' * e[1].λ) / (v' * e[1].λ)
  b = - (probe' * e[3].λ) / (v' * e[3].λ)
  u = sort([a,b])
  return floor(Int,u[1]):ceil(Int,u[2])
end

"""
  `get_probe_intervall(Δ::Polygon, i::Integer, probe::Vector{<:Integer})`

Get range of `k` for which `probe + k*[0 -1; 1 0]*Δ.edges[i].λ` is a reasonable probe in the polygon Δ.
"""
function get_probe_intervall(Δ::Polygon, i::Integer, probe::Vector{<:Integer})
  if 1<i<length(Δ.edges)
    e = Δ.edges[i-1:i+1]
  elseif i == 1
    e = Δ.edges[[end,1,2]]
  else
    e = Δ.edges[[i-1,i,1]]
  end
  get_probe_intervall(e, probe)
end

"""
  `get_probe_range(Δ::Polygon, edge_index::Integer, probe::Vector{<:Integer})`

Get a polygon describing which fibres in Δ can be displaced by a probe shot from edge `Δ.edges[edge_index]` in direction `probe`.
"""
function get_probe_range(Δ::Polygon, edge_index::Integer, probe::Vector{<:Integer})::Polygon
  l = Δ.vertices[edge_index]
  r = Δ.vertices[edge_index%end + 1]
  e = Δ.edges[edge_index]
  λ = [0 -1;1 0]*probe
  range = Polygon(Edge(λ,-λ'*r), Edge(-λ, λ'*l), Δ.edges...)
  M = [[0 -1;1 0]*e.λ//1 probe//1]
  offset = - e.λ * e.c//sum(t->t^2, e.λ)
  offset + (M * [1 0; 0 1//2] * M^(-1)) * (-offset + range)
end

"""
  `get_probe_ranges(Δ::Polygon, edge_index::Integer)`

Get a list of polygons describing which fibres in `Δ` can be displaced by shooting reasonable probes from `Δ.edges[edge_index]`.
"""
function get_probe_ranges(Δ::Polygon, edge_index::Integer)::Vector{Polygon}
  e = Δ.edges[edge_index]
  ed = [0 -1; 1 0] * e.λ
  dp = get_default_probe(e)
  ranges = Polygon[]
  for i in get_probe_intervall(Δ, edge_index, dp)
    probe = dp + i*ed
    push!(ranges, get_probe_range(Δ, edge_index, probe))
  end
  ranges
end

"""
  `get_probe_ranges(Δ::Polygon)`

Get a list for every edge `e` of `Δ` of polygons describing which fibres in `Δ` can be displaced by shooting reasonable probes from `e`.
"""
get_probe_ranges(Δ::Polygon)::Vector{Vector{Polygon}} = [get_probe_ranges(Δ,i) for i in eachindex(Δ.edges)]

"""
  `get_branch_cut_line(e1::Edge, e2::Edge)`

Get an `Edge` representing the branch cut line if there was an almost toric corner at the intersection of `e1` and `e2`, aswell as the number of nodes.
"""
function get_branch_cut_line(e1::Edge, e2::Edge)
  v, orientation = e1 ∩ e2
  if isnothing(v)
    return (nothing, nothing)
  end
  if orientation == -1
    e1,e2 = e2,e1
  end
  λ = e2.λ - e1.λ; kdist = gcd(λ...); λ = λ .÷ kdist;
  dist = - λ[1]*e1.λ[2] + λ[2]*e1.λ[1]
  if kdist % dist != 0
    @warn "I don't think this is a potential almost toric corner"
  end
  (Edge(λ, - λ'*v), kdist ÷ dist)
end
"""
  `get_branch_cut_line(Δ::Polygon, vertex_index::Integer)`

Get an `Edge` representing the branch cut line if there was an almost toric corner at the vertex `Δ.vertices[vertex_index]`, aswell as the number of nodes.
"""
get_branch_cut_line(Δ::Polygon, vertex_index::Integer) = 
  get_branch_cut_line(vertex_index == 1 ? last(Δ.edges) : Δ.edges[vertex_index - 1],
                      Δ.edges[vertex_index])

"""
  `mutate(Δ::Polygon, branch_cut_line::Tuple{Edge, <:Integer})`

Get a polygon obtained by appling a shear matrix `M` corresponding to `branch_cut_line = (e, k)` to the positive side of `Δ` wrt. `e`, where `k` determines the power of `M`. If `k` is negative it will instead be applied to the negative side of `Δ` wrt. `e`.
"""
function mutate(Δ::Polygon, branch_cut_line::Tuple{Edge, <:Integer})
  # 1. Insert missing vertices at intersections of branch_cut_line with polygon:
  intersections = []
  for (i,e) in enumerate(Δ.edges) # collect intersections
    int,_ = e ∩ branch_cut_line[1]
    if !isnothing(int)
      push!(intersections, (int,i))
    end
  end
  for e in Δ.edges # discard intersections outside
    filter!(i->(ev(e,i[1]) >= 0),intersections)
  end
  offset = intersections[1][1]
  uint = unique(i->(i[1]), intersections)
  intersections = filter!(int1->(count(int2->(int2[1] == int1[1]), intersections) == 1),
                          uint) # select intersections that occur exactly once
  sort!(intersections, by=i->(i[2]), rev=true) # reverse order by the edge index
  vertices = copy(Δ.vertices)
  for (int, edge_index) in intersections
    insert!(vertices, edge_index+1, int)
  end

  # 2. Apply mutation
  p,q = [0 1;-1 0]*branch_cut_line[1].λ; k=branch_cut_line[2]
  M = [1+k*p*q -k*p^2; k*q^2 1-k*p*q]
  for i in eachindex(vertices)
    l = ev(branch_cut_line[1], vertices[i])
    if l*k > 0
      vertices[i] = offset + M*(vertices[i]-offset)
    end
  end
  Polygon(vertices..., ensure_convex=true)
end

"""
  `mutate(Δ::Polygon, vertex_index::Integer, k::Integer=1)`

Assume the corner `Δ.vertices[vertex_index]` is almost toric, and try to perform `k` mutations.
"""
function mutate(Δ::Polygon, vertex_index::Integer, k::Integer=1)
  bcl = get_branch_cut_line(Δ, vertex_index)
  if abs(k) > bcl[2]
    @warn "k is too large." k bcl
  end
  k = (sign(atan(bcl[1].λ[2],bcl[1].λ[1])) == 1) ? k : -k
  mutate(Δ, (bcl[1],k))
end


function interact(Δ::Polygon; button_size = 30, button_color = (:black, 0.1))
  fig = Figure()
  ax = Axis(fig[1,1], autolimitaspect = 1.0)

  Δ = Observable(Δ)

  probe_range_colors = Integer[]
  probe_ranges = lift(get_probe_ranges, Δ)
  probe_range_colors = @lift(Integer[i for (i,edge_ranges) in enumerate($probe_ranges) for _ in edge_ranges])
  flat_probe_ranges = @lift begin
    [Point2.(Δ.vertices) for Δ in Iterators.flatten($probe_ranges)]
  end

  Δ_plt = poly!(Δ,
                color=(:black,0),
                strokecolor=(:black,1),
                strokewidth=2,
                inspectable=false)
  probe_ranges_plot = poly!(flat_probe_ranges,
                            color=probe_range_colors,
                            colormap=(:rainbow,0.1),
                            strokewidth=0.1,
                            strokecolor=:black,
                            inspectable=false
                           )
  vertex_buttons = scatter!(Δ,
                            overdraw = true,
                            color=button_color,
                            markersize=button_size)


  bcl_line = Observable(Point2[])
  linesegments!(bcl_line,
                color=:black,
                linestyle=[0,5,10],
                linewidth=1.0,
                overdraw=true
               )

  on(events(fig).mousebutton) do event
    if event.button == Mouse.left && event.action == Mouse.press
      plt, idx = pick(fig)
      if plt == vertex_buttons
        k = Keyboard.left_shift in events(fig).keyboardstate ? -1 : 1
        Δ[] = mutate(Δ[],idx, k)
      end
    end
  end

  on(events(fig).mouseposition) do event
    plt, idx = pick(fig)
    if plt == vertex_buttons
      bcl = get_branch_cut_line(Δ[], idx)
      bcl_line[] = Point2.(bcl[1] ∩ Δ[])
    elseif !isempty(bcl_line[])
      bcl_line[] = Point2[]
    end
  end

  on(events(fig).keyboardbutton) do event
    if event.action == Keyboard.press
      if event.key == Keyboard.right
        Δ[] = [1 1;0 1] * Δ[]
      elseif event.key == Keyboard.left
        Δ[] = [1 -1;0 1] * Δ[]
      elseif event.key == Keyboard.up
        Δ[] = [1 0;1 1] * Δ[]
      elseif event.key == Keyboard.down
        Δ[] = [1 0;-1 1] * Δ[]
      end
    end
  end

  return fig
end

CP2 = Polygon([-1,-1],[2,-1],[-1,2]);
CP2_1 = Polygon([-1,0],[0,-1],[2,-1],[-1,2]);
CP2_2 = Polygon([-1,0],[0,-1],[1,-1],[1,0],[-1,2]);
CP2_3 = Polygon([-1,0],[0,-1],[1,-1],[1,0],[0,1],[-1,1]);
Δ1 = Polygon([0,0],[26,0],[26,1],[22,9],[21,10],[20,10]);
