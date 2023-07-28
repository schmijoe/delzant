using LinearAlgebra
using StaticArrays, CircularArrays
using GLMakie

Base.rationalize(x::Rational) = x
Base.rationalize(x::Integer) = x//1

struct Edge
  λ::SVector{2,<:Integer}
  c::Rational
  """
    `Edge(λ::Vector{<:Real}, c::<:Real)`

  Construct edge/halfspace satisfying equation `λ'*x + c >= 0`.
  """
  function Edge(λ::AbstractVector{<:Real}, c::T) where T<:Real
    λ = rationalize.(λ); c = rationalize(c)
    if λ[1] == λ[2] == 0
      @error "Normal Vector must be non-zero" λ
    end
    d = gcd(λ);
    new(SVector{2}(λ .÷ d), c//d)
  end
end
"""
  `Edge(p1::Vector{<:Real}, p2::Vector{<:Real})`

Consturct edge/hafspace containing points `p1` and `p2`.
"""
function Edge(p1::AbstractVector{<:Real}, p2::AbstractVector{<:Real})
  λ = [0 -1; 1 0] * (p2-p1)
  Edge(λ, -λ'*p1)
end
Base.hash(e::Edge, h::UInt) = hash((e.λ,e.c),h)
Base.isequal(e1::Edge,e2::Edge) = (e1.λ == e2.λ && e1.c == e2.c)
Base.:(==)(e1::Edge,e2::Edge) = Base.isequal(e1,e2)
Base.:(-)(e::Edge) = Edge(-e.λ,-e.c)
Base.:*(M::AbstractMatrix{<:Rational}, e::Edge) = Edge(M'^-1 * e.λ, e.c)
Base.:+(u::AbstractVector{<:Rational}, e::Edge) = Edge(e.λ, e.c - e.λ · u)
Base.:+(o::Rational, e::Edge) = Edge(e.λ, e.c + o)
Base.show(io::IO, e::Edge) = print(io, "⟨$(e.λ),⋅⟩ + $(e.c) ≥ 0")

(e::Edge)(p::AbstractVector{<:Real}) = p'*e.λ + e.c

Base.big(e::Edge) = Edge(big.(e.λ), big(e.c))


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
  edges::CircularVector{Edge}
  vertices::CircularVector{<:Union{SVector{2,<:Rational}, Nothing}}
end
#Polygon(edges::AbstractVector{Edge}, vertices::AbstractVector{<:Union{SVector{2,<:Rational}, Nothing}}) = Polygon(CircularVector(edges), CircularVector(vertices))

Base.big(Δ::Polygon) = Polygon(big.(Δ.edges), map(v->big.(v), Δ.vertices))

function drop_colinear!(vertices::AbstractVector{<:SVector{2,<:Rational}})
  vertices = CircularVector(vertices)
  tobedeleted = []
  for i in eachindex(vertices)
    u = vertices[i + 1] - vertices[i]
    v = vertices[i - 1] - vertices[i]
    if u[1]*v[2] - u[2]*v[1] == 0 && (u != [0,0])
      push!(tobedeleted, i)
    end
  end
  deleteat!(vertices, tobedeleted)
end
"""
Construct polygon given by `vertices`. They must be in counter-clockwise order.
"""
function Polygon(vertices::SVector{2,<:Real}...)
  vertices = drop_colinear!(collect(vertices))
  edges = Edge[]
  for i in eachindex(vertices)
    push!(edges, Edge(vertices[i], vertices[i + 1]))
  end
  Polygon(CircularVector(edges), vertices)
end
Polygon(vertices::AbstractVector{<:Real}...) = Polygon((SVector{2,Rational}(rationalize.(v)) for v in vertices)...)

function Polygon(edges::Edge...)
  edges = sort!(collect(edges),by=e->(atan(e.λ[2],e.λ[1]), e.c))
  edges = unique!(e->e.λ, edges) |> CircularVector
  while !isempty(edges)
    vertices = CircularVector(Vector{Union{Nothing, SVector{2,Rational}}}(nothing, length(edges)))
    vertices[0] = first(edges[-1] ∩ edges[0])
    tobedeleted = Int[]
    for i in eachindex(edges)
      vertices[i] = first(edges[i-1] ∩ edges[i])
      n1 = vertices[i-1] === nothing
      n2 = vertices[i] === nothing
      if (n1 && n2) || (!(n1 || n2) && det([(vertices[i]-vertices[i-1])' ; edges[i-1].λ']) < 0)
        push!(tobedeleted, i-1)
      end
    end
    if isempty(tobedeleted)
      return Polygon(edges, vertices)
    end
    deleteat!(edges, tobedeleted)
  end
  return Polygon(CircularVector(Edge[]), CircularVector(SVector{2,Rational}[]))
end

function Base.:*(M::AbstractMatrix{<:Rational}, Δ::Polygon)
  vertices = map(v->SVector{2}(M*v), Δ.vertices)
  edges = map(e->M*e, Δ.edges)
  if det(M) < 0
    reverse!(vertices)
    reverse!(edges)
  end
  Polygon(edges, vertices)
end
Base.:+(u::AbstractVector{<:Rational}, Δ::Polygon) = Polygon(map(e->u+e, Δ.edges),map(v->u + v, Δ.vertices))
Base.:+(o::Rational, Δ::Polygon) = Polygon(map(e->o+e, Δ.edges)...)
Base.isempty(Δ::Polygon) = isempty(Δ.edges)

(Δ::Polygon)(p::AbstractVector{<:Real}) = minimum(e->e(p), Δ.edges)


get_affine_edge_length(Δ::Polygon, i::Int) = gcd(Δ.vertices[i+1] - Δ.vertices[i])

function Base.show(io::IO, Δ::Polygon)
  for i in eachindex(Δ.vertices)
    w = Rational[0 -1;1 0] * (Δ.edges[i-1].λ-Δ.edges[i].λ)
    k = gcd(w)
    area = det([Δ.edges[i-1].λ' ; Δ.edges[i].λ'])
    println(io::IO, "$(sprint(show, Δ.vertices[i])); [q,p] = $(w .÷ k); k = $k; area = $area")
  end
  for i in eachindex(Δ.edges)
    println(io::IO, "$(sprint(show,Δ.edges[i])); Length: $(get_affine_edge_length(Δ,i))")
  end
end

"""
  `intersect(e::Edge, Δ::Polygon)`

Intersections of `e` with of edges of `Δ`.
"""
function intersect(e::Edge, Δ::Polygon)
  intersections = []
  for (i,e1) in enumerate(Δ.edges)
    int = e1 ∩ e
    if !isnothing(int[1])
      dir = Δ.vertices[i+1] - Δ.vertices[i]
      t = dir ⋅ (int[1] - Δ.vertices[i])
      if 0 <= t <= dir⋅dir
        push!(intersections,int)
      end
    end
  end
  unique!(i->i[1], intersections)
  sort!(intersections, by=i->i[2])
  first.(intersections)
end

"""
  `refine!(Δ::Polygon, e::Edge)`

Get vector of vertices of `Δ` with added vertices for intersections with `e`
"""
function refine!(Δ::Polygon, halfspaces::Edge...)
  intersections = []
  for e in halfspaces
    for (i,e1) in enumerate(Δ.edges)
      int = e1 ∩ e
      if !isnothing(int[1])
        dir = Δ.vertices[i+1] - Δ.vertices[i]
        t = dir ⋅ (int[1] - Δ.vertices[i])
        if 0 < t < dir⋅dir
          push!(intersections,(int..., i, t))
        end
      end
    end
  end
  sort!(intersections, by=i->(i[3],i[4]), rev=true)
  for int in intersections
    insert!(Δ.vertices.data, int[3] + 1, int[1])
    insert!(Δ.edges.data, int[3] + 1, Δ.edges[int[3]])
  end
  Δ
end

"""
  `slice(Δ::Polygon, halfspaces::Edge...)`

Get Polygon of intersection of `Δ` with `halfspaces`
"""
function slice(Δ::Polygon, halfspaces::Edge...)
  vertices = refine!(deepcopy(Δ), halfspaces...).vertices
  tobedeleted = Int[]
  for (i,v) in enumerate(vertices)
    for e in halfspaces
      if e(v) < 0
        push!(tobedeleted, i)
        break
      end
    end
  end
  Polygon(deleteat!(vertices, tobedeleted)...)
end

"""
Get some possible probe direction from `e`
"""
get_default_probe(e::Edge) = collect(gcdx(e.λ...)[2:3])

"""
  `get_probe_intervall(e::AbstractVector{Edge}, probe::Vector{<:Integer})`

Get range of `k` for which `probe + k*[0 -1; 1 0]*e[2].λ` is a reasonable probe in the tripple of edges `e[1]`, `e[2]` ,`e[3]`.

"""
function get_probe_intervall(e::AbstractVector{Edge}, probe::Vector{<:Integer})
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
get_probe_intervall(Δ::Polygon, i::Integer, probe::Vector{<:Integer}) = get_probe_intervall(Δ.edges[i-1:i+1], probe)

"""
  `get_probe_range(Δ::Polygon, edge_index::Integer, probe::Vector{<:Integer})`

Get a polygon describing which fibres in Δ can be displaced by a probe shot from edge `Δ.edges[edge_index]` in direction `probe`.
"""
function get_probe_range(Δ::Polygon, edge_index::Integer, probe::Vector{<:Integer})::Polygon
  l = Δ.vertices[edge_index]
  r = Δ.vertices[edge_index%end + 1]
  e = Δ.edges[edge_index]
  λ = [0 -1;1 0]*probe

  range = slice(Δ, Edge(λ,-λ'*r), Edge(-λ, λ'*l))
  M = [[0 -1;1 0]*e.λ//1 probe//1]
  l + (M * [1 0; 0 1//2] * M^(-1)) * (-l + range)
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
  λ = e2.λ - e1.λ; kdist = gcd(λ); λ = λ .÷ kdist;
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
  get_branch_cut_line(Δ.edges[vertex_index-1], Δ.edges[vertex_index])

"""
  `mutate!(Δ::Polygon, branch_cut_line::Tuple{Edge, <:Integer})`

Get a polygon obtained by appling a shear matrix `M` corresponding to `branch_cut_line = (e, k)` to the positive side of `Δ` wrt. `e`, where `k` determines the power of `M`. If `k` is negative it will instead be applied to the negative side of `Δ` wrt. `e`.
"""
function mutate!(Δ::Polygon, branch_cut_line::Tuple{Edge, <:Integer})
  # 1. Insert missing vertices at intersections of branch_cut_line with polygon:
  refine!(Δ, branch_cut_line[1])
  offset = - branch_cut_line[1].λ * branch_cut_line[1].c//sum(t->t^2, branch_cut_line[1].λ)

  # 2. Apply mutation
  p,q = [0 1;-1 0]*branch_cut_line[1].λ; k=branch_cut_line[2]
  M = Rational[1+k*p*q -k*p^2; k*q^2 1-k*p*q]
  l1 = l2 = branch_cut_line[1](Δ.vertices[1])
  for i in eachindex(Δ.vertices)
    l1 = l2
    l2 = branch_cut_line[1](Δ.vertices[i+1])
    if l1*k > 0
      Δ.vertices[i] = offset + M*(Δ.vertices[i]-offset)
    end
    if l1*k > 0 || l2*k > 0
      Δ.edges[i] = offset + M*(-offset + Δ.edges[i])
    end
  end

  # Remove duplicate edges
  for i in reverse(eachindex(Δ.edges))
    if Δ.edges[i] == Δ.edges[i-1]
      deleteat!(Δ.edges, i)
      deleteat!(Δ.vertices, i)
    end
  end
  Δ
end

"""
  `mutate!(Δ::Polygon, vertex_index::Integer, k::Integer=1)`

Assume the corner `Δ.vertices[vertex_index]` is almost toric, and try to perform `k` mutations.
"""
function mutate!(Δ::Polygon, vertex_index::Integer, k::Integer=1)
  bcl = get_branch_cut_line(Δ, vertex_index)
  if abs(k) > bcl[2]
    @warn "k is too large." k bcl
  end
  k = (sign(atan(bcl[1].λ[2],bcl[1].λ[1])) == 1) ? k : -k
  mutate!(Δ, (bcl[1],k))
end

function mutate_with!(range::Polygon, Δ::Polygon, vertex_index, k::Integer=1)
  bcl = get_branch_cut_line(Δ, vertex_index)
  if abs(k) > bcl[2]
    @warn "k is too large." k bcl
  end
  k = (sign(atan(bcl[1].λ[2],bcl[1].λ[1])) == 1) ? k : -k
  mutate!(range, (bcl[1],k))
end

# Enable Makie to plot Polygon
Makie.convert_arguments(P::PointBased, Δ::Polygon) = (decompose(Point2f, Point2f.(Δ.vertices)),)

function draw(Δ::Polygon;
  fig=nothing,
  ax=nothing, 
  outline_width = 4,
  outline_color = (:black, 1),
  draw_contours = true,
  contour_width = 0.2,
  contour_density::Rational = 1//5,
  contour_color = (:black, 1)
  )
  if fig === nothing
    fig = Figure()
  end
  if ax===nothing
    ax = Axis(fig[1,1], autolimitaspect = 1.0)
  end
  current_axis!(ax)

  if draw_contours
    for inset in 0:contour_density:10
      cont = (-inset + Δ)
      if isempty(cont)
        break
      end
      poly!(Point2f.(cont.vertices.data),
            color=(:black,0),
            strokecolor=contour_color,
            strokewidth=contour_width)
    end
  end

  poly!(Point2f.(Δ.vertices.data),
        color=(:black, 0),
        strokecolor=outline_color,
        strokewidth=outline_width)

  return fig,ax
end

function interact(Δ::Polygon; button_size = 30, button_color = (:black, 0.1), show_probe_ranges = true)
  fig = Figure()
  ax = Axis(fig[1,1], autolimitaspect = 1.0)

  Δ = Observable(Δ)
  on(Δ) do Δ
    println(Δ)
  end
  if show_probe_ranges
    probe_ranges = Observable(get_probe_ranges(Δ[]))
    probe_range_colors = @lift(Integer[i for (i,edge_ranges) in enumerate($probe_ranges) for _ in edge_ranges])
    flat_probe_ranges = @lift begin
      [Point.(Δ.vertices.data) for Δ in Iterators.flatten($probe_ranges)]
    end

    probe_ranges_plot = poly!(flat_probe_ranges,
                              color=probe_range_colors,
                              colormap=(:rainbow,0.1),
                              strokewidth=0.1,
                              strokecolor=:black,
                              inspectable=false
                             )
  end
  Δ_plt = poly!(Δ,
                color=(:black,0),
                strokecolor=(:black,1),
                strokewidth=4,
#                inspector_label = (plt, idx, pos) -> "Affine Lenght $(get_affine_edge_length(Δ[], i))"
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
        if show_probe_ranges
          for edge_range in probe_ranges[]
            for i in eachindex(edge_range)
              edge_range[i] = mutate_with!(edge_range[i], Δ[], idx, k)
            end
          end
        end
        Δ[] = mutate!(Δ[], idx, k)
        show_probe_ranges && notify(probe_ranges)
      end
    end
  end

  on(events(fig).mouseposition) do event
    plt, idx = pick(fig)
    if plt == vertex_buttons
      bcl = get_branch_cut_line(Δ[], idx)
      bcl_line[] = Point.(intersect(bcl[1], Δ[]))
    elseif !isempty(bcl_line[])
      bcl_line[] = Point2[]
    end
  end

  on(events(fig).keyboardbutton) do event
    if event.action == Keyboard.press
      if event.key == Keyboard.right
        M = [1 1;0 1//1]
        if show_probe_ranges
          for edge_range in probe_ranges[]
            for i in eachindex(edge_range)
              edge_range[i] = M*edge_range[i]
            end
          end; notify(probe_ranges)
        end
        Δ[] = M * Δ[]
      elseif event.key == Keyboard.left
        M = [1 -1;0 1//1]
        if show_probe_ranges
          for edge_range in probe_ranges[]
            for i in eachindex(edge_range)
              edge_range[i] = M*edge_range[i]
            end
          end; notify(probe_ranges)
        end
        Δ[] = M * Δ[]
      elseif event.key == Keyboard.up
        M = [1 0;1 1//1]
        if show_probe_ranges
          for edge_range in probe_ranges[]
            for i in eachindex(edge_range)
              edge_range[i] = M*edge_range[i]
            end
          end; notify(probe_ranges)
        end
        Δ[] = M * Δ[]
      elseif event.key == Keyboard.down
        M = [1 0;-1 1//1]
        if show_probe_ranges
          for edge_range in probe_ranges[]
            for i in eachindex(edge_range)
              edge_range[i] = M*edge_range[i]
            end
          end; notify(probe_ranges)
        end
        Δ[] = M * Δ[]
      elseif event.key == Keyboard.p && show_probe_ranges
        probe_ranges[] = get_probe_ranges(Δ[])
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
Δ2 = Polygon([0,0],[11,0],[11,1],[10,3],[9,4],[8,4]);
Δ3 = Polygon([-3,0],[-2,-2],[-1,-3],[0,-3],[2,-2],[3,-1],[3,0],[2,2],[1,3],[0,3],[-2,2],[-3,1]);
Δ4 = Polygon([-4,0],[-2,-4],[-1,-5],[0,-5],[4,-3],[5,-2],[5,-1],[3,3],[2,4],[1,4],[-3,2],[-4,1]);
Δ5 = Polygon([-5,-5],[-3,-7],[-1,-7],[5,-5],[7,-3],[7,-1],[5,5],[3,7],[1,7],[-5,5],[-7,3],[-7,1]);
