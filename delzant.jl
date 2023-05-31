using LinearAlgebra
using Plots

Base.rationalize(x::Rational) = x
Base.rationalize(x::Integer) = x//1

struct Edge
  λ::Vector{<:Integer}
  c::Rational
  function Edge(λ::Vector{<:Real}, c::T) where T<:Real
    λ = rationalize.(λ); c = rationalize(c)
    d = gcd(λ...);
    new(λ .÷ d, c//d)
  end
end
function Edge(p1::Vector{<:Real}, p2::Vector{<:Real})
  λ = [0 -1; 1 0] * (p2-p1)
  Edge(λ, -λ'*p1)
end
Base.hash(e::Edge, h::UInt) = hash((e.λ,e.c),h)
Base.isequal(e1::Edge,e2::Edge) = (e1.λ == e2.λ && e1.c == e2.c)
Base.:(==)(e1::Edge,e2::Edge) = Base.isequal(e1,e2)

function ev(e::Edge, p::Vector{<:Real})
  return p'*e.λ + e.c
end


function intersect(e1::Edge, e2::Edge)
  e1.λ==e2.λ || e1.λ==-e2.λ ? nothing : -[e1.λ'//1; e2.λ'//1]^-1*[e1.c; e2.c]
end
e1::Edge ∩ e2::Edge = intersect(e1,e2)


struct Polygon
  edges::Vector{<:Edge}
  vertices::Vector{Vector{<:Rational}}
end
function Polygon(edges::Edge...)
  edges = unique(edges) # remove duplicate edges
  intersections = []
  # Collect all intersections
  for (n,e1) in enumerate(edges) # iterate over all unordered pairs
    for e2 in Iterators.drop(edges,n)
      i = e1 ∩ e2
      if !isnothing(i)
        # intersections are labled by their edges
        push!(intersections, (i, Set((e1,e2))))
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
    return nothing
  end

  # Walk through intersections to get them sorted:
  current_intersection = pop!(intersections)
  vertices = [current_intersection[1]]
  e = first(current_intersection[2])
  edges = [e]
  while !isempty(intersections)
    tobedeleted = []
    for (i,intersection) in enumerate(intersections)
      if e in intersection[2]
        current_intersection=intersection
        push!(vertices,current_intersection[1])
        e = only(setdiff(current_intersection[2], Set((e,))))
        push!(edges, e)
        push!(tobedeleted, i)
      end
    end
    if isempty(tobedeleted)
      error("Something went wrong.")
    end
    deleteat!(intersections, tobedeleted)
  end

  # Test if they are in clockwise order, and reverse the order
  if ev(edges[1],vertices[3]) < 0
    reverse!(vertices)
    push!(edges, popfirst!(edge)) |> reverse!
  end

  Polygon(edges, vertices)
end
function Polygon(vertices::Vector{<:Real}...)
  vertices = [rationalize.(v) for v in vertices]
  edges = Edge[]
  for i in 1:length(vertices)
    push!(edges, Edge(vertices[i], vertices[i%end+1]))
  end
  Polygon(edges, vertices)
end

function Base.:*(M::Matrix{<:Real}, Δ::Polygon)
  vertices = [M * v for v in Δ.vertices]
  Polygon(vertices...)
end
function Base.:+(u::Vector{<:Real}, Δ::Polygon)
  vertices = [u + v for v in Δ.vertices]
  Polygon(vertices...)
end

function draw(Δ::Polygon; fill=:red, falpha=0.0, border=:black, balpha=1.0, lw=2.0)
  s = Shape(first.(Δ.vertices), last.(Δ.vertices))
  plot!(s,
        aspect_ratio=:equal,
        legend=false,
        fillcolor=plot_color(fill,falpha),
        linecolor=plot_color(border,balpha),
        linewidth=lw,
       )
end

function draw_probe_ranges(Δ::Polygon; color=:rainbow, alpha=0.15)
  grad = cgrad(color)
  n = length(Δ.edges)
  for i in 1:n
    for range in get_probe_ranges(Δ,i)
      draw(range,
           fill=grad[(i-1)/(n-1)],
           falpha=alpha,
           lw=0.2,
          )
    end
  end
  draw(Δ)
end


function get_default_probe(e::Edge)
  collect(gcdx(e.λ...)[2:3])
end

function get_probe_intervall(e::Vector{Edge}, probe::Vector{<:Integer})
  v = [0 -1; 1 0] * e[2].λ
  a = - (probe' * e[1].λ) / (v' * e[1].λ)
  b = - (probe' * e[3].λ) / (v' * e[3].λ)
  u = sort([a,b])
  return floor(Int,u[1]):ceil(Int,u[2])
end
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

function get_probe_range(Δ::Polygon, edge_index::Integer, probe::Vector{<:Integer})
  l = Δ.vertices[edge_index]
  r = Δ.vertices[edge_index%end + 1]
  e = Δ.edges[edge_index]
  λ = [0 -1;1 0]*probe
  range = Polygon(Edge(λ,-λ'*r), Edge(-λ, λ'*l), Δ.edges...)
  M = [[0 -1;1 0]*e.λ//1 probe//1]
  offset = - e.λ * e.c//sum(t->t^2, e.λ)
  offset + (M*[1 0; 0 1//2]*M^(-1)) * (-offset + range)
end

function get_probe_ranges(Δ::Polygon, edge_index::Integer)
  e = Δ.edges[edge_index]
  ed = [0 -1; 1 0] * e.λ
  dp = get_default_probe(e)
  I = get_probe_intervall(Δ, edge_index, dp)
  ranges = []
  for i in I
    probe = dp + i*ed
    push!(ranges, get_probe_range(Δ, edge_index, probe))
  end
  ranges
end

CP2 = Polygon([-1,-1], [1,-1], [-1,1])
Δ1 = Polygon([0,0],[26,0],[26,1],[22,9],[21,10],[20,10])
