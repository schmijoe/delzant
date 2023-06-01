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
#Base.:(<)(e1::Edge,e2::Edge) = (det([e1.λ e2.λ]) < 0) # counter-clockwise ordering at intersection

function ev(e::Edge, p::Vector{<:Real})
  return p'*e.λ + e.c
end


function intersect(e1::Edge, e2::Edge) # Intersection with sign
  A = [e1.λ'//1; e2.λ'//1]; detA = det(A);
  (detA == 0 ? nothing : -A^-1*[e1.c; e2.c], sign(detA))
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
  sort!(intersections, by=i->(-atan(i[2][2].λ...)))
  vertices=first.(intersections)
  edges=[int[2][2] for int in intersections]
  Polygon(edges, vertices)
end
function Polygon(vertices::Vector{<:Real}...; ensure_convex=false)
  vertices = [rationalize.(v) for v in vertices]
  edges = Edge[]
  for i in eachindex(vertices)
    push!(edges, Edge(vertices[i], vertices[i%end+1]))
  end
  ensure_convex ? Polygon(edges...) : Polygon(edges, vertices)
end

Base.:*(M::Matrix{<:Real}, Δ::Polygon) = Polygon([M * v for v in Δ.vertices]...)
Base.:+(u::Vector{<:Real}, Δ::Polygon) = Polygon([u + v for v in Δ.vertices]...)


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
  draw(Δ, lw=0.8)
end


get_default_probe(e::Edge) = collect(gcdx(e.λ...)[2:3])

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
  offset + (M * [1 0; 0 1//2] * M^(-1)) * (-offset + range)
end

function get_probe_ranges(Δ::Polygon, edge_index::Integer)
  e = Δ.edges[edge_index]
  ed = [0 -1; 1 0] * e.λ
  dp = get_default_probe(e)
  ranges = []
  for i in get_probe_intervall(Δ, edge_index, dp)
    probe = dp + i*ed
    push!(ranges, get_probe_range(Δ, edge_index, probe))
  end
  ranges
end

function get_branch_cut_line(Δ::Polygon, vertex_index::Integer)
  v = Δ.vertices[vertex_index]
  e1 = Δ.edges[vertex_index]
  e2 = vertex_index == 1 ? Δ.edges[end] : Δ.edges[vertex_index - 1]
  
  pq = e1.λ - e2.λ; kdist = gcd(pq...); pq = pq .÷ kdist;
  dist = e1.λ[2]*pq[1] - e1.λ[1]*pq[2]
  if kdist % dist != 0
    @error "I don't think this is a potential almost toric corner"
  end
  (Edge(pq, pq'*v), kdist ÷ dist)
end

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
  sort!(intersections)
  uint = unique(i->(i[1]), intersections)
  intersections = filter!(int1->(count(int2->(int2[1] == int1[1]), intersections) == 1),
                          uint) # select intersections that occur exactly once
  sort!(intersections, by=i->(i[2]), rev=true) # reverse order by the edge index
  vertices = copy(Δ.vertices)
  for (int, edge_index) in intersections
    insert!(vertices, edge_index+1, int)
  end

  # 2. Apply mutation
  p,q = [0 -1;1 0]*branch_cut_line[1].λ; k=branch_cut_line[2]
  M = [1+k*p*q -k*p^2; k*q^2 1-k*p*q]
  offset = intersections[1][1]
  for i in eachindex(vertices)
    l = ev(branch_cut_line[1], vertices[i])
    if l*k > 0
      vertices[i] = offset + M*(vertices[i]-offset)
    end
  end
  Polygon(vertices..., ensure_convex=true)
end
mutate(Δ::Polygon, vertex_index::Integer) = mutate(Δ, get_branch_cut_line(Δ, vertex_index))

CP2 = Polygon([-1,-1],[2,-1],[-1,2]);
CP2_1 = Polygon([-1,0],[0,-1],[2,-1],[-1,2]);
CP2_2 = Polygon([-1,0],[0,-1],[1,-1],[1,0],[-1,2]);
CP2_3 = Polygon([-1,0],[0,-1],[1,-1],[1,0],[0,1],[-1,1]);
Δ1 = Polygon([0,0],[26,0],[26,1],[22,9],[21,10],[20,10]);
