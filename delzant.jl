using LinearAlgebra
using PrettyPrint
using Luxor

Base.rationalize(x::Rational) = x
Base.rationalize(x::Integer) = x//1

struct Edge
  λ::Vector{<:Integer}
  c::Rational
  function Edge(λ::Vector{<:Number}, c::T) where T<:Number
    λ = rationalize.(λ); c = rationalize(c)
    d = gcd(λ...);
    new(λ .÷ d, c//d)
  end
end
function Edge(p1::Vector{<:Number}, p2::Vector{<:Number})
  λ = [0 -1; 1 0] * (p2-p1)
  Edge(λ, -λ'*p1)
end
Base.isequal(e1::Edge,e2::Edge) = (e1.λ == e2.λ && e1.c == e2.c)
Base.:(==)(e1::Edge,e2::Edge) = Base.isequal(e1,e2)

function ev(e::Edge, p::Vector{<:Number})
  return p'*e.λ + e.c
end


function intersect(e1::Edge, e2::Edge)
  e1.λ==e2.λ ? nothing : -[e1.λ'//1; e2.λ'//1]^-1*[e1.c; e2.c]
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
        e = only(setdiff(current_intersection[2], e))
        push!(edges, e)
        push!(tobedeleted, i)
      end
    end
    if isempty(tobedeleted) # The Polygon seems to have redundant edges.
      error("AHHHHHHHHHHH")
    end
    deleteat!(intersections, tobedeleted)
  end
  Polygon(edges, vertices)
end
function Polygon(vertices::Vector{<:Number}...)
  edges = Edge[]
  for i in 1:length(vertices)
    push!(edges,Edge(vertices[i],vertices[i%end+1]))
  end
  Polygon(edges, collect(vertices))
end

function draw(Δ::Polygon)
  @png begin
    pts = 50.0*[Point(p[1],-p[2]) for p in Δ.vertices]
    poly(pts, :stroke, close=true)
  end
end


function get_default_probe(e::Edge)
  collect(gcdx(e.λ...)[2:3])
end

function get_probe_range(e::Vector{Edge}, probe::Vector{<:Integer})
  v = [0 -1; 1 0] * e[2].λ
  a = - (probe' * e[1].λ) / (v' * e[1].λ)
  b = - (probe' * e[3].λ) / (v' * e[3].λ)
  u = sort([a,b])
  return floor(Int,u[1]):ceil(Int,u[2])
end
