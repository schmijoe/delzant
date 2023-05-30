using LinearAlgebra
using Luxor
using Polyhedra

struct Polygon
  edges::Vector{<:HalfSpace}
  vertices::Vector{Vector{<:Rational}}
end

function Polygon(Δ::Polyhedron)
  start_idx = first(eachindex(points(vrep(Δ))))
  v_idx = start_idx
  vertices = [first(points(vrep(Δ)))]
  edges = typeof(first(halfspaces(hrep(Δ))))[]
  h_idx = nothing
  for i in 1:10
    for h_idx1 in incidenthalfspaceindices(Δ,v_idx)
      if h_idx1 != h_idx
        h_idx = h_idx1
        break
      end
    end
    push!(edges, get(Δ,h_idx))
    for v_idx1 in incidentpointindices(Δ,h_idx)
      if v_idx1 != v_idx
        v_idx = v_idx1
        break
      end
    end
    if v_idx == start_idx
      break
    end
    push!(vertices, get(Δ,v_idx))
  end

  Polygon(edges, vertices)
end

function Polygon(points...)
  Polygon(polyhedron(convexhull(points...)))
end

function draw(Δ::Polygon)
  @png begin
    pts = 50.0*[Point(p[1],-p[2]) for p in Δ.vertices]
    poly(pts, :stroke, close=true)
  end
end


function get_default_probe(e::HalfSpace)
  [ gcdx(e.a...)[i] for i in 2:3 ]
end

function get_probe_range(e::Vector{HalfSpace}, probe::Vector{<:Integer})
  v = [0 -1; 1 0] * e[2].a
  a = - (probe' * e[1].a) / (v' * e[1].a)
  b = - (probe' * e[3].a) / (v' * e[3].a)
  u = sort([a,b])
  return floor(Int,u[1]):ceil(Int,u[2])
end

function get_probe_polygons(Δ::Polyhedron)
  for h_idx in eachindex(halfspaces(Δ))
    e0 = get(p, h_idx)
    vertices = incidentpoints(p, h_idx)
    neighbors = [ e
                 for e in incidenthalfspaces(Δ,p_idx)
                 for p_idx in incidentpointindices(Δ,h_idx)
                 if e != e0
                ]

    probe0 = get_default_probe(e0)
    probe_range = get_probe_range([neighbors[1],e0,neighbors[2]],probe0)
    v = [0 -1; 1 0] * e[2].a
  end
end

