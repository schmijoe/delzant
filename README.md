# Delzant
Calculations with 2D Delzant and Almost Toric Moment Polytopes

## Usage

In [Julia](https://julialang.org):

```julia
include("delzant.jl")
P = Polygon([0,0],[1,0],[1,1],[0,1]) # vertices in counter-clockwise order
interact(P)
```

There are also some predefined Polygons like `CP2` and `CP2_1`, `CP2_2`, `CP2_3` being its 1,2,3 point blowups.
