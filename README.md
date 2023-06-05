# Delzant
Calculations with 2D Delzant and almost toric moment Polytopes

## Requirements

Uses the `GLMakie package.

## Usage

In [Julia](https://julialang.org):

```julia
include("delzant.jl")
P = Polygon([0,0],[1,0],[1,1],[0,1]) # vertices in counter-clockwise order
interact(P)
```
Then you can click a corner to do a mutation there, arrow keys apply a shear matrix, and ctrl+click resizes the viewport to fit the polygon.

There are also some predefined Polygons like `CP2` and `CP2_1`, `CP2_2`, `CP2_3` being its 1,2,3 point blowups.
