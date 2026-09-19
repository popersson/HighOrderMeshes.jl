# HighOrderMeshes.jl

HighOrderMeshes.jl is a light-weight Julia library for high-order
unstructured meshes on simplex and block elements, with a single element
type and a single polynomial degree per mesh.

## Installation

```julia
import Pkg
Pkg.add(url="https://github.com/popersson/HighOrderMeshes.jl")
```

## Hello world

```julia
using HighOrderMeshes
msh = ex1mesh(eg=Simplex{2}())
msh = set_degree(msh, 3)
```

See the [Quick start](@ref) page for a slightly longer walkthrough.
