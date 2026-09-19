# Quick start

This walks through building a mesh, raising its polynomial degree, plotting
it, and saving it to disk.

```julia
using HighOrderMeshes

msh = mshsquare(8)          # 8x8 quad mesh on [0,1]^2, degree 1
msh = set_degree(msh, 3)    # raise to degree 3

using CairoMakie            # or GLMakie for an interactive window
plot(msh)

savemesh("square.hom", msh)
msh2 = loadmesh("square.hom")
```
