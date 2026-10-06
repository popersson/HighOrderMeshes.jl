###
### Convection in 2D using upwind DG
###
### ∇·(v u) = 0 on the unit square with v = (1, 2x). u = 1 enters through the
### bottom boundary and u = 0 through the left one; right and top are outflow
### boundaries. The two states meet along the streamline y = x², where the
### DG solution oscillates as expected for a discontinuity.

using HighOrderMeshes
using GLMakie

m = mshsquare(16, p=3)                             # boundaries: 1 left, 2 right, 3 bottom, 4 top
d = DGData(m)
vel(x) = (1.0, 2x[1])
bc = Dict(1 => x -> 0.0, 3 => x -> 1.0)
A, b = dg_convection(d, vel, bc)
U = reshape(A \ vec(b), size(m.el))
plot(m, U, mesh_edges=true)
