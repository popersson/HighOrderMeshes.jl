###
### Poisson in 2D using the symmetric interior penalty DG method
###
### The problem of ex_cgpoisson.jl, -∇²u = x² on the unit circle with u = 0 on
### the boundary, with the boundary condition imposed weakly.

using HighOrderMeshes
using GLMakie

m = mshcircle(2, p=3)
d = DGData(m)
A, b = dg_laplace(d, Dict(1 => xy -> 0.0))         # boundary 1 (the circle): u = 0
F = elvec_source(d, xy -> xy[1]^2) + b
U = reshape(A \ vec(F), size(m.el))                # DG solution, ns × nel
plot(m, U, mesh_edges=true)
