###
### Poisson in 2D using CG finite elements
###
### Solve -∇²u = x² on the unit circle with u = 0 on the boundary.

using HighOrderMeshes
using GLMakie

m = mshcircle(2, p=3)
d = CGData(m)                                      # reference operators and metric terms
A = assemble_matrix(m.el, elmats(elmat_laplace!, d))
f = assemble_vector(m.el, elvec_source(d, xy -> xy[1]^2))
A, f = strong_dirichlet(A, f, boundary_nodes(m))
u = A \ f
plot(m, u, mesh_edges=true)
