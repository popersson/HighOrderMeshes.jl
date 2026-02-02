###
### Poisson in 2D using CG finite elements
###

using HighOrderMeshes
using GLMakie

include(joinpath(@__DIR__, "assemble_utils.jl"))

# Solve and plot -∇²u = x² with zero Dirichlet boundary conditions on the unit circle
m = mshcircle(2,p=3)
pc = FEM_precomp(m)
u,A,f = cg_poisson(m, pc, xy->xy[1]^2)
#plot(m,u, contours=10, mesh_edges=true) # Contours not supported in GLMakie
plot(m,u, mesh_edges=true)
