###
### Convection-diffusion in 2D using upwind DG with interior penalty diffusion
###
### -ε∇²u + ∇·(v u) = 0 on the unit square with v = (1, 2x). The inflow data of
### ex_dgconvection.jl (u = 1 on the bottom, u = 0 on the left) are Dirichlet
### conditions for both operators; the outflow boundaries on the right and top
### get the natural conditions, homogeneous Neumann for the diffusion.

using HighOrderMeshes
using GLMakie

m = mshsquare(16, p=3)                             # boundaries: 1 left, 2 right, 3 bottom, 4 top
d = DGData(m)
ε = 0.01
vel(x) = (1.0, 2x[1])
bc = Dict(1 => x -> 0.0, 3 => x -> 1.0)
C, bC = dg_convection(d, vel, bc)
K, bK = dg_laplace(d, bc)
U = reshape((C + ε*K) \ vec(bC + ε*bK), size(m.el))
plot(m, U, mesh_edges=true)
