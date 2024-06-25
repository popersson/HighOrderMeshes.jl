###
### Example of half-closed reference nodes (Radau)
### Change nodes of mesh, evaluate Vandermondes at nodes + at endpoints
###

using HighOrderMeshes
using Plots
using FastGaussQuadrature

nel = 10
porder = 4

m = mshline(nel)

# New reference nodes
ξ,w = gaussradau(porder+1)
m = change_ref_nodes(m, ξ)

# Shape functions and gradient at new reference nodes
V = eval_shapefcns(m.fe, ξ)                   # Should be identity
dV = eval_shapefcns(m.fe, ξ, gradient=true)

# Shape functions and gradient at left/right of element
Vendpoints = eval_shapefcns(m.fe, [0.0;1.0])
dVendpoints = eval_shapefcns(m.fe, [0.0;1.0], gradient=true)

x = dg_nodes(m)
utest(x) = exp(-(x-0.5).^2 / 0.1^2)
u = utest.(x)

plot(m, u)
