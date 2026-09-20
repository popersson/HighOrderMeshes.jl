# Solver node sets

A `HighOrderMesh` always carries a conforming node set, so that nodes
coincide across element faces and the mesh fits the compact `x`, `el`, `nb`
form. A solver is free to work on other nodes. Discontinuous Galerkin
spectral element methods (DG-SEM), for example, often use Gauss-Legendre or
Gauss-Radau points, which do not contain the element vertices and therefore
cannot be mesh nodes. The recipe is: keep the mesh as it is, build a second
reference element for the solver's nodes, and move data between the two with
[`interpolate`](@ref).

## Building the solver element

The `FiniteElement` constructors validate node sets by default, because mesh
nodes must be conforming. `check=false` skips the validation for a solver
element that is deliberately non-conforming:

```julia
using HighOrderMeshes

m  = set_degree(mshsquare(8), 3)          # mesh with conforming degree-3 nodes
p  = porder(m)
fe_sol = FiniteElement(Block{2}(), gauss_legendre01_nodes(p+1); check=false)
```

Any unisolvent node set of the right size works: a 1D line for blocks, as
above, or an explicit `nnodes × D` matrix for simplices. Half-closed nodes are
available as `gauss_radau01_nodes(p+1)`.

## Geometry and metric terms at the solver nodes

The mesh's shape functions evaluate the geometry map anywhere in the
reference element, so the coordinates and Jacobians at the solver nodes are
two lines:

```julia
ξ = ref_nodes(fe_sol)
x = interpolate(shapefcns(m.fe, ξ),  dg_nodes(m))    # nξ × nel × D      coordinates
J = interpolate(dshapefcns(m.fe, ξ), dg_nodes(m))    # nξ × nel × D × D  J[:,:,i,j] = ∂x_i/∂ξ_j
```

With Gauss-Legendre nodes the quadrature weights come from the same rule:
`gauss_legendre01_quadrature(p+1)` in 1D, or `quadrature(Block{2}(), 2p+1)`,
whose points coincide with `ξ`, for the tensor-product weights.

If the solver degree `q` is lower than the mesh degree, DG-SEM schemes
usually want the geometry represented in the solver's own polynomial space so
that the discrete metric identities hold. On the mesh side that is simply
`set_degree(m, q)` before taking the metric terms.

## Plotting and exporting a solver field

`interpolate(fe_sol, m.fe, u_sol)` re-nodalizes a field from the solver
nodes to the mesh nodes. When both elements have the same degree this is
exact, because the Lagrange interpolants through either node set are the same
polynomial. The result is an ordinary DG field on the mesh nodes, so it plots
and exports like any other:

```julia
u_sol  = x[:, :, 1] .^ 2                       # any field on the solver nodes
u_mesh = interpolate(fe_sol, m.fe, u_sol)      # nnodes(m.fe) × nel
plot(m, u_mesh)                                # with a Makie backend loaded
vtkwrite("solution.vtk", m, u_mesh)
```

Plotting can also evaluate the solver's interpolant directly, which skips the
transfer and also covers a solver degree different from the mesh degree:

```julia
plot(m, u_sol; fe=fe_sol)
viz_solution(m, u_sol; fe=fe_sol)              # backend-independent plot data
```
