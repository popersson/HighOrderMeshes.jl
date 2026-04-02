"""
    HighOrderMeshes

A Julia package for high-order unstructured mesh generation, manipulation, and
finite element assembly on simplex (triangles, tetrahedra) and block (quads,
hexahedra) elements at arbitrary polynomial order.

**Submodule layout** (all exported into the top-level namespace):

| Directory      | Contents                                              |
|:-------------- |:----------------------------------------------------- |
| `mesh/`        | Element geometry, mesh struct, refinement, generators |
| `basis/`       | Legendre polynomials, quadrature, reference elements  |
| `fem/`         | FEM assembly, CG solvers, DG utilities                |
| `viz/`         | Backend-agnostic mesh and solution visualization data |
| `io/`          | Binary `.hom` format, Gmsh import, VTK export         |

Visualization backends are loaded as package extensions:
- `using Makie` (or a Makie backend) enables `plot(m)` and `plot(m, u)` via Makie.
- `using Plots, TriplotRecipes` enables the same via Plots.jl.
"""
module HighOrderMeshes

using LinearAlgebra, SparseArrays, StaticArrays, Serialization

# mesh/
export ElementGeometry, Simplex, Block, dim, name, nvertices, nfaces, nedges, facemap, edgemap, subgeom
export HighOrderMesh, dg_nodes, elgeom, dim, porder, nnodes, nel
export el2nb, set_ref_nodes, set_degree, set_lobatto_nodes, mkface2nodes
export uniref, boundary_nodes, set_bnd_numbers!, set_bnd_periodic!, unique_mesh_nodes
export blockmesh_hypercube, mshhypercube, mshcube, mshsquare, mshline, mshcircle

# basis/
export legendre_poly, legendre01_poly, multivar_legendre01_poly, multivar_monomial_poly, eval_poly
export gauss_legendre_nodes, gauss_legendre01_nodes, gauss_legendre_quadrature, gauss_legendre01_quadrature
export gauss_lobatto_nodes, gauss_lobatto01_nodes, gauss_lobatto_quadrature, gauss_lobatto01_quadrature
export equispaced, ref_nodes, quadrature
export FiniteElement, elgeom, dim, porder, nbr_ho_nodes, corner_nodes, name
export eval_shapefcns, eval_field

# fem/
export mkldgswitch, align_with_ldgswitch!
export FEM_precomp, eval_gϕx
export elmat_mass, elmat_laplace, elres_source
export assemble_matrix, assemble_vector, strong_dirichlet!
export cg_mass, cg_poisson

# viz/
export viz_mesh, viz_solution, mesh_function_type

# io/
export savemesh, loadmesh, savemeshtxt, loadmeshtxt
export mshto3dg, gmsh2msh, rungmsh2msh, gmshstr2msh, vtkwrite

# mesh/ (element topology, no basis dependency)
include("mesh/element_geometry.jl")

# basis/ (polynomials, quadrature, reference element — depends on element_geometry)
include("basis/simplex_quadrature.jl")
include("basis/poly_tools.jl")
include("basis/finite_element.jl")

# mesh/ (mesh struct and utilities — depends on basis)
include("mesh/high_order_mesh.jl")
include("mesh/mesh_utils.jl")
include("mesh/basic_meshes.jl")

# fem/
include("fem/dg_utils.jl")
include("fem/assembly.jl")

# viz/
include("viz/post_processing.jl")

# io/
include("io/converters.jl")
include("io/io.jl")

include("../examples/sample_meshes.jl")

end
