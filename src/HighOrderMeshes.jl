"""
    HighOrderMeshes

A Julia package for high-order unstructured mesh generation, manipulation, and
finite element assembly on simplex (triangles, tetrahedra) and block (quads,
hexahedra) elements at arbitrary polynomial order.

**Submodule layout** (all exported into the top-level namespace):

| Directory      | Contents                                              |
|:-------------- |:----------------------------------------------------- |
| `mesh/`        | Element geometry, mesh struct, refinement             |
| `meshgen/`     | Mesh generators: basic, sample and airfoil meshes     |
| `basis/`       | Legendre polynomials, quadrature, reference elements  |
| `fem/`         | Precomputed data, CG assembly, DG operators           |
| `viz/`         | Backend-agnostic mesh and solution visualization data |
| `io/`          | Binary `.hom` format, Gmsh import, VTK export         |

Visualization is loaded as a package extension: `using Makie` (or a Makie
backend such as GLMakie or CairoMakie) enables `plot(m)` and `plot(m, u)`.

**Mesh generation cheat sheet** (`p` is the degree, `eg` the element geometry):

```julia
msh = mshsquare(8)                              # 8x8 quads on [0,1]^2
msh = mshsquare(8; eg=Simplex{2}(), p=3)        # 128 triangles, degree 3
msh = mshcube(4; eg=Simplex{3}())               # 384 tetrahedra on [0,1]^3
msh = mshcircle(4; eg=Simplex{2}(), p=3)        # curved triangles on the unit disk
msh = gmsh_sample(:sphere; h=0.3, p=2)          # gmsh: :square, :circle, :cube, :sphere
msh = mshairfoil(:naca0012; aoa=5)              # gmsh: quads around an airfoil

using DistMesh                                  # straight mesh from DistMesh, then curved
fd(p) = sqrt(sum(p.^2)) - 1
msh = HighOrderMesh(distmesh2d(fd, huniform, 0.2, ((-1,-1), (1,1))))
msh = curve_boundary(set_degree(msh, 4), fd)
```

See the help of each function for more examples.
"""
module HighOrderMeshes

using LinearAlgebra, SparseArrays, StaticArrays

# mesh/
export ElementGeometry, Simplex, Block, dim, nvertices, nfaces, nedges, nnodes, vertices
export facemap, edgemap, subgeom, symmetries
export HighOrderMesh, dg_nodes, elgeom, porder, nel
export Neighbor, isboundary, bndtag
export set_ref_nodes, set_degree, set_lobatto_nodes, mkface2nodes
export refine, refine_with_parents, uniref, bndlayer_refine, bndlayer_refine_with_elements
export boundary_nodes, boundary_distance, set_bnd_numbers!, set_bnd_periodic!, unique_mesh_nodes

# meshgen/
export mshhypercube, mshcube, mshsquare, mshline, mshcircle
export ex1mesh, ex1solution, gmsh_sample
export project_points, interp_arclength, curve_boundary
export mshairfoil, airfoil_geo, airfoil_coordinates, naca4

# basis/
export jacobi, djacobi, legendre, dlegendre, legendre01, dlegendre01
export polybasis, dpolybasis, multiindices, equispaced_nodes, tensor_nodes
export gauss_legendre_nodes, gauss_legendre01_nodes, gauss_legendre_quadrature, gauss_legendre01_quadrature
export gauss_lobatto_nodes, gauss_lobatto01_nodes, gauss_lobatto_quadrature, gauss_lobatto01_quadrature
export gauss_radau_nodes, gauss_radau01_nodes, gauss_radau_quadrature, gauss_radau01_quadrature
export quadrature
export FiniteElement, ref_nodes, subelement, corner_nodes, check_conforming
export shapefcns, dshapefcns, interpolate

# fem/
export mkldgswitch, align_with_ldgswitch!
export RefOps, Metric, FaceMetric, FEMData, CGData, DGData, npoints, jacobians, face_refops
export elmat_laplace!, elmat_mass!, elmats, elvec_source, laplace_residual!
export assemble_matrix, assemble_vector, strong_dirichlet
export dg_laplace, dg_convection

# viz/
export viz_mesh, viz_solution, mesh_function_type

# io/
export savemesh, loadmesh
export mshto3dg, mshfrom3dg, gmsh2msh, rungmsh2msh, gmshstr2msh, gmsh_physical_names, vtkwrite

# mesh/ (element topology, no basis dependency)
include("mesh/element_geometry.jl")

# basis/ (polynomials, quadrature, reference element — depends on element_geometry)
include("basis/poly_tools.jl")
include("basis/finite_element.jl")

# mesh/ (mesh struct and utilities — depends on basis)
include("mesh/high_order_mesh.jl")
include("mesh/mesh_utils.jl")
include("mesh/refinement.jl")

# meshgen/
include("meshgen/basic_meshes.jl")
include("meshgen/sample_meshes.jl")
include("meshgen/curved_boundary.jl")
include("meshgen/airfoil.jl")

# fem/
include("fem/dg_utils.jl")
include("fem/fem_data.jl")
include("fem/assembly.jl")
include("fem/dg_operators.jl")

# viz/
include("viz/post_processing.jl")

# io/
include("io/converters.jl")
include("io/io.jl")

end
