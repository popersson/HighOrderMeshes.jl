module HighOrderMeshes

export ElementGeometry, Simplex, Block, dim, name, nvertices, nfaces, nedges, facemap, edgemap, subgeom

export legendre_poly, legendre01_poly, multivar_legendre01_poly, multivar_monomial_poly, eval_poly
export gauss_legendre_nodes, gauss_legendre01_nodes, gauss_legendre_quadrature, gauss_legendre01_quadrature
export gauss_lobatto_nodes, gauss_lobatto01_nodes, gauss_lobatto_quadrature, gauss_lobatto01_quadrature
export equispaced, ref_nodes, quadrature

export FiniteElement, elgeom, dim, porder, nbr_ho_nodes, corner_nodes, name
export eval_shapefcns, eval_fcn

export HighOrderMesh, dg_nodes, elgeom, dim, porder, el2nb, set_ref_nodes, set_degree, set_lobatto_nodes

export uniref, mkface2nodes, mkldgswitch, boundary_nodes, align_with_ldgswitch!, set_bnd_numbers!, set_bnd_periodic!, unique_mesh_nodes

export blockmesh_hypercube, mshhypercube, mshcube, mshsquare, mshline, mshcircle
export viz_mesh, viz_solution, mesh_function_type
export savemesh, loadmesh, savemeshtxt, loadmeshtxt
export mshto3dg, gmsh2msh, rungmsh2msh, gmshstr2msh, vtkwrite

include("element_geometry.jl")
include("simplex_quadrature.jl")
include("poly_tools.jl")
include("finite_element.jl")
include("high_order_mesh.jl")
include("mesh_utils.jl")
include("basic_meshes.jl")
include("post_processing.jl")
include("converters.jl")
include("io.jl")
include("../examples/sample_meshes.jl")

end
