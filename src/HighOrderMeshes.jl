module HighOrderMeshes

export ElementGeometry, Simplex, Block, dim, name, nvertices, nfaces, nedges, facemap, edgemap, subgeom
export legendre_poly, multivar_legendre_poly, multivar_monomial_poly, eval_poly
export equispaced, ref_nodes
export FiniteElement, elgeom, dim, porder, nbr_ho_nodes, corner_nodes, name
export eval_shapefcns, eval_fcn
export HighOrderMesh, dg_nodes, elgeom, dim, porder, el2nbor, change_ref_nodes, change_degree
export write_matrix, read_matrix!, savemesh, loadmesh, uniref, gmsh2msh, rungmsh2msh, vtkwrite, mkface2nodes, mkldgswitch, boundary_nodes, align_with_ldgswitch!, set_bnd_numbers!, unique_mesh_nodes
export quadrature
export blockmesh_hypercube, mshhypercube, mshcube, mshsquare, mshline, mshcircle
export viz_mesh, viz_solution, mesh_function_type

include("misc_utils.jl")
include("element_geometry.jl")
include("poly_tools.jl")
include("finite_element.jl")
include("high_order_mesh.jl")
include("mesh_utils.jl")
include("basic_meshes.jl")
include("post_processing.jl")
include("simplex_quadrature_autogen.jl")
include("quadrature.jl")
include("../examples/sample_meshes.jl")

end
