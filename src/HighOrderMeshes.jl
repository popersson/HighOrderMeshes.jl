module HighOrderMeshes

export ElementGeometry, Simplex, Block, dim, name, nvertices, nfaces, nedges, facemap, edgemap, subgeom
export legendre_poly_neg1pos1, legendre_poly, multivar_legendre_poly, multivar_monomial_poly, eval_poly
export equispaced, ref_nodes
export FiniteElement, elgeom, dim, porder, nbr_ho_nodes, corner_nodes, name
export eval_shapefcns, eval_fcn
export HighOrderMesh, dg_nodes, elgeom, dim, porder, el2nbor, change_ref_nodes, change_degree, change_to_lobatto_nodes
export write_matrix, read_matrix!, savemesh, loadmesh, savemeshtxt, loadmeshtxt, uniref, mkface2nodes, mkldgswitch, boundary_nodes, align_with_ldgswitch!, set_bnd_numbers!, set_bnd_periodic!, unique_mesh_nodes
export quadrature, gauss_legendre_quadrature, gauss_lobatto_nodes, gauss_lobatto_quadrature
export FEM_precomp, eval_gϕx, elmat_mass, elmat_poisson, elres_rhs, assemble_matrix, assemble_vector, strong_dirichlet!, cg_mass, cg_poisson
export blockmesh_hypercube, mshhypercube, mshcube, mshsquare, mshline, mshcircle
export viz_mesh, viz_solution, mesh_function_type
export mshto3dg, gmsh2msh, rungmsh2msh, gmshstr2msh, vtkwrite

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
include("assemble_utils.jl")
include("converters.jl")
include("../examples/sample_meshes.jl")

end
