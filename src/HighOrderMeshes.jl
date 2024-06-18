module HighOrderMeshes

export ElementGeometry, Simplex, Block, dim, name, nvertices, nfaces, nedges, facemap, edgemap, linegeom, facegeom
export legendre_poly, multivar_legendre_poly, multivar_monomial_poly, eval_poly
export equispaced, ref_nodes
export FiniteElement, elgeom, dim, porder, nbr_ho_nodes, corner_nodes, name
export eval_shapefcns, eval_fcn
export HighOrderMesh, dg_nodes, elgeom, dim, porder, el2nbor, change_ref_nodes, change_degree
export write_matrix, read_matrix!, savemesh, loadmesh, uniref, gmsh2msh, rungmsh2msh, vtkwrite
export blockmesh_hypercube, mshcube, mshsquare, mshline
export plot

include("element-geometry.jl")
include("finite-element.jl")
include("high-order-mesh.jl")
include("mesh-utils.jl")
include("basic-meshes.jl")

using Plots, TriplotRecipes

include("plotting.jl")

include("../examples/sample-meshes.jl")

end
