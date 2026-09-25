###
### NACA 0012 quad mesh with boundary layer refinement
###
### Julia version of 3DG's mknaca1msh.m. The coarse mesh is generated from
### naca.geo, which needs gmsh on the PATH.
###

using HighOrderMeshes

"""
    mknaca1msh(ref=0, nbndlayers=4)

Mesh the NACA 0012 airfoil with p=3 quads using Gmsh, refine it uniformly
`ref` times and then add `nbndlayers` boundary layer refinements along the
airfoil (boundary 1).
"""
function mknaca1msh(ref=0, nbndlayers=4)
    msh = rungmsh2msh(joinpath(@__DIR__, "naca.geo"); porder=3)
    msh = set_lobatto_nodes(msh)
    msh = uniref(msh, ref)
    bndlayer_refine(msh, 1, nbndlayers)
end

msh = mknaca1msh(1, 3)
println(msh)

# using GLMakie
# plot(msh)
