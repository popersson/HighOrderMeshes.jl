###
### Airfoil quad meshes for RANS, with a wake line
###
### Elements stretched along the wake, and boundary layers on the airfoil
### that continue into the wake as in a C-mesh and thin out downstream.
### Needs gmsh on the PATH.
###

using HighOrderMeshes

# NACA 0012 with the wake at 5 degrees and 8 boundary layers
msh = mshairfoil(:naca0012; aoa=5, nbndlayers=8)
println(msh)

# RAE 2822 (AGARD case 9 angle), finer at the leading edge, 10 layers
msh = mshairfoil(:rae2822; aoa=2.79, nfoil=64, hle=0.002, nbndlayers=10)
println(msh)

# Cambered NACA 2412, all layers continued to the far field
msh = mshairfoil(naca4("2412"); aoa=4, nbndlayers=4, layerlength=Inf, layerratio=1)
println(msh)

# using GLMakie
# plot(msh)
