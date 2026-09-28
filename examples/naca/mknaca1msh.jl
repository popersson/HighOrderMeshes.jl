###
### NACA 0012 quad mesh for wall-resolved LES
###
### Isotropic elements around the airfoil and in the near wake, and three
### boundary layers on the airfoil that end right behind the trailing edge.
### Julia version of 3DG's mknaca1msh.m. Needs gmsh on the PATH.
###

using HighOrderMeshes

msh = mshairfoil(:naca0012;
                 nfoil=56, hle=0.005, hte=0.005, tband=0.01,   # airfoil surface
                 hnear=0.0175, dnear=0.15,                     # isotropic zone around it
                 gwake=0.03, awake=1, Lwake=4, dwake=0.2,      # isotropic near wake
                 growth=0.3, hmax=15,
                 nbndlayers=3, layerlength=0)
println(msh)

# using GLMakie
# plot(msh)
