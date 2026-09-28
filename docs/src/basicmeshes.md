# Basic meshes

A gallery of the mesh constructors built into the package, each returning a
`HighOrderMesh` ready to use, plot, or refine further.

## Square

```julia
msh = mshsquare(8; p=2)   # 8x8 quads on [0,1]^2, degree 2
```

## Cube

```julia
msh = mshcube(4; p=1)     # 4x4x4 hexahedra on [0,1]^3
```

## Circle

```julia
msh = mshcircle(3; p=3, shape=:full)   # transfinite quad mesh on the unit disk
```

## Curved sample mesh

```julia
msh = ex1mesh(eg=Simplex{2}(), nref=2)   # sinusoidally curved triangle mesh
```

## Gmsh sphere

```julia
msh = gmsh_sphere(hmax=0.3, porder=2)    # requires gmsh on the PATH
```

## Airfoils

`mshairfoil` meshes the region around an airfoil with curved quads (it
requires gmsh on the PATH). A thin structured band of quads follows the
airfoil and a structured wedge follows a wake line from the trailing edge to
the far field; the rest is unstructured. Boundary layers are added by
isoparametric refinement and continue into the wake as in a C-mesh.

```julia
msh = mshairfoil(:naca0012; aoa=5, nbndlayers=8)            # RANS
msh = mshairfoil(:rae2822; aoa=2.79, hle=0.002)             # sample airfoils
msh = mshairfoil(naca4("2412"); aoa=4)                      # NACA 4-digit
msh = mshairfoil("myfoil.dat")                              # Selig or Lednicer file
```

The sizes of the gmsh mesh are keyword arguments; for example `hnear`,
`dnear`, `awake=1` and `dwake` give the isotropic elements that wall-resolved
LES needs (see `examples/naca/mknaca1msh.jl`). `airfoil_geo` returns the
generated `.geo` file for a look in the gmsh GUI.

## Periodic meshes

```julia
msh = mshsquare(8; periodic_dirs=(1,2))  # periodic in both x and y
```

## Refinement

Any 2D quad or triangle mesh can be refined, curved or not; the children
follow the parent's polynomial geometry.

```julia
msh = uniref(ex1mesh(eg=Simplex{2}()), 2)    # split every element into 4, twice

marked = falses(size(msh.nb))                # one flag per element edge
marked[:, 1:10] .= true
msh = refine(msh, marked)                    # local refinement, closed conformingly
```

Quad meshes can also be refined anisotropically towards a boundary, here
three boundary layers along an airfoil:

```julia
msh = mshairfoil(nbndlayers=0)               # requires gmsh on the PATH
msh = bndlayer_refine(uniref(msh), 1, 3)     # boundary tag 1, 3 layers
```
