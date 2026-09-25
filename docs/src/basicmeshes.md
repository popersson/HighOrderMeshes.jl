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

Quad meshes can also be refined anisotropically towards a boundary. The NACA
example in `examples/naca/mknaca1msh.jl` adds three boundary layers along the
airfoil:

```julia
msh = rungmsh2msh("examples/naca/naca.geo"; porder=3)   # requires gmsh on the PATH
msh = set_lobatto_nodes(msh)
msh = bndlayer_refine(uniref(msh), 1, 3)     # boundary tag 1, 3 layers
```
