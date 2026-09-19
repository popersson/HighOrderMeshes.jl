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
