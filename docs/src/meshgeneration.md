# Mesh generation

Every generator returns a `HighOrderMesh`, ready to use, plot, refine or
save. Each line below can be copied as is; `p` is the polynomial degree and
`eg` the element geometry (blocks by default).

## Cheat sheet

```julia
using HighOrderMeshes

msh = mshsquare(8)                              # 8x8 quads on [0,1]^2
msh = mshsquare(8; eg=Simplex{2}(), p=3)        # 128 triangles, degree 3
msh = mshcube(4; eg=Simplex{3}())               # 384 tetrahedra on [0,1]^3
msh = mshline(10; p=4)                          # 10 line elements on [0,1]
msh = mshcircle(4; p=3)                         # curved quads on the unit disk
msh = mshcircle(4; eg=Simplex{2}(), p=3)        # curved triangles on the unit disk
msh = ex1mesh(eg=Simplex{2}())                  # sinusoidally curved sample mesh

# gmsh on the PATH:
msh = gmsh_sample(:circle; h=0.1, p=3)          # unstructured curved triangles
msh = gmsh_sample(:sphere; h=0.3, p=2)          # unstructured curved tetrahedra
msh = mshairfoil(:naca0012; aoa=5)              # curved quads around an airfoil
msh = rungmsh2msh("mygeometry.geo"; p=3)        # your own gmsh geometry

# DistMesh installed (] add DistMesh):
using DistMesh
fd(p) = sqrt(sum(p.^2)) - 1                     # signed distance to the unit circle
msh = HighOrderMesh(distmesh2d(fd, huniform, 0.2, ((-1,-1), (1,1))))
msh = curve_boundary(set_degree(msh, 4), fd)    # degree 4, boundary on the circle
```

To look at a mesh, load a Makie backend:

```julia
using CairoMakie            # or GLMakie for an interactive window
plot(msh)
```

The boundary faces carry numbers (`bndtag(msh.nb[j, iel])`) that solvers
use for boundary conditions; each generator below lists its numbering.

## Squares, cubes and lines

`mshsquare(m, n)`, `mshcube(m, n, o)` and `mshline(m)` mesh the unit
square, cube and interval with `m × n (× o)` cells. With `eg=Simplex{2}()`
every cell is split into 2 triangles, and with `eg=Simplex{3}()` into 6
tetrahedra. The boundaries are numbered `1..2D`: `xmin, xmax, ymin, ymax,
zmin, zmax`.

```julia
msh = mshsquare(8, 4; p=2)                      # 8x4 quads, degree 2
msh = mshsquare(8; eg=Simplex{2}())             # 128 triangles
msh = mshcube(3; eg=Simplex{3}(), p=2)          # 162 tetrahedra, degree 2
msh = mshsquare(8; periodic_dirs=(1,2))         # periodic in both x and y
```

## Disk

`mshcircle(n)` meshes the unit disk, or with `shape=:half` or
`shape=:quarter` the upper half or first quadrant. The nodes on the
circle are exactly on it, at equal angles. Quads come from a transfinite map
of 12 patches (`12 n^2` elements), triangles from a hexagon of 6 sectors
(`6 n^2` elements) mapped smoothly onto the disk. The circle is boundary 1;
`:half` numbers `y=0, arc` and `:quarter` numbers `y=0, x=0, arc`.

```julia
msh = mshcircle(3; p=3)                         # 108 curved quads
msh = mshcircle(4; eg=Simplex{2}(), p=3)        # 96 curved triangles
msh = mshcircle(2; shape=:quarter, p=4)         # 12 quads on the quarter disk
```

## Gmsh samples

`gmsh_sample(shape; h, p, eg)` makes unstructured meshes of size about `h`
with gmsh, which must be on the system `PATH`. On the circle and the sphere
gmsh places the high-order nodes on the exact geometry.

| `shape`   | Domain    | Boundary numbers                        |
|:--------- |:--------- |:--------------------------------------- |
| `:square` | `[0,1]^2` | `xmin, xmax, ymin, ymax`, as `mshsquare` |
| `:circle` | unit disk | 1: circle                                |
| `:cube`   | `[0,1]^3` | `xmin, xmax, ..., zmax`, as `mshcube`    |
| `:sphere` | unit ball | 1: sphere                                |

```julia
msh = gmsh_sample(:square; h=0.1)                       # triangles
msh = gmsh_sample(:circle; h=0.2, p=3, eg=Block{2}())   # curved quads
msh = gmsh_sample(:cube; h=0.3, eg=Block{3}())          # hexahedra
msh = gmsh_sample(:sphere; h=0.3, p=2)                  # curved tetrahedra
```

For other geometries, write a `.geo` file and import it with
[`rungmsh2msh`](@ref), see [Gmsh import](@ref).

## Curved meshes from a distance function

Any straight-sided 2D mesh of triangles or quads becomes a curved
high-order mesh in two steps: [`set_degree`](@ref) adds the high-order
nodes, and [`curve_boundary`](@ref) moves the boundary faces onto the zero
level set of a signed distance function `fd` (negative inside), the way
DistMesh describes geometries. The nodes along each boundary face are
placed at equal arclengths (or at the arclength fractions of the reference
nodes, for example Gauss-Lobatto nodes), and the elements next to the
boundary bend with their faces.

With [DistMesh](https://github.com/JuliaGeometry/DistMesh.jl), which is
optional, `HighOrderMesh(dm)` converts its meshes:

```julia
using HighOrderMeshes, DistMesh

# Unit circle
fd(p) = sqrt(sum(p.^2)) - 1                     # or dcircle(p)
fh(p) = 1.0                                     # or huniform(p)
msh = HighOrderMesh(distmesh2d(fd, fh, 0.2, ((-1,-1), (1,1))))
msh = curve_boundary(set_degree(msh, 4), fd)

# Square plate with a hole, finer at the hole; the straight sides do not move
fd(p) = ddiff(drectangle(p, -1, 1, -1, 1), dcircle(p; r=0.5))
fh(p) = 0.05 + 0.3*dcircle(p; r=0.5)
corners = ((-1,-1), (-1,1), (1,-1), (1,1))
msh = HighOrderMesh(distmesh2d(fd, fh, 0.05, ((-1,-1), (1,1)), corners))
msh = curve_boundary(set_degree(msh, 3), fd)

# Ellipse: any function with the right zero level set works
fd(p) = (p[1]/2)^2 + p[2]^2 - 1
msh = HighOrderMesh(distmesh2d(fd, huniform, 0.2, ((-2,-1), (2,1))))
msh = curve_boundary(set_degree(msh, 3), fd)
```

`fd` is called with a point `p` as an `SVector`, so `p[1]`, `p[2]`,
`norm(p)` and the DistMesh distance functions all work. Without DistMesh,
pass node coordinates and an element table to the constructor,
`HighOrderMesh(x, el)` (`x` is `nnodes × 2`, `el` is `3 × nel` or `4 × nel`).

To curve only some boundaries, pass their numbers; here the arc of a
quarter disk with straight edges:

```julia
msh = set_degree(mshcircle(2; shape=:quarter), 3)   # p=1 mesh: straight edges
msh = curve_boundary(msh, p -> sqrt(sum(p.^2)) - 1, 3)
```

The building blocks are available on their own: [`project_points`](@ref)
projects points onto the zero level set, and [`interp_arclength`](@ref)
redistributes points along a polyline by arclength.

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
