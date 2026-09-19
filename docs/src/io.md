# File I/O

Meshes are saved and loaded with [`savemesh`](@ref) and [`loadmesh`](@ref) in
a custom binary `.hom` format; solutions can additionally be exported to VTK
for viewing in ParaView or similar tools.

## The `.hom` format

A `.hom` file is a sequence of named arrays ("records") stored as raw
little-endian binary. It is self-describing and readable from any language
with a few lines of code:

```
[8 bytes]         magic string "HOMESHv1"
[Int64]           number of records
then, for each record:
  [Int64]         length of the name in bytes
  [bytes]         name (ASCII)
  [Int64]         element type code (see hom_type_codes)
  [Int64]         number of dimensions (0 for a scalar)
  [Int64 × ndims] shape
  [data]          elements in column-major order
```

Type codes: 1 Float64, 2 Float32, 3 Float16, 11 Int64, 12 Int32, 13 Int16, 14 UInt8.

Records written by `savemesh` (`D` = spatial dimension, `ns` = nodes per
element, `nf` = faces per element, `T` = coordinate type):

| Record      | Type    | Shape          | Contents                                             |
|:----------- |:------- |:-------------- |:----------------------------------------------------- |
| `elgeom`    | Int64   | scalar         | `1` = Simplex, `2` = Block                            |
| `ref_nodes` | `T`     | `ns × D`       | reference nodes of the `FiniteElement`                |
| `x`         | `T`     | `nnodes × D`   | node coordinates                                      |
| `el`        | Int64   | `ns × nel`     | element-to-node table (one-based)                     |
| `nb_el`     | Int32   | `nf × nel`     | neighbor element (`<= 0`: boundary face with tag `-value`) |
| `nb_face`   | Int16   | `nf × nel`     | neighbor face index                                   |
| `nb_perm`   | Int16   | `nf × nel`     | neighbor face permutation                             |

The polynomial degree follows from the number of reference nodes. Readers
must ignore records they do not know, so that records can be added later.

```julia
msh = mshsquare(8)
savemesh("square.hom", msh)
msh2 = loadmesh("square.hom")
```

## VTK export

```julia
msh = mshcircle(2, p=2)
u = hcat(msh.x[:,1].^2, msh.x[:,2], -msh.x[:,1])   # 3 components
vtkwrite("out.vtk", msh, u, umap=[1, 2:3])         # scalar u1, vector u23
```
