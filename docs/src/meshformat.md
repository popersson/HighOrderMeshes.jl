# Mesh format

A `HighOrderMesh` is three arrays: `x` (global nodes, `nnodes × D`), `el`
(element-to-node table, `nnodes_per_elem × nel`), and `nb` (per-face
[`Neighbor`](@ref), `nfaces × nel`). DG nodes are obtained as `x[el, :]`, or
directly via [`dg_nodes`](@ref)`(m)`.

## Node order

The reference nodes of an element are in canonical order: multi-indices
`(i_1, ..., i_D)` with `i_1` varying fastest, `0 <= i_d <= p`, and for
simplices kept only when `sum(i) <= p`. This order never changes; every
import/export table (Gmsh, VTK, 3DG) is expressed relative to it.

## Face convention

Vertices of a `Block{D}` are in lexicographic tensor order (`i_1` varying
fastest); face `2d-1` is the side `ξ_d = 0` and face `2d` is the side
`ξ_d = 1`. For the unit square (`Block{2}`):

```
         face 4 (ξ_2 = 1)
    v3 ─────────────── v4
     │                  │
face │                  │ face 2
  1  │                  │ (ξ_1 = 1)
(ξ_1 │                  │
 = 0)│                  │
    v1 ─────────────── v2
         face 3 (ξ_2 = 0)
```

Faces of a `Simplex{D}` are numbered by the opposite vertex instead.
Conversions to the Gmsh and VTK orderings live only in `src/io/`.
