# Gmsh import

Meshes can be generated with [Gmsh](https://gmsh.info) and imported as a
`HighOrderMesh`. Boundary tags come from Gmsh's physical groups, so every
`.geo` file needs a `Physical` line per boundary region you want to identify.

```
// square.geo
Point(1) = {0, 0, 0}; Point(2) = {1, 0, 0};
Point(3) = {1, 1, 0}; Point(4) = {0, 1, 0};
Line(1) = {1, 2}; Line(2) = {2, 3}; Line(3) = {3, 4}; Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Physical Curve("left", 1) = {4};
Physical Curve("rest", 2) = {1, 2, 3};
Physical Surface("domain", 1) = {1};
```

```julia
using HighOrderMeshes
msh = rungmsh2msh("square.geo"; porder=3)
```

`rungmsh2msh` runs the `gmsh` executable (it must be on the system `PATH`)
and imports the result; `gmsh2msh` imports an existing `.msh` v2.2 file
directly and prints the boundary tags found (`boundary 1: "left"`, ...)
unless called with `verbose=false`. Use [`gmsh_physical_names`](@ref) to get
the tag-to-name mapping programmatically.
