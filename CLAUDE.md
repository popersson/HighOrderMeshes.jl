# HighOrderMeshes.jl

Light-weight Julia library for high-order unstructured meshes on simplex and
block elements. Single element type and single polynomial degree per mesh.

## Commands

- Run the tests: `julia --project -e 'using Pkg; Pkg.test()'` (about 90 s).
- Plotting tests (needs CairoMakie, slow): see the header of `test/plotting/runtests.jl`.
- Hello world:
  ```julia
  using HighOrderMeshes
  msh = ex1mesh(eg=Simplex{2}())
  msh = set_degree(msh, 3)
  ```
- Gmsh is installed on the development machine and used through the command line.

## Layout

- `src/mesh/`: element geometry, the mesh struct, refinement, basic meshes.
- `src/basis/`: polynomials, quadrature, the reference element.
- `src/io/`: `.hom` binary format, Gmsh import, VTK and 3DG export.
- `src/viz/`: backend-independent plotting geometry; Makie extension in `ext/`.
- `src/fem/`: prototype CG assembly. Do not redesign; only rename what a core change requires.
- `docs/dev/redesign-plan.md`: the current work plan with per-step instructions.

## Design rules (from the author; do not reverse)

- No new package dependencies. Optional integrations go in package extensions.
- The mesh format is `x` (global nodes, nnodes × D), `el` (element-to-node table,
  nnodes_per_elem × nel), `nb` (per-face neighbor data, nfaces × nel). Short names stay.
  DG nodes are `x[el, :]`. Do not add fields to the mesh without a strong reason.
- The mesh node set is always continuous and conforming: reference node sets
  must be symmetric so nodes coincide across shared faces, edges and corners.
  Solver-specific node sets belong in the solver, not in the mesh.
- Vertex order is lexicographic tensor order; faces of blocks are numbered per
  coordinate (face 2d-1 is ξ_d = 0, face 2d is ξ_d = 1); simplex faces are
  numbered by the opposite vertex. Conversions to Gmsh and VTK orderings live
  only in `src/io/`.
- The canonical node order (multi-indices with the first index varying
  fastest, simplices filtered by `sum(i) <= p`) must never change; every
  import and export table depends on it.
- 1D meshes always use `Block{1}`. `Simplex{1}` is not supported.
- Type parameters are for dispatch only: `{D, G, T}`. The polynomial degree is
  a value, never a type parameter.
- Functions return one type. Never use a keyword flag that changes the return
  type; write a separate function instead (`shapefcns` and `dshapefcns`,
  `legendre` and `dlegendre`). Derivatives use a `d` prefix.
- Prefer scalar functions plus broadcasting or comprehensions over
  MATLAB-style vectorized functions in the polynomial tools. Performance is
  not critical there.
- Coordinates are stored with points as rows and the coordinate index last
  (`x` is nnodes × D, `dg_nodes` is ns × nel × D, gradients put D last).
- Meshes are values: functions like `set_degree` return new meshes and share
  no arrays with the input. The bang functions `set_bnd_numbers!`,
  `set_bnd_periodic!`, `align_with_ldgswitch!` mutate `nb` or `el` in place.
- Keep the `.hom` format trivially readable from C: named records of raw
  little-endian arrays, documented in the header of `src/io/io.jl`.

## Style

- Docstrings on every exported function, in the existing style.
- Tests in `test/runtests.jl` inside the existing testsets or a new one per
  feature; every numerical routine gets a correctness test, not only a smoke test.
- Commit messages: `Step N: <title>` while working through the redesign plan.
