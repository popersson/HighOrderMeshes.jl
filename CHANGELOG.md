# Changelog

## Unreleased

Isoparametric h-refinement of curved 2D quad and triangle meshes, ported from
3DG's `qmshrefine`, `qmshuniref` and `qmshbndlayer`.

- `refine(m, marked)` refines the edges marked in an `nfaces × nel` `Bool`
  matrix, extending the marks so the result is conforming. Quads split into
  2, 3 or 4 children, triangles into 2 (bisection) or 4.
  `refine_with_parents` also returns each element's parent index.
- `uniref(m, nref)` now works for any degree and node set and refines curved
  elements isoparametrically (it used to require `p=1`). Element and node
  order of the result changed.
- `bndlayer_refine(m, bnds, nlayers)` refines quad meshes anisotropically
  towards the boundaries `bnds`; `bndlayer_refine_with_elements` also returns
  the indices of the boundary layer elements.
- Example `examples/naca/mknaca1msh.jl`: NACA 0012 mesh with boundary layers.
- Fix: Gmsh physical names containing spaces were cut at the first space.

## v0.3.1

Support for solvers that work on their own, non-conforming node sets, such
as DG-SEM on Gauss-Legendre or Gauss-Radau points. The mesh keeps its
conforming nodes; the solver builds a second reference element and moves
data between the two.

- `interpolate(fe_from, fe_to, u)` re-nodalizes a field between two reference
  elements on the same geometry; exact when the degrees agree.
- `FiniteElement(eg::Block, s1; check=false)` builds a block element from a
  non-conforming 1D line, matching the existing `check` keyword of the
  node-matrix constructor.
- `viz_solution(m, u; fe)` and the Makie solution plot attribute `fe` plot a
  field straight from a solver's node layout.
- `gauss_radau_nodes`, `gauss_radau01_nodes`, `gauss_radau_quadrature` and
  `gauss_radau01_quadrature` (left Gauss-Radau rules, exact to degree `2n-2`).
- New documentation page "Solver node sets" with the complete workflow.

## v0.3.0

Breaking redesign of the core types and the polynomial tools. The mesh
format (`x`, `el`, `nb`), the basic mesh constructors and the hello world
are unchanged; most other public names changed. The last release before the
redesign is tagged `v0.2.0`.

### Types

- `HighOrderMesh{D,G,T}` and `FiniteElement{D,G,T}`: the polynomial degree is
  a value (`porder(m)`), no longer a type parameter, so `set_degree` and the
  element constructors are type-stable and no code is recompiled per degree.
- `FiniteElement` stores the degree, the full `D`-dimensional reference node
  matrix and the inverse Vandermonde matrix. Sub-dimensional nodes and
  elements are derived on demand with `ref_nodes(fe, d)` and
  `subelement(fe, d)`.
- Reference node sets are validated by `check_conforming`: the vertices must
  be present, the set must be invariant under `symmetries(eg)`, unisolvent,
  and conforming on the faces. Non-symmetric node lines are rejected, so
  meshes are always continuous across shared faces, edges and corners.
- `Neighbor` struct with fields `el`, `face`, `perm` replaces the neighbor
  tuple, with `isboundary` and `bndtag` accessors.
- 1D meshes always use `Block{1}`; `Simplex{1}` is no longer supported, and
  `subgeom` of a simplex returns `Block{d}` for `d <= 1`.

### Polynomial tools

- Scalar `jacobi`, `djacobi`, `legendre`, `dlegendre`, `legendre01` and
  `dlegendre01` replace the vectorized `legendre_poly` family; Vandermonde
  matrices are comprehensions. Exact number types are preserved.
- `polybasis` and `dpolybasis` replace `eval_poly` and the `multivar_*`
  functions. Simplices use an orthogonal Koornwinder basis instead of
  monomials (Vandermonde condition number at degree 5 down from 4e4 to 22 on
  triangles and from 1.5e5 to 73 on tetrahedra).
- `equispaced_nodes(eg, p)`, `tensor_nodes(eg, s1)` and `multiindices(eg, p)`
  define the canonical node order, which is unchanged.
- `quadrature(eg, p, T=Float64)`: `p` is the degree of exactness for both
  geometries (blocks previously took points per direction), and weights sum
  to the reference volume for both (simplex weights previously summed to 1).
  Simplex rules live in `src/basis/simplex_quadrature_rules.txt`.
- The number type is a positional argument, `gauss_legendre_nodes(n, Float32)`.

### Shape functions

- `shapefcns(fe, ξ)`, `dshapefcns(fe, ξ)` and `interpolate(N, u)` replace
  `eval_shapefcns` and `eval_field`. No keyword changes a return type.

### Mesh utilities

- `set_ref_nodes` accepts a `FiniteElement`, a node matrix, or for blocks a
  1D node line, and returns a mesh that shares no arrays with the input.
- `facemap` and `edgemap` for blocks are generated from the tensor
  structure. In 2D `edgemap == facemap`; 3D block edges are numbered per
  coordinate direction.
- `vertices(eg)`, `nnodes(eg, p)` and `symmetries(eg)` added.
- `unique_mesh_nodes(x, el; tol)` exposes the deduplication tolerance.

### I/O

- The `.hom` format is redefined as self-describing named records with the
  magic string `HOMESHv1`; files written by v0.2.0 cannot be read. The text
  format is removed.
- `gmsh_physical_names(fname)` added; `gmsh2msh(fname; verbose=true)` prints
  the boundary tags it found.

### Visualization

- The Makie extension is rewritten on the recipe system, so `plot`, `plot!`,
  `plot(m, u)` and `plot!(m, u)` all work, and it carries a precompile
  workload. The Plots.jl extension is removed.

### Removed

- `legendre_poly`, `legendre01_poly`, `multivar_legendre01_poly`,
  `multivar_monomial_poly`, `eval_poly`, `eval_shapefcns`, `eval_field`,
  `nbr_ho_nodes`, `SimplexQuadRule`, `savemeshtxt`, `loadmeshtxt`,
  `ex_plot_mesh`, `ex_plot_solution`, `ex_vtkexport`, `ex_read_gmsh`.
- No longer exported: `name`, `el2nb`, `equispaced`, `blockmesh_hypercube`.

### Documentation

- Documenter site with quick start, mesh format, basic meshes, Gmsh import,
  file I/O and API reference pages, deployed by the `Docs` workflow.

## v0.2.0

Last release before the v0.3 redesign: reorganized directory structure,
binary `.hom` format, Makie and Plots extensions, 3DG import and export.

## v0.1.0

Original 2024 to 2025 code base.
