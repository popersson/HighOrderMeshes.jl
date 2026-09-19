# HighOrderMeshes.jl v0.3 redesign plan

This document is the working plan for the v0.3 redesign. It is written so that
an implementing model can execute the open steps one at a time without further
design decisions. Read `CLAUDE.md` in the repository root first; it states the
project conventions that every step must respect.

Decisions in this plan were made together with the package author and are
final. Do not reopen them. If a step is ambiguous, choose the simplest
interpretation that keeps the tests green and mention the choice in the commit
message.

## Rules of engagement for the implementing model

1. Work on the branch `redesign-v0.3`. Do one step per session. Never combine
   steps in one commit.
2. Before starting a step, run the test suite once to confirm it is green:
   `julia --project -e 'using Pkg; Pkg.test()'` (about 90 s).
3. Implement the step exactly as specified, including the listed tests.
4. Run the test suite again. It must be green before committing. If a test
   in an unrelated area breaks, fix the cause, not the test.
5. Commit with the message `Step N: <title from this document>` and a short
   body listing anything you had to decide yourself.
6. Update the status table below in the same commit.
7. Do not add package dependencies. Do not touch `src/fem/` beyond what a
   rename requires. Do not push.

## Status

| Step | Title | Owner | Status |
|---|---|---|---|
| 0 | Tag v0.2.0, create branch, write plan and CLAUDE.md | Fable | done |
| 1 | Polynomial tools: scalar functions, Koornwinder basis, quadrature semantics | Fable | done |
| 2 | FiniteElement redesign: node matrix, conformity check, sub-elements, shapefcns/interpolate | Fable | done |
| 3 | Degree out of the type parameters, API propagation, .hom v1 redefined, Plots extension removed | Fable | done |
| 4 | Neighbor struct | implementer | done |
| 5 | Generated Block face and edge tables, edge renumbering, uniref update | implementer | done |
| 6 | Simplex quadrature rules as a text data file | implementer | done |
| 7 | Makie extension on the recipe system, precompile workload | implementer | done |
| 8 | Gmsh physical names, export cleanup, sample meshes into src | implementer | done |
| 9 | Documenter skeleton and docs workflow | implementer | done |
| 10 | Remaining tests and small cleanups | implementer | done |
| 11 | Version bump, changelog, final review | Fable | open |

## Target API after the redesign

Names below are final. Anything not listed keeps its current name.

### Element geometry (`src/mesh/element_geometry.jl`)

```julia
abstract type ElementGeometry{D} end
struct Simplex{D} <: ElementGeometry{D} end     # D >= 2 only; 1D meshes are Block{1}
struct Block{D}   <: ElementGeometry{D} end     # D >= 0

dim(eg), name(eg), nvertices(eg), nfaces(eg), nedges(eg)
vertices(eg)          # nvertices × D matrix of reference vertex coordinates
facemap(eg)           # nfv × nfaces, column j = local vertex indices of face j
edgemap(eg)           # 2 × nedges
subgeom(eg, d)        # sub-geometry of dimension d; Simplex gives Block{d} for d <= 1
midpoint(eg)          # centroid in reference coordinates
nnodes(eg, p)         # number of Lagrange nodes at degree p
symmetries(eg)        # vector of (A, b): affine maps ξ -> Aξ + b of the reference element onto itself
find_elgeom(D, nv)    # internal
plot_face_order(eg)   # internal, used by viz
```

Vertex order is lexicographic tensor order for `Block` (first coordinate
varies fastest) and origin followed by the unit vectors for `Simplex`. Faces
of `Block{D}` are numbered per coordinate: face `2d-1` is the side `ξ_d = 0`
and face `2d` is the side `ξ_d = 1`. Faces of `Simplex{D}` are numbered by the
opposite vertex.

### Polynomial tools (`src/basis/poly_tools.jl`)

```julia
jacobi(n, α, β, x)      # scalar Jacobi polynomial P_n^{(α,β)}(x) on [-1,1]
djacobi(n, α, β, x)     # its derivative
legendre(n, x)          # P_n(x) on [-1,1]
dlegendre(n, x)
legendre01(n, x)        # shifted to [0,1]
dlegendre01(n, x)

polybasis(eg, ξ, p)     # nξ × nbasis matrix; ξ is nξ × D
dpolybasis(eg, ξ, p)    # nξ × nbasis × D
                        # Block: tensor product of shifted Legendre polynomials
                        # Simplex: orthogonal (unnormalized) Koornwinder basis

equispaced(p)                    # (0:p) // p
equispaced_nodes(eg, p)          # nnodes × D matrix of rationals, in the canonical node order
tensor_nodes(eg::Block, s1)      # tensor product of a 1D node line

gauss_legendre_nodes(n, T=Float64), gauss_legendre_quadrature(n, T=Float64), and the 01 and lobatto variants
quadrature(eg, p, T=Float64)     # (ξ, w) exact for polynomials of degree p, for both geometries
                                 # (the number type is a positional argument: a type passed as a keyword
                                 #  is not visible to inference)
```

Canonical node order: iterate multi-indices `(i_1, ..., i_D)` with `i_1`
varying fastest, `0 <= i_d <= p`, and for simplices keep only `sum(i) <= p`.
Node coordinates are `i ./ p`. This order must never change; the Gmsh, VTK
and 3DG permutation tables are expressed relative to it.

### Reference element (`src/basis/finite_element.jl`)

```julia
struct FiniteElement{D,G<:ElementGeometry{D},T}
    p::Int              # polynomial degree
    nodes::Matrix{T}    # nnodes × D reference nodes
    coeff::Matrix{T}    # inverse Vandermonde: shapefcns(ξ) = polybasis(ξ) * coeff
end

FiniteElement(eg, nodes::AbstractMatrix; check=true)   # explicit node set (degree inferred)
FiniteElement(eg, p::Int, T=Float64)                    # equispaced nodes
FiniteElement(eg::Block, s1::AbstractVector)            # tensor product of a 1D line

elgeom(fe), dim(fe), porder(fe), nnodes(fe), name(fe)
ref_nodes(fe)          # the node matrix
ref_nodes(fe, d)       # nodes on the canonical d-dimensional sub-face ξ_{d+1..D} = 0, as an n × d matrix
subelement(fe, d)      # FiniteElement of subgeom(eg, d) with those nodes
corner_nodes(fe)       # indices of the reference vertices within the node set
check_conforming(eg, nodes)   # throws ArgumentError if the node set cannot be used for a conforming mesh

shapefcns(fe, ξ)       # nξ × nnodes
dshapefcns(fe, ξ)      # nξ × nnodes × D
shapefcns(eg, ξ), dshapefcns(eg, ξ)     # same, for the p=1 element on eg
interpolate(N, u)      # N is nξ × ns from shapefcns; u is ns × ...; result nξ × ...
interpolate(dN, u)     # dN is nξ × ns × D from dshapefcns; result nξ × ... × D
```

A node set is conforming when it has the right number of nodes for some
degree, contains the reference vertices, is invariant under every map in
`symmetries(eg)`, is unisolvent, and its restriction to the canonical
`(D-1)`-face is itself conforming for the sub-geometry. The constructor
checks this and throws `ArgumentError` otherwise.

### Mesh (`src/mesh/high_order_mesh.jl`)

```julia
struct HighOrderMesh{D,G<:ElementGeometry{D},T}
    fe::FiniteElement{D,G,T}
    x::Matrix{T}            # nnodes × D
    el::Matrix{Int}         # nnodes_per_elem × nel
    nb::Matrix{Neighbor}    # nfaces × nel   (Neighbor struct arrives in step 4)
end

HighOrderMesh(fe, x, el; bndexpr), HighOrderMesh(x, el; bndexpr)
elgeom(m), dim(m), porder(m), nnodes(m), nel(m), dg_nodes(m)
set_ref_nodes(m, fe), set_ref_nodes(m, nodes::AbstractMatrix), set_ref_nodes(m::Block mesh, s1::AbstractVector)
set_degree(m, p), set_lobatto_nodes(m)
mkface2nodes(m), mkface2nodes(fe)
```

`set_ref_nodes` and its wrappers return a new mesh that shares nothing with
the old one.

### Neighbor (step 4)

```julia
struct Neighbor
    el::Int32      # > 0: neighbor element; <= 0: boundary face with tag -el
    face::Int16    # local face index in the neighbor element, 0 for boundary faces
    perm::Int16    # one-based face permutation (1 in 1D and 2D), 0 for boundary faces
end
isboundary(nb), bndtag(nb)
```

### File format (`src/io/io.jl`)

`savemesh(fname, m)` and `loadmesh(fname)` use the `.hom` format defined in
the header comment of `io.jl`: a magic string, a record count, then named
arrays each stored as name, type code, rank, shape, raw little-endian data.
The reader ignores unknown records. There is no text format any more.

## Steps

### Step 1. Polynomial tools (done)

Rewrote `poly_tools.jl` around scalar functions. `jacobi` uses the standard
three-term recurrence, `djacobi(n, α, β, x) = (n+α+β+1)/2 * jacobi(n-1, α+1, β+1, x)`.
The block basis is the tensor product of shifted Legendre polynomials in the
canonical multi-index order. The simplex basis is the orthogonal Koornwinder
basis in collapsed coordinates on the `[0,1]` simplex, following the
Hesthaven and Warburton formulas without the normalization constants so that
exact number types keep working. Gradients use the singularity-free form of
the same reference. `quadrature(eg::Block, p)` now uses `cld(p+1, 2)` Gauss
points per direction so that `p` means degree of exactness for both
geometries. The tabulated simplex weights sum to 1; `quadrature` rescales
them by `1/factorial(D)` so that weights sum to the reference volume for
both geometries, which is what `FEM_precomp` assumes when it multiplies by
the Jacobian. `find_elgeom(D, P, nnodes)` was removed. The number type is a
positional argument (`gauss_legendre_nodes(n, Float32)`) because a type
passed as a keyword is invisible to inference.

Findings about the tabulated simplex rules, recorded for step 6: the
triangle rule for degree 19 is accurate to about 1.4e-12 while all others
reach 1e-13; the tetrahedron rules for degrees 14 and 15 (330 points) contain
a negative weight of about -1.2 relative to a total of 1, which is a poor
rule. Replacing those two with positive-weight rules from the literature is
a good follow-up once the rules live in a data file, but it is not part of
step 6.

### Step 2. FiniteElement redesign (done)

Implemented as specified in the target API. `symmetries(eg)` generates the
hyperoctahedral group for blocks and the barycentric permutation group for
simplices. Sub-elements are derived on demand from the node matrix. The
constructor validates conformity by default; `check=false` skips it for
internally derived elements.

### Step 3. Degree out of the type, API propagation, .hom, Plots removal (done)

`HighOrderMesh{D,G,T}` and `FiniteElement{D,G,T}`. All call sites moved to
`shapefcns`, `dshapefcns`, `interpolate`, `ref_nodes`, `subelement`. The
`.hom` format was redefined as named records and the example mesh was
converted. The Plots extension, its weak dependencies and its tests were
removed. `Simplex{1}` methods were removed; `subgeom` of a simplex returns
`Block{d}` for `d <= 1`.

### Step 4. Neighbor struct

Files: `src/mesh/high_order_mesh.jl`, `src/mesh/mesh_utils.jl`,
`src/io/io.jl`, `src/io/converters.jl`, `src/viz/post_processing.jl`,
`src/fem/dg_utils.jl`, `src/HighOrderMeshes.jl`, `test/runtests.jl`.

1. In `high_order_mesh.jl` replace `const NeighborData = Tuple{Int32,Int16,Int16}`
   with the `Neighbor` struct from the target API, a docstring, and
   ```julia
   isboundary(nb::Neighbor) = nb.el <= 0
   bndtag(nb::Neighbor) = -Int(nb.el)
   Base.show(io::IO, nb::Neighbor) = isboundary(nb) ?
       print(io, "Neighbor(boundary $(bndtag(nb)))") :
       print(io, "Neighbor(el=$(nb.el), face=$(nb.face), perm=$(nb.perm))")
   ```
   The default constructor already converts integer arguments, so
   `Neighbor(2, 1, 3)` works.
2. Replace every tuple construction `(a, b, c)` stored into an `nb` matrix by
   `Neighbor(a, b, c)`, and every positional access `nb[1]`, `nb[2]`, `nb[3]`
   or destructuring `jel, _, _ = m.nb[j,iel]` by field access `.el`, `.face`,
   `.perm`. Use `isboundary` where the code tests `nb[1] < 1` or
   `nb[1] <= 0`, and `bndtag` where it negates `nb[1]`. Grep for `nb[`,
   `.nb`, `first.(`, `last.(` and `NeighborData` to find all sites, including
   the tests.
3. In `io.jl` nothing changes in the file layout; the writer builds the three
   integer arrays with `getfield.(m.nb, :el)` and so on, and the reader
   builds `Neighbor.(nb_el, nb_face, nb_perm)`.
4. Export `Neighbor`, `isboundary`, `bndtag`.
5. Tests to add in a new testset "Neighbor":
   ```julia
   @test Neighbor(2, 1, 3) == Neighbor(2, 1, 3)
   @test isboundary(Neighbor(-2, 0, 0)) && bndtag(Neighbor(-2, 0, 0)) == 2
   @test !isboundary(Neighbor(5, 1, 1))
   @test sizeof(Neighbor) == 8
   m = mshsquare(2)
   @test count(isboundary, m.nb) == 8
   ```
   Update the existing tests that compare `nb` entries with tuples.

### Step 5. Generated Block face and edge tables, edge renumbering, uniref update

Files: `src/mesh/element_geometry.jl`, `src/mesh/mesh_utils.jl`, `test/runtests.jl`.

1. Before changing anything, copy the current hand-written `facemap` and
   `edgemap` matrices for `Block{1}`, `Block{2}`, `Block{3}`, `Simplex{2}`,
   `Simplex{3}` into the test file as literal expected values (see item 6).
2. Replace the hand-written `facemap(::Block{D})` methods by one generated
   method, built by extrusion:
   ```julia
   function facemap(::Block{D}) where {D}
       D == 0 && return zeros(Int, 1, 0)
       D == 1 && return [1 2]                     # face 1: ξ=0, face 2: ξ=1 (a 1 × 2 matrix)
       lower = facemap(Block{D-1}())              # nfv_lower × 2(D-1)
       nvl   = 2^(D-1)
       extruded = vcat(lower, lower .+ nvl)        # each lower face extruded along ξ_D
       bottom = _oriented_block_face(D, 0)         # side ξ_D = 0
       top    = _oriented_block_face(D, 1)         # side ξ_D = 1
       hcat(extruded, bottom, top)
   end
   ```
   where `_oriented_block_face(D, side)` returns the vertex list of
   `Block{D-1}` in tensor order, shifted by `side * 2^(D-1)`, with the last
   local axis flipped when the orientation determinant is negative. The
   orientation determinant is `det([n; a_1; ...; a_{D-1}])` where `n` is the
   outward unit normal of that side and `a_k` is the direction of local axis
   `k` of the face, all as rows. Flipping the last local axis means reversing
   the order of the `2^(D-2)`-blocks of the vertex list. This rule reproduces
   the existing tables for D = 2 and 3 exactly; the regression test in item 6
   enforces it.
3. Replace the hand-written `edgemap(::Block{D})` methods:
   ```julia
   edgemap(eg::ElementGeometry{2}) = facemap(eg)     # in 2D edges are faces
   edgemap(::Block{1}) = [1 2]
   ```
   and for `D >= 3` generate edges along direction `d = 1, ..., D`, for each
   position of the other `D-1` coordinates in lexicographic order with the
   lowest coordinate varying fastest, oriented from `ξ_d = 0` to `ξ_d = 1`.
   For `Block{3}` this gives
   `[1 3 5 7 1 2 5 6 1 2 3 4; 2 4 6 8 3 4 7 8 5 6 7 8]`.
   `Simplex` edge tables stay hand-written; `Simplex{2}` already satisfies
   `edgemap == facemap` and `Simplex{3}` keeps the lexicographic pairs.
4. Keep `facemap(::Simplex{2})` and `facemap(::Simplex{3})` hand-written; the
   `Simplex{3}` order matches the 3DG convention and must not change. Add a
   test that every simplex face is the complement of its opposite vertex and
   is outward oriented (the sign of `det([a_1; ...; a_{D-1}; n])` with
   `a_k = v_{k+1} - v_1` and `n` the outward normal of that face is
   positive).
5. Update the `Block` branch of `uniref` in `mesh_utils.jl`, which currently
   assumes the cyclic edge order. With vertices `v1..v4` in tensor order,
   face midpoints `m1..m4` (face `j` midpoint is `mapmid[j, it]`) and centroid
   `c`, the four sub-quads in tensor order are
   ```
   [v1, m3, m1, c]   [m3, v2, c, m2]   [m1, c, v3, m4]   [c, m2, m4, v4]
   ```
   Keep the `Simplex` branch as is.
6. Tests. Add a testset "Generated geometry tables":
   - `facemap(Block{2}()) == [3 2 1 4; 1 4 2 3]` and
     `facemap(Block{3}()) == [3 2 1 4 3 5; 1 4 2 3 4 6; 7 6 5 8 1 7; 5 8 6 7 2 8]`
     (these are the current tables; copy them verbatim from the file before
     editing).
   - `edgemap(Block{2}()) == facemap(Block{2}())`, `edgemap(Simplex{2}()) == facemap(Simplex{2}())`.
   - `edgemap(Block{3}())` equals the matrix in item 3.
   - For `D in 1:4`: `size(facemap(Block{D}())) == (2^(D-1), 2D)` and every
     column lists exactly the vertices with `ξ_d = side`.
   - Orientation test for every face of `Block{2}`, `Block{3}`, `Simplex{2}`,
     `Simplex{3}` as described in item 4 (for blocks use the vertex list as
     `v_1, v_2, ..., v_{2^{D-1}}` and axes `a_1 = v_2 - v_1`, `a_2 = v_3 - v_1`).
   - `uniref(mshsquare(2), 1)` gives 16 elements, all with positive
     orientation, and `set_degree(uniref(mshsquare(2), 2), 2)` has the same
     boundary node count as `set_degree(mshsquare(8), 2)`.

### Step 6. Simplex quadrature rules as a text data file

Files: `src/basis/simplex_quadrature.jl` (to be deleted), a new
`src/basis/simplex_quadrature_rules.txt`, `src/basis/poly_tools.jl`,
`src/HighOrderMeshes.jl`, `test/runtests.jl`.

1. Write a one-off script (do not commit it) that loads the package, iterates
   `D in 2:3` and `p` from 1 upward until `simplex_quadrature` errors, and
   writes each rule to the text file in this format:
   ```
   # Symmetric quadrature rules on the reference simplex {ξ >= 0, sum(ξ) <= 1}.
   # Each block: "rule D p n" followed by n lines "ξ_1 ... ξ_D w".
   rule 2 1 1
   0.333333333333333 0.333333333333333 0.5
   rule 2 2 3
   ...
   ```
   Write numbers with `repr` so that no digits are lost. Rules for `D = 1`
   are not written. Write the weights as `quadrature(Simplex{D}(), p)`
   returns them, i.e. summing to the reference simplex volume `1/factorial(D)`
   (the current Julia tables sum to 1 and `quadrature` rescales them; after
   this step no rescaling happens in `quadrature`).
2. Delete `simplex_quadrature.jl` and its include. In `poly_tools.jl` add
   ```julia
   const _simplex_rules = Dict{Tuple{Int,Int},Tuple{Matrix{Float64},Vector{Float64}}}()
   function _load_simplex_rules!()   # parses the text file once, fills the Dict
   function quadrature(::Simplex{D}, p::Integer, ::Type{T}=Float64) where {D,T}
       isempty(_simplex_rules) && _load_simplex_rules!()
       haskey(_simplex_rules, (D, p)) || error("No simplex quadrature rule for D=$D, p=$p; maximum degree is $(max_quadrature_degree(Simplex{D}()))")
       ξ, w = _simplex_rules[(D, p)]
       T.(ξ), T.(w)
   end
   max_quadrature_degree(::Simplex{D}) where {D} = maximum(p for (d, p) in keys(_simplex_rules) if d == D)
   ```
   Use `joinpath(@__DIR__, "simplex_quadrature_rules.txt")` to locate the
   file. The parser is a plain loop over `eachline` with `split` and
   `parse(Float64, ...)`.
3. Tests. Replace the quadrature part of the "Polynomials" testset:
   - For `D in 2:3` and every available `p`: integrate every monomial
     `ξ^α` with `sum(α) <= p` against the exact value
     `prod(factorial.(big.(α))) / factorial(big(sum(α) + D))` with `atol = 1e-11`
     (the tabulated rules are accurate to about 1e-12 in absolute terms).
   - For `Block{D}`, `D in 1:3`, `p in 1:12`: same monomial test with exact
     value `prod(1 ./ (α .+ 1))`, and `length(w) == cld(p+1, 2)^D`.
   - `eltype(quadrature(Simplex{2}(), 4, Float32)[1]) == Float32`.
   - `@test_throws ErrorException quadrature(Simplex{2}(), 1000)`.
   - `Base.return_types(quadrature, (Simplex{2}, Int))[1] == Tuple{Matrix{Float64}, Vector{Float64}}`.

### Step 7. Makie extension on the recipe system, precompile workload

Files: `ext/HighOrderMeshesMakieExt.jl`, `test/plotting/runtests.jl`,
`test/plotting/Project.toml`.

1. Define two recipes with `Makie.@recipe`:
   - `MeshPlot` with argument `msh` and attributes `labels=()`,
     `reltol=1e-3`, `abstol=Inf`, `maxref=6`,
     `colors=("#CCFFCC", :black, :blue, :darkgray, :darkblue)`. Its `plot!`
     method calls `viz_mesh` and adds `poly!`, `lines!`, and optionally
     `scatter!` and `text!` exactly as the current code does.
   - `SolutionPlot` with arguments `msh, u` and attributes `nsub=nothing`,
     `mesh_edges=false`. Its `plot!` method calls `viz_solution` and adds
     `lines!` in 1D or `mesh!` in 2D, plus the edge overlay.
   Set `Makie.plottype(::HighOrderMesh) = MeshPlot` and
   `Makie.plottype(::HighOrderMesh, ::AbstractArray) = SolutionPlot` so that
   `plot(msh)`, `plot!(msh)`, `plot(msh, u)` and `plot!(msh, u)` all work.
2. Keep the current user-facing behaviour that `plot(msh)` and `plot(msh, u)`
   return a `Figure` with `DataAspect()` in 2D and a `Colorbar` for
   solutions. Implement this with `Makie.plot(m::HighOrderMesh; kw...)` and
   `Makie.plot(m::HighOrderMesh, u::AbstractArray; kw...)` methods that
   create `Figure()` and `Axis(f[1,1]; aspect)` and call the recipe's bang
   form, then add the colorbar with `limits = extrema(u)`. Remove
   `init_figure` and the axis-reuse behaviour.
3. Add a precompile workload at the end of the extension module, guarded so
   that it only runs while precompiling and without any new dependency:
   ```julia
   if ccall(:jl_generating_output, Cint, ()) == 1
       let m = ex1mesh(nref=1), u = ex1solution(m)
           plot(m); plot(m, u); plot(set_degree(mshline(3), 2), zeros(3, 3))
       end
   end
   ```
   Wrap the workload in `try ... catch end` so that a failure can never
   break loading.
4. Remove the Plots part from `test/plotting/runtests.jl` and its
   `Project.toml` if any remains, and add checks that `plot!(m)` on an
   existing axis adds plots and that `plot(m, u)` produces a `Colorbar`.
   Run the plotting tests locally with the command in that file's header.

### Step 8. Gmsh physical names, export cleanup, sample meshes into src

Files: `src/io/converters.jl`, `src/HighOrderMeshes.jl`,
`examples/sample_meshes.jl` (moved to `src/mesh/sample_meshes.jl`), `README.md`.

1. Add `gmsh_physical_names(fname) -> Dict{Int,String}` that parses only the
   `$PhysicalNames` section (reuse `parse_gmsh`) and maps the physical tag
   of every `(D-1)`-dimensional group to its name, with quotes stripped.
2. Give `gmsh2msh` a keyword `verbose=true`; when true it prints one line
   per boundary tag, `boundary 1: "Circle"`, using the names when present.
3. Move `examples/sample_meshes.jl` to `src/mesh/sample_meshes.jl`, include
   it from the main module after `basic_meshes.jl`, and delete
   `ex_vtkexport`, `ex_plot_mesh`, `ex_plot_solution`, `ex_read_gmsh` (they
   referenced `plot` from inside the package, which never worked). Keep
   `ex1mesh`, `ex1solution`, `gmsh_sphere`, and document each with a
   copy-paste example.
4. Exports: remove `name`, `el2nb`, `equispaced`, `blockmesh_hypercube` from
   the export list. Keep everything else that a user calls. Check with
   `grep -n "^export" src/HighOrderMeshes.jl` that each exported name is
   defined.
5. Replace `if fid == -1` in `parse_gmsh` with nothing; `open` throws.
6. Tests: `gmsh_physical_names(joinpath(pkgdir(HighOrderMeshes), "examples/gmsh/circle_tris.msh")) == Dict(1 => "Circle")`
   and `gmsh2msh(...; verbose=false)` still returns the same mesh.

### Step 9. Documenter skeleton and docs workflow

Files: `docs/Project.toml`, `docs/make.jl`, `docs/src/*.md`,
`.github/workflows/Docs.yml`, `README.md`, `.gitignore`.

1. `docs/Project.toml` with `Documenter` and `HighOrderMeshes` (added with
   `Pkg.develop(path="..")` when building). Documenter is a docs-only
   dependency and does not enter the package `Project.toml`.
2. `docs/make.jl` with `makedocs(sitename="HighOrderMeshes.jl", modules=[HighOrderMeshes], pages=[...])`
   and `deploydocs(repo="github.com/popersson/HighOrderMeshes.jl.git", push_preview=true)`.
3. Pages, each a placeholder with a heading, two sentences, and where
   indicated a runnable snippet:
   - `index.md`: what the package is, installation, the three-line hello
     world.
   - `quickstart.md`: build a mesh, raise the degree, plot, save.
   - `meshformat.md`: `x`, `el`, `nb`, the node order, the face convention
     with the tensor-order figure as a text diagram, `dg_nodes`.
   - `basicmeshes.md`: gallery with one snippet each for `mshsquare`,
     `mshcircle`, `mshcube`, `ex1mesh`, `gmsh_sphere`, and periodic meshes.
   - `gmsh.md`: how to write a `.geo`, run `rungmsh2msh`, boundary tags.
   - `io.md`: `.hom` layout copied from the `io.jl` header, VTK export.
   - `api.md`: `@autodocs` of the module.
4. `Docs.yml`: run on push to `main` and on pull requests, Julia 1, with
   `permissions: contents: write, pull-requests: read, statuses: write`,
   using `julia-actions/julia-docdeploy@v1` and `GITHUB_TOKEN`.
5. Add `docs/build/` to `.gitignore`. Update `README.md`: remove the Plots
   section, mention the docs link, keep the quick start.
6. Verify locally with `julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate(); include("docs/make.jl")'`
   and fix every doctest or missing docstring warning.

### Step 10. Remaining tests and small cleanups

Files: `test/runtests.jl`, `src/mesh/mesh_utils.jl`, `src/mesh/high_order_mesh.jl`.

1. Gmsh affine round trip, guarded by `Sys.which("gmsh") !== nothing`: for
   tri, quad, tet and hex meshes at `p = 1:5` generated with `gmshstr2msh`
   from straight-sided geometries, check that the high-order nodes equal the
   `p = 1` interpolation of the corner nodes to `1e-10`. Use the geometry
   strings from `docs/dev/gmsh_affine_test.geo.jl` if present, otherwise a
   unit square (tris; and quads via `Transfinite Surface` and `Recombine
   Surface`) and a unit cube (tets via `Box`; hexes via a transfinite
   extrusion with `Layers` and `Recombine`).
2. `unique_mesh_nodes(x, el; tol=sqrt(eps(T)))` keyword, threaded through
   `snap`, documented as relative to the largest coordinate magnitude.
3. `Base.show` of `HighOrderMesh` says "nodes" instead of "vertices".
4. Test `set_ref_nodes` on a Block mesh with a non-symmetric 1D line throws
   `ArgumentError`, and `FiniteElement(Simplex{2}(), nodes)` with a node
   moved off its symmetric position throws `ArgumentError`.
5. Test that `set_degree(m, 4)` followed by `set_degree(m, 2)` reproduces
   `set_degree(m, 2)` to `1e-12` for `m = mshsquare(2)`, and the same for a
   triangle mesh from `ex1mesh(eg=Simplex{2}(), nref=1)` at degree 3.

### Step 11. Version bump, changelog, final review

Fable: bump `Project.toml` to 0.3.0, write `CHANGELOG.md` with the breaking
changes and the new API, review the whole branch, merge into `main`.
