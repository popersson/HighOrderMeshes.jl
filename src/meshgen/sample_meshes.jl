sample_mesh(::Block{2}) = [0 0; 1 0; 0 1; 1.2 0.9],
                          [1,2,3,4][:,:]
sample_mesh(::Simplex{2}) = [0 0; 1 0; 0 1; 1 1; .7 .5],
                            [[1,2,5] [2,4,5] [4,3,5] [3,1,5]]
sample_fcn(x) = sin.(-sum((x.-0.5).^2/0.2^2,dims=ndims(x)))

"""
    ex1mesh(; nref=3, eg=Block{2}())

Creates a sample high-order curved mesh by applying a sinusoidal transformation
to a refined unit hypercube mesh of the specified `ElementGeometry`.

```julia
using HighOrderMeshes
msh = ex1mesh(nref=2, eg=Simplex{2}())
```
"""
function ex1mesh(; nref=3, eg=Block{2}())
    x,el = sample_mesh(eg)
    m = HighOrderMesh(x,el)
    m = uniref(m,nref)
    m = set_degree(m,3)
    x,y = m.x[:,1], m.x[:,2]
    fx = @. 0.1*sin(2pi*x)
    y2 = @. fx + (1-fx)*y
    m.x[:,2] .= y2
    m
end

"""
    ex1solution(m; dg=true)

Generates a sample solution field (Gaussian-like) on the given mesh `m`.
If `dg=true`, the solution is returned on the DG nodes.

```julia
using HighOrderMeshes
msh = ex1mesh()
u = ex1solution(msh)
```
"""
function ex1solution(m; dg=true)
    x = dg ? dg_nodes(m) : m.x
    u = sample_fcn(x)
    dg && (u[:,size(u,2)÷3,:] .+= 0.5)
    u
end

"""
    gmsh_sample(shape=:circle; h=0.25, p=1, eg=Simplex{D}(), verbose=false)

Unstructured mesh of degree `p` and element size about `h` of a simple
domain, generated with gmsh (which must be on the system `PATH`). For curved
domains gmsh places the high-order boundary nodes on the exact geometry.

| `shape`   | Domain               | Boundary regions                                |
|:--------- |:-------------------- |:----------------------------------------------- |
| `:square` | `[0,1]^2`            | `1..4`: `xmin, xmax, ymin, ymax`, as `mshsquare` |
| `:circle` | unit disk            | `1`: circle                                      |
| `:cube`   | `[0,1]^3`            | `1..6`: `xmin, xmax, ..., zmax`, as `mshcube`    |
| `:sphere` | unit ball            | `1`: sphere                                      |

The element geometry `eg` defaults to triangles (2D) or tetrahedra (3D).
`eg=Block{2}()` recombines the triangles into quads (Frontal-Delaunay with
blossom full-quad recombination, which avoids the nearly flat corners at the
boundary of gmsh's quad algorithm on the circle), and `eg=Block{3}()`
splits every tetrahedron into 4 hexahedra (from a tetrahedral mesh of size
`2h`, so the hexahedra have size about `h`). With `verbose`, gmsh prints its
log and the boundary names.

```julia
msh = gmsh_sample(:circle; h=0.2, p=3)                  # curved triangles
msh = gmsh_sample(:circle; h=0.2, p=3, eg=Block{2}())   # curved quads
msh = gmsh_sample(:square; h=0.1)                       # unstructured triangles
msh = gmsh_sample(:sphere; h=0.3, p=2)                  # curved tetrahedra
msh = gmsh_sample(:cube; h=0.5, eg=Block{3}())          # unstructured hexahedra
```
"""
function gmsh_sample(shape::Symbol=:circle; h=0.25, p=1, eg=nothing, verbose=false)
    haskey(_GMSH_SAMPLES, shape) ||
        throw(ArgumentError("unknown shape :$shape; choose from $(keys(_GMSH_SAMPLES))"))
    D, geo = _GMSH_SAMPLES[shape]
    eg = something(eg, Simplex{D}())
    eg isa ElementGeometry{D} ||
        throw(ArgumentError("element geometry $eg does not match the $(D)D shape :$shape"))
    opts = eg isa Simplex ? "" :
           D == 2 ? "Mesh.RecombineAll = 1;\nMesh.Algorithm = 6;\nMesh.RecombinationAlgorithm = 3;" :
                    "Mesh.SubdivisionAlgorithm = 2;"
    hgmsh = eg isa Block{3} ? 2h : h
    str = """
        SetFactory("OpenCASCADE");
        $geo
        Mesh.MeshSizeMin = $hgmsh;
        Mesh.MeshSizeMax = $hgmsh;
        $opts
        """
    gmshstr2msh(str; p, verbose)
end

# Dimension and OpenCASCADE geometry with physical groups of each sample shape
const _GMSH_SAMPLES = Dict(
    :square => (2, """
        Rectangle(1) = {0, 0, 0, 1, 1};
        Physical Curve("xmin", 1) = {4};
        Physical Curve("xmax", 2) = {2};
        Physical Curve("ymin", 3) = {1};
        Physical Curve("ymax", 4) = {3};
        Physical Surface("domain", 1) = {1};"""),
    :circle => (2, """
        Disk(1) = {0, 0, 0, 1};
        Physical Curve("circle", 1) = {1};
        Physical Surface("domain", 1) = {1};"""),
    :cube => (3, """
        Box(1) = {0, 0, 0, 1, 1, 1};
        Physical Surface("xmin", 1) = {1};
        Physical Surface("xmax", 2) = {2};
        Physical Surface("ymin", 3) = {3};
        Physical Surface("ymax", 4) = {4};
        Physical Surface("zmin", 5) = {5};
        Physical Surface("zmax", 6) = {6};
        Physical Volume("domain", 1) = {1};"""),
    :sphere => (3, """
        Sphere(1) = {0, 0, 0, 1};
        Physical Surface("sphere", 1) = {1};
        Physical Volume("domain", 1) = {1};"""),
)
