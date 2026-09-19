export ex1mesh, ex1solution, gmsh_sphere

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
    gmsh_sphere(; hmax=0.5, porder=1)

Generates a first-order mesh of a unit sphere using Gmsh's OpenCASCADE kernel.
Requires `gmsh` to be on the system `PATH`.

```julia
using HighOrderMeshes
msh = gmsh_sphere(hmax=0.3, porder=2)
```
"""
function gmsh_sphere(; hmax=0.5, porder=1)
    gmsh = """
        SetFactory("OpenCASCADE");
        Sphere(1) = {0, 0, 0, 1};
        Mesh.MeshSizeMax = $hmax;
    """
    m1 = gmshstr2msh(gmsh; porder)
end
