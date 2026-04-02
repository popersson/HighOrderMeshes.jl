export ex1mesh, ex1solution, ex_vtkexport, ex_plot_mesh, ex_plot_solution, ex_read_gmsh, gmsh_sphere

sample_mesh(::Block{2}) = [0 0; 1 0; 0 1; 1.2 0.9],
                          [1,2,3,4][:,:]
sample_mesh(::Simplex{2}) = [0 0; 1 0; 0 1; 1 1; .7 .5],
                            [[1,2,5] [2,4,5] [4,3,5] [3,1,5]]
sample_fcn(x) = sin.(-sum((x.-0.5).^2/0.2^2,dims=ndims(x)))

"""
    ex1mesh(; nref=3, eg=Block{2}())

Creates a sample high-order curved mesh by applying a sinusoidal transformation 
to a refined unit hypercube mesh of the specified `ElementGeometry`.
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
"""
function ex1solution(m; dg=true)
    x = dg ? dg_nodes(m) : m.x
    u = sample_fcn(x)
    dg && (u[:,size(u,2)÷3,:] .+= 0.5)
    u
end

"""
    ex_vtkexport(fname="ex1.vtk")

Exports a sample high-order mesh and solution to a VTK file.
"""
function ex_vtkexport(fname="ex1.vtk")
    m = ex1mesh()
    u = ex1solution(m)
    vtkwrite(fname, m, u)
end

"""
    ex_plot_mesh()

Plots a sample high-order mesh with node and element labels. 
Requires a visualization backend like `Makie` or `Plots` with `TriplotRecipes`.
"""
function ex_plot_mesh()
    m = ex1mesh()
    plot(m, labels=(:nodes, :elements))
end

"""
    ex_plot_solution()

Plots a sample solution on a high-order mesh.
Requires a visualization backend like `Makie` or `Plots` with `TriplotRecipes`.
"""
function ex_plot_solution()
    m = ex1mesh()
    u = ex1solution(m)
    plot(m, u)
end

"""
    ex_read_gmsh(fname=joinpath(pkgdir(HighOrderMeshes), "examples/gmsh/circle_quads.msh"))

Reads a Gmsh file and plots the resulting mesh.
Requires a visualization backend like `Makie` or `Plots` with `TriplotRecipes`.
"""
function ex_read_gmsh(fname=joinpath(pkgdir(HighOrderMeshes), "examples/gmsh/circle_quads.msh"))
    m = gmsh2msh(fname)
    plot(m, labels=(:nodes, :elements))
end

"""
    gmsh_sphere()

Generates a first-order mesh of a unit sphere using Gmsh's OpenCASCADE kernel.
"""
function gmsh_sphere(; hmax=0.5, porder=1)
    gmsh = """
        SetFactory("OpenCASCADE");
        Sphere(1) = {0, 0, 0, 1};
        Mesh.MeshSizeMax = $hmax;
    """
    m1 = gmshstr2msh(gmsh; porder)
end

