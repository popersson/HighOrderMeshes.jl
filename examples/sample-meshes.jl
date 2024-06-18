export ex1mesh, ex1solution

sample_mesh(::Block{2}) = [0 0; 1 0; 0 1; 1.2 0.9],
                          [1,2,3,4][:,:]
sample_mesh(::Simplex{2}) = [0 0; 1 0; 0 1; 1 1; .7 .5],
                            [[1,2,5] [2,4,5] [4,3,5] [3,1,5]]
sample_fcn(x) = sin.(-sum((x.-0.5).^2/0.2^2,dims=ndims(x)))

function ex1mesh(nref=3, eg=Block{2}())
    x,el = sample_mesh(eg)
    m = HighOrderMesh(x,el)
    m = uniref(m,nref)
    m = change_degree(m,3)
    x,y = m.x[:,1], m.x[:,2]
    fx = @. 0.1*sin(2pi*x)
    y2 = @. fx + (1-fx)*y
    m.x[:,2] .= y2
    m
end

function ex1solution(m; dg=true)
    x = dg ? dg_nodes(m) : m.x
    u = sample_fcn(x)
    dg && (u[:,size(u,2)÷3,:] .+= 0.5)
    u
end

# Export to VTK:
# m = ex1mesh()
# u = ex1solution(m)
# vtkwrite("ex1.vtk", m, u)

# Plotting:
# m = ex1mesh()
# u = ex1solution(m)
# plot(m, labels=(:nodes, :elements)) # Plot mesh
# plot(m,u)  # Plot solution

# Read gmsh:
# m = gmsh2msh(joinpath(pkgdir(HighOrderMeshes), "examples/gmsh/circle_quads.msh"))
# plot(m, labels=(:nodes, :elements)) # Plot mesh
