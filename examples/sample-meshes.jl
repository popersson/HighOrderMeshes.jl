function test_tri()
    x = [0 0; 1 0; 0 1; 1 1; .7 .5]
    el = [[1,2,5] [2,4,5] [4,3,5] [3,1,5]]
    m = HighOrderMesh(x,el)
end

function test_quad()
    x = [0 0; 1 0; 0 1; 1.2 0.9]
    el = [1,2,3,4][:,:]
    m = HighOrderMesh(x,el)
end

function test_tri3(nref=1)
    m = test_tri()
    m = uniref(m,nref)
    m = change_degree(m,3)
    x,y = m.x[:,1], m.x[:,2]
    fx = @. 0.1*sin(2pi*x)
    y2 = @. fx + (1-fx)*y
    m.x[:,2] .= y2
    m
end

function test_quad3(nref=1)
    m = test_quad()
    m = uniref(m,nref)
    m = change_degree(m,3)
    x,y = m.x[:,1], m.x[:,2]
    fx = @. 0.1*sin(2pi*x)
    y2 = @. fx + (1-fx)*y
    m.x[:,2] .= y2
    m
end

testfcn(x) = sin.(-sum((x.-0.5).^2/0.2^2,dims=ndims(x)))

function ex1(; dim=2, dg=true, nref=3, vtk=false)
    m = dim == 1 ? mshline(2^nref, porder=3) : test_quad3(nref)
    x = dg ? dg_nodes(m) : m.x
    u = testfcn(x)
    dg && (u[:,size(u,2)÷3,:] .+= 0.5)
    vtk && vtkwrite("ex1.vtk", m, u)
    plot(m,u)
end

function ex2()
    m = gmsh2msh("gmsh/circle_quads.msh")
    m = change_degree(m, 5)
    plot(m, labels=(:nodes, :elements))
end

function ex3()
    m = gmsh2msh("gmsh/circle_quads.msh")
    m = change_degree(m, 5)
    plot(m, sin.(2pi*(m.x[:,1] + m.x[:,2])))
end
