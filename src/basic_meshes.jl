"""
    blockmesh_hypercube(dims::NTuple{D,Int}, T=Float64) where D

TBW
"""
function blockmesh_hypercube(dims::NTuple{D,Int}; T=Float64) where D
    x = zeros(T, 1, 0)
    for d in dims
        xx = (0:d) / T(d)
        x = hcat(repeat(x,length(xx),1),
                 reshape(repeat(xx',size(x,1),1),:,1))
    end

    el = [1]
    off = cumprod((0,dims...) .+ 1)
    for i = 1:D
        el2 = repeat(vcat(el, el .+ off[i]),1,1,dims[i]) .+ 
              off[i]*reshape(0:(dims[i]-1),1,1,dims[i])
        el2 = reshape(el2, 2*size(el,1), size(el,2)*dims[i])
        el = el2
    end

    x,el
end

"""
    mshhypercube(dims::NTuple{D,Int}, T=Float64) where D

TBW
"""
function mshhypercube(dims::NTuple{D,Int}; p=1, T=Float64) where D
    bndexpr(x) = [x'; 1 .- x'][:]
    change_degree(HighOrderMesh(blockmesh_hypercube(dims; T=T)..., bndexpr=bndexpr), p)
end

mshcube(m=5, n=m, o=n; kwargs...) = mshhypercube((m,n,o); kwargs...)
mshsquare(m=5, n=m; kwargs...) = mshhypercube((m,n); kwargs...)
mshline(m=5; kwargs...) = mshhypercube((m,); kwargs...)

function mshcircle(n=1; p=1, shape=:full)
    n1 = n * p
    s = (0:n1) ./ n1
    e = ones(n1+1)
    z = zeros(n1+1)
    a = (2 + π/8) / (4 + √2) # All coarsest quads same average (curved) side length

    ix = reshape(1:(n1+1)^2, n1+1, n1+1)
    q0 = fill(0, p+1, p+1, n, n)
    for j = 1:n, i = 1:n
        ix0 = 0:p
        i0,j0 = (i-1)*p + 1, (j-1)*p + 1
        q0[:,:,i,j] .= ix[i0 .+ ix0, j0 .+ ix0]
    end
    q0 = reshape(q0, (p+1)^2, n^2)

    # Boundaries of block
    phi = π/4 * (0:n1)/n1
    x1 = a*e
    y1 = a*s
    x2 = cos.(phi)
    y2 = sin.(phi)
    x3 = range(a, stop=1, length=n1+1)
    y3 = z
    x4 = range(a, stop=1/√2, length=n1+1)
    y4 = x4

    # Transfinite interpolation
    X1 = @. (1-s)*x1' + s*x2'
    Y1 = @. (1-s)*y1' + s*y2'
    X2 = @. x3*(1-s') + x4*s'
    Y2 = @. y3*(1-s') + y4*s'
    X12 = @. (1-s)*(1-s')*X1[1,1]+s*(1-s')*X1[end,1]+(1-s)*s'*X1[1,end]+s*s'*X1[end,end]
    Y12 = @. (1-s)*(1-s')*Y1[1,1]+s*(1-s')*Y1[end,1]+(1-s)*s'*Y1[1,end]+s*s'*Y1[end,end]
    X = @. X1 + X2 - X12
    Y = @. Y1 + Y2 - Y12
    p1 = [X[:] Y[:]]
    X2,Y2 = X[:,end:-1:1], Y[:,end:-1:1]
    p2 = [Y2[:] X2[:]]

    # Middle box
    a1 = a*s
    xx = a*(s.*e')
    yy = a*(e.*s')
    p0 = [xx[:] yy[:]]

    # All three quads
    pp = [p1; p2; p0]
    el = [q0 q0 .+ (n1+1)^2 q0 .+ 2(n1+1)^2]
    pp,el = unique_mesh_nodes(pp,el)

    fe = FiniteElement(Block{2}(), p)
    if shape == :quarter
        return HighOrderMesh(fe, pp, el, bndexpr=x->[x[2],x[1],0])
    else
        # Rotate 90 degrees
        p1 = pp
        p2 = [-pp[:,2] pp[:,1]]
        el = hcat(el, el .+ size(p1,1))
        pp = vcat(p1,p2)
        pp,el = unique_mesh_nodes(pp,el)
        if shape == :half
            return HighOrderMesh(fe, pp, el, bndexpr=x->[x[2],0])
        else
            # Rotate 180 degrees
            p1 = pp
            p2 = [-pp[:,1] -pp[:,2]]
            el = hcat(el, el .+ size(p1,1))
            pp = vcat(p1,p2)
            pp,el = unique_mesh_nodes(pp,el)
            if shape == :full
                return HighOrderMesh(fe, pp, el, bndexpr=p->[0])
            else
                throw("`shape` must be :quarter, :half, or :full")
            end
        end
    end
end
