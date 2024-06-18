function blockmesh_hypercube(dims::NTuple{D,Int}, T=Float64) where D
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

mshhypercube(dims::NTuple{D,Int}, T=Float64) where D =
    HighOrderMesh(blockmesh_hypercube(dims,T)...)

mshcube(m=5, n=m, o=n, T=Float64) = mshhypercube((m,n,o), T)
mshsquare(m=5, n=m, T=Float64) = mshhypercube((m,n), T)
mshline(m=5, T=Float64) = mshhypercube((m,), T)
