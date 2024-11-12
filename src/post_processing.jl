"""
    viz_mesh(m::HighOrderMesh{2,G,P,T};
                  reltol=1e-3, abstol=Inf, maxref=6) where {G,P,T}

TBW
"""
function viz_mesh(m::HighOrderMesh{2,G,P,T};
                  reltol=1e-3, abstol=Inf, maxref=6) where {G,P,T}
    F = Float64
    f2n = mkface2nodes(m.fe)
    xdg = dg_nodes(m)
    order = plot_face_order(G())
    nf, nel = size(m.nbor)
    curves = fill(zeros(F,0,0), nf, nel)
    for iel = 1:nel, iface = 1:nf
        jel, jface = m.nbor[iface, iel]

        is_main_side = true # iel < jel || jel < 1  # removed to allow for periodic boundaries
        if is_main_side
            xy = xdg[f2n[:, iface], iel, :]
            evaledge(nsub) = F.(eval_fcn(m.fe, xy, (0:nsub) ./ T(nsub)))
            rowwise2norm(X) = sqrt.(sum(X .^ 2, dims=2))
            function piecewise_linear(x, nref)
                y = copy(x)
                for iref = 1:nref
                    w = zeros(size(y, 1) * 2 - 1, size(y, 2))
                    w[1:2:end, :] = y
                    w[2:2:end-1, :] = (y[1:end-1, :] + y[2:end, :]) / 2
                    y = w
                end
                y
            end

            # If curved, find subdivision based on reltol/abstol
            X = evaledge(1)
            L0 = rowwise2norm(X[[end], :] - X[[1], :])[1]
            if P > 1
                X = evaledge(2^maxref)
                for nref = 0:maxref-1
                    # Find error in straight line approximation
                    cX = X[1:2^(maxref-nref):end, :]
                    dX = X - piecewise_linear(cX, maxref - nref)
                    dXmax = maximum(rowwise2norm(dX))
                    if dXmax < abstol && dXmax < L0 * reltol
                        X = cX
                        break
                    end
                end
            end

            curves[iface, iel] = X
            if !is_main_side
                curves[jface, jel] = reverse(X, dims=1)
            end
        end
    end

    nan_nan = F.([NaN NaN])
    elem_lines = [nan_nan]
    int_lines = [nan_nan]
    bnd_lines = [nan_nan]
    for iel = 1:nel
        for iface = 1:nf
            cxy = curves[order[iface], iel]
            iface > 1 && (cxy = cxy[2:end, :])
            push!(elem_lines, cxy)

            jel, jface = m.nbor[iface, iel]
            0 < jel < iel && continue

            cxy = curves[iface, iel]
            if iel < jel
                push!(int_lines, cxy, nan_nan)
            elseif jel < 1
                push!(bnd_lines, cxy, nan_nan)
            end
        end
        push!(elem_lines, nan_nan)
    end
    elem_lines = vcat(elem_lines...)
    int_lines = vcat(int_lines...)
    bnd_lines = vcat(bnd_lines...)

    mid = midpoint(G())
    elem_mid = eval_fcn(m.fe, xdg, mid')[1,:,:]
    
    elem_lines, int_lines, bnd_lines, elem_mid
end

subelement_mesh(::ElementGeometry) = throw("Not implemented")

function subelement_mesh(::Block{1}, n, T=Float64)
    x = collect((0:n)[:,[1]] ./ T(n))
    el = collect(hcat(1:n, 2:n+1)')
    x,el
end

subelement_mesh(::Simplex{1}, n, T=Float64) = subelement_mesh(Block{1}(), n, T)

function subelement_mesh(::Block{2}, n, T=Float64)
    s = (0:n) ./ T(n)
    sx = s .+ 0*s'
    sy = 0*s .+ s'
    x = hcat(sx[:], sy[:])
    el = hcat( ([i,i+1,i+n+1,i+n+2] for i = 1:(n+1)^2 - n-1 if mod(i,n+1)!=0 )... )
    el = hcat(el[[1,2,3],:], el[[2,4,3],:])
    x,el
end

function subelement_mesh(::Simplex{2}, n, T=Float64)
    x,el = subelement_mesh(Block{2}(), n, T)
    keep = (sum(x, dims=2) .<= 1)[:]
    x = x[keep,:]
    el = el[:,all(keep[el], dims=1)[:]]
    map = similar(keep, Int)
    map[keep] = 1:sum(keep)
    el = map[el]
    x,el
end

function mesh_function_type(m::HighOrderMesh, u::Array)
    if size(u,1) == size(m.x,1)
        return :cg
    elseif size(u,1) == size(m.el,1)
        return :dg
    else
        throw("Solution field does not match CG or DG type for the mesh")
    end
end

"""
    viz_solution(m::HighOrderMesh{D,G,P,T}, u::Array{T}; nsub=nothing) where {D,G,P,T}

TBW
"""
function viz_solution(m::HighOrderMesh{D,G,P,T}, u::Array{T}; nsub=nothing) where {D,G,P,T}
    F = Float64

    is_cg = mesh_function_type(m,u) == :cg
    is_cg && (u = u[m.el,:])
    
    if ndims(u) > 2
        u = u[:,:,1]
    end

    if isnothing(nsub)
        P == 1 ? nsub = 1 : nsub = 3P
    end
    
    x1,el1 = subelement_mesh(G(), nsub, T)
    xdg = dg_nodes(m)

    allx,allu = Matrix{F}[], F[]
    allel = Matrix{Int}[]
    for iel = 1:size(u,2)
        xy1 = F.(eval_fcn(m.fe, xdg[:,iel,:], x1))
        u1 = F.(eval_fcn(m.fe, u[:,iel], x1))
        push!(allx, xy1)
        append!(allu, u1[:,1])
        push!(allel, el1 .+ (iel-1)*size(x1,1))
    end
    allel = hcat(allel...)
    allx = vcat(allx...)

    if is_cg
        allx,allel,ix = unique_mesh_nodes(allx,allel,output_ix=true)
        allu = allu[ix]
    end
    
    allx, allu, allel
end
