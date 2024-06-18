const meshgreen = Plots.RGB(0.8, 1, 0.8)
function Plots.plot(m::HighOrderMesh{2,G,P,T};
    labels=(), reltol=1e-3, abstol=Inf, maxref=6,
    colors=(meshgreen, :black, :blue, :darkgray, :darkblue)
    # colors: elements, int edges, bnd edges, ho-nodes, vertices
) where {G,P,T}

    F = Float64
    f2n = mkface2nodes(m.fe)
    xdg = dg_nodes(m)
    order = plot_face_order(G())
    nf, nel = size(m.nbor)
    curves = fill(zeros(F,0,0), nf, nel)
    for iel = 1:nel, iface = 1:nf
        jel, jface = m.nbor[iface, iel]
        if iel < jel || jel < 1
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
            if jel > 0
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

    h = Plots.plot(aspect_ratio=:equal, leg=false)
    Plots.plot!(elem_lines[:, 1], elem_lines[:, 2], seriestype=:shape,
        linecolor=nothing, fillcolor=colors[1])
    Plots.plot!(int_lines[:, 1], int_lines[:, 2], linewidth=1, color=colors[2])
    Plots.plot!(bnd_lines[:, 1], bnd_lines[:, 2], linewidth=2, color=colors[3])

    ### Labels

    if isa(labels, Symbol)
        labels = (labels,)
    end

    if :nodes in labels
        scatter!(xdg[:,:,1], xdg[:,:,2], marker=(:circle, 2, 1.0, colors[4]))
    end

    if :elements in labels
        mid = midpoint(G())
        if :nodes in labels && any(all(m.fe.ref_nodes.vol .- mid' .≈ 0, dims=2))
            mid .*= (P - 1) / P
        end
        mids = eval_fcn(m.fe, xdg, mid')[1,:,:]
        for iel = 1:nel
            annotate!(mids[iel,1], mids[iel,2], text("$iel", 8, :center))
        end
    end

    return h
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

function Plots.plot(m::HighOrderMesh{D,G,P,T}, u::Array{T}; nsub=nothing) where {D,G,P,T}
    F = Float64

    if size(u,1) == size(m.x,1)
        # Assume CG solution
        u = u[m.el,:]
    elseif size(u,1) == size(m.el,1)
        # Assume DG solution
    end
    if ndims(u) > 2
        u = u[:,:,1]
    end

    if isnothing(nsub)
        P == 1 ? nsub = 1 : nsub = 3P
    end
    
    x1,el1 = subelement_mesh(G(), nsub, T)
    xdg = dg_nodes(m)

    allx,allu = Matrix{F}[], F[]
    allel1 = Matrix{Int}[]
    for iel = 1:size(u,2)
        xy1 = F.(eval_fcn(m.fe, xdg[:,iel,:], x1))
        u1 = F.(eval_fcn(m.fe, u[:,iel], x1))
        push!(allx, xy1)
        append!(allu, u1[:,1])
        push!(allel1, el1 .+ (iel-1)*size(x1,1))
    end
    allel1 = hcat(allel1...)
    allx = vcat(allx...)

    if D == 1
        return Plots.plot(allx[allel1], allu[allel1], color=:black, legend=false)
    elseif D == 2
        return Plots.tripcolor(allx[:,1], allx[:,2], allu, allel1, color=:viridis, aspect_ratio=:equal)
    else
        throw("Dimension not implemented")
    end
end
