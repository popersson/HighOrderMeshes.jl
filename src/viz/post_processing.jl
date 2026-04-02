###########################################################################
## Mesh visualization geometry

"""
    viz_mesh(m::HighOrderMesh{2}; reltol=1e-3, abstol=Inf, maxref=6)

Compute plotting geometry for a 2D mesh. Returns `(elem_lines, int_lines, bnd_lines, elem_mid)`:

- `elem_lines`: `(npts × 2)` matrix of `(x,y)` coordinates tracing all element
  boundaries in order, separated by `NaN` rows. Suitable for a single polyline plot.
- `int_lines`: same format, interior edges only (shared between two elements).
- `bnd_lines`: same format, boundary edges only.
- `elem_mid`:  `(nel × 2)` matrix of element centroid coordinates (for labels).

Curved edges (high-order meshes, `P > 1`) are adaptively subdivided until the
piecewise-linear approximation error is below `reltol * edge_length` or `abstol`,
using up to `2^maxref` sub-intervals.
"""
function viz_mesh(m::HighOrderMesh{2,G,P,T};
                  reltol=1e-3, abstol=Inf, maxref=6) where {G,P,T}
    F      = Float64
    f2n    = mkface2nodes(m.fe)
    xdg    = dg_nodes(m)
    order  = plot_face_order(G())
    nf, nel = size(m.nb)

    # --- Step 1: compute a refined polyline for every face of every element ---
    curves = fill(zeros(F,0,0), nf, nel)
    for iel = 1:nel, iface = 1:nf
        xy = xdg[f2n[:, iface], iel, :]

        evaledge(nsub)    = F.(eval_field(m.fe, xy, (0:nsub) ./ T(nsub)))
        rowwise2norm(X)   = sqrt.(sum(X .^ 2, dims=2))

        # Refine a coarse polyline by inserting midpoints nref times
        function piecewise_linear(x, nref)
            y = copy(x)
            for _ = 1:nref
                w = zeros(size(y,1)*2 - 1, size(y,2))
                w[1:2:end, :]     = y
                w[2:2:end-1, :]   = (y[1:end-1,:] + y[2:end,:]) / 2
                y = w
            end
            y
        end

        X  = evaledge(1)
        L0 = rowwise2norm(X[[end],:] - X[[1],:])[1]  # chord length

        if P > 1
            # Evaluate at finest resolution, then coarsen until tolerance is met
            X = evaledge(2^maxref)
            for nref = 0:maxref-1
                cX    = X[1:2^(maxref-nref):end, :]
                dX    = X - piecewise_linear(cX, maxref - nref)
                dXmax = maximum(rowwise2norm(dX))
                if dXmax < abstol && dXmax < L0 * reltol
                    X = cX
                    break
                end
            end
        end
        curves[iface, iel] = X
    end

    # --- Step 2: assemble NaN-separated polyline arrays for each line category ---
    nan_nan   = F.([NaN NaN])
    elem_lines = [nan_nan]
    int_lines  = [nan_nan]
    bnd_lines  = [nan_nan]

    for iel = 1:nel
        for iface = 1:nf
            # Element outline: concatenate faces in plotting order
            cxy = curves[order[iface], iel]
            iface > 1 && (cxy = cxy[2:end, :])  # drop shared endpoint
            push!(elem_lines, cxy)

            jel, _, _ = m.nb[iface, iel]
            0 < jel < iel && continue  # interior edge already added from other side

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
    int_lines  = vcat(int_lines...)
    bnd_lines  = vcat(bnd_lines...)

    # Element centroids for optional element-number labels
    mid      = midpoint(G())
    elem_mid = eval_field(m.fe, xdg, mid')[1,:,:]

    elem_lines, int_lines, bnd_lines, elem_mid
end

###########################################################################
## Solution visualization helpers

# Build a uniform reference sub-mesh on the element geometry for smooth
# solution plotting. Returns (x, el): reference coords and triangle connectivity.
subelement_mesh(::ElementGeometry, n, T=Float64) = error("subelement_mesh not implemented")

function subelement_mesh(::Block{1}, n, T=Float64)
    x  = collect((0:n)[:,[1]] ./ T(n))
    el = collect(hcat(1:n, 2:n+1)')
    x, el
end

subelement_mesh(::Simplex{1}, n, T=Float64) = subelement_mesh(Block{1}(), n, T)

function subelement_mesh(::Block{2}, n, T=Float64)
    s  = (0:n) ./ T(n)
    sx = s .+ 0*s'
    sy = 0*s .+ s'
    x  = hcat(sx[:], sy[:])
    # Split each quad into two triangles
    el = hcat( ([i,i+1,i+n+1,i+n+2] for i = 1:(n+1)^2-n-1 if mod(i,n+1) != 0)... )
    el = hcat(el[[1,2,3],:], el[[2,4,3],:])
    x, el
end

function subelement_mesh(::Simplex{2}, n, T=Float64)
    x, el  = subelement_mesh(Block{2}(), n, T)
    keep   = (sum(x, dims=2) .<= 1)[:]
    x      = x[keep,:]
    el     = el[:, all(keep[el], dims=1)[:]]
    remap  = similar(keep, Int)
    remap[keep] = 1:sum(keep)
    el = remap[el]
    x, el
end

"""
    mesh_function_type(m::HighOrderMesh, u::Array) -> :cg or :dg

Infer whether `u` is a CG field (one value per global node, `size(u,1) == nnodes`)
or a DG field (one value per local element node, `size(u,1) == nnodes_per_elem`).
Errors if neither matches.
"""
function mesh_function_type(m::HighOrderMesh, u::Array)
    size(u,1) == size(m.x,1)  && return :cg
    size(u,1) == size(m.el,1) && return :dg
    error("Solution field does not match CG or DG node count for the mesh")
end

function convert_3dg_solution(m::HighOrderMesh, u)
    u = permutedims(u, (1,3,2))
    u[node_order_3dg(m),:,:] = u
    u
end

###########################################################################
## Solution visualization

"""
    viz_solution(m::HighOrderMesh, u; nsub=nothing)

Evaluate the FEM solution `u` on a refined sub-mesh for smooth plotting.
Returns `(allx, allu, allel)`:

- `allx`:  `(npts × D)` physical coordinates of all sub-mesh nodes.
- `allu`:  `(npts,)` solution values at those nodes.
- `allel`: `(3 × ntris)` triangle connectivity into `allx` (triangles for both
           simplex and block elements).

`u` may be a CG field (`nnodes` values) or a DG field (`nnodes_per_elem × nel`).
For CG fields, coincident sub-mesh nodes are deduplicated so that `allx` and
`allu` have no repeated entries.

`nsub` controls the number of sub-intervals per element edge (default: `1` for
`p=1`, `3p` for higher order).
"""
function viz_solution(m::HighOrderMesh{D,G,P,T}, u::Array{T}; nsub=nothing) where {D,G,P,T}
    F     = Float64
    is_cg = mesh_function_type(m, u) == :cg

    # Normalize to DG format: nnodes_per_elem × nel
    if is_cg
        u = u[m.el,:]                             # gather to element-local DOFs
    elseif size(u,3) == size(m.el,2)              # 3DG format: reorder axes
        u = convert_3dg_solution(m, u)
    end
    ndims(u) > 2 && (u = u[:,:,1])               # take first component if vector field

    nsub = isnothing(nsub) ? (P == 1 ? 1 : 3P) : nsub

    x1, el1 = subelement_mesh(G(), nsub, T)
    xdg     = dg_nodes(m)

    # Evaluate physical coords and solution on the sub-mesh for each element
    allx  = Matrix{F}[]
    allu  = F[]
    allel = Matrix{Int}[]
    for iel = 1:size(u,2)
        xy1 = F.(eval_field(m.fe, xdg[:,iel,:], x1))
        u1  = F.(eval_field(m.fe, u[:,iel], x1))
        push!(allx, xy1)
        append!(allu, u1[:,1])
        push!(allel, el1 .+ (iel-1)*size(x1,1))
    end
    allel = hcat(allel...)
    allx  = vcat(allx...)

    # For CG fields, merge coincident nodes from adjacent elements
    if is_cg
        allx, allel, ix = unique_mesh_nodes(allx, allel, output_ix=true)
        allu = allu[ix]
    end

    allx, allu, allel
end
