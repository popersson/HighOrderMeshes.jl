module HighOrderMeshesPlotsExt

using HighOrderMeshes
using Plots, TriplotRecipes

const meshgreen = Plots.RGBX(0.8, 1, 0.8)

"""
    plot(m::HighOrderMesh{2}; labels=(), reltol=1e-3, abstol=Inf, maxref=6,
         colors=(meshgreen, :black, :blue, :darkgray, :darkblue))

Plot a 2D mesh. Element interiors, interior edges, and boundary edges are drawn
in `colors[1..3]`. Optional `labels` may include `:nodes` (scatter the global
nodes) and/or `:elements` (annotate each element with its index).
`reltol`, `abstol`, `maxref` are forwarded to `viz_mesh` for curved-edge refinement.
"""
function Plots.plot(m::HighOrderMesh{2,G,P,T};
    labels=(),
    reltol=1e-3, abstol=Inf, maxref=6,
    colors=(meshgreen, :black, :blue, :darkgray, :darkblue)
    # colors: (element fill, interior edges, boundary edges, nodes, vertices)
) where {G,P,T}
    elem_lines, int_lines, bnd_lines, elem_mid =
        viz_mesh(m; reltol, abstol, maxref)

    h = Plots.plot(aspect_ratio=:equal, leg=false)
    Plots.plot!(elem_lines[:,1], elem_lines[:,2], seriestype=:shape,
                linecolor=nothing, fillcolor=colors[1])
    Plots.plot!(int_lines[:,1], int_lines[:,2], linewidth=1, color=colors[2])
    Plots.plot!(bnd_lines[:,1], bnd_lines[:,2], linewidth=2, color=colors[3])

    labels isa Symbol && (labels = (labels,))
    :nodes    in labels && scatter!(m.x[:,1], m.x[:,2],
                                    marker=(:circle, 2, 1.0, colors[4]))
    :elements in labels && for (iel, mid) in enumerate(eachrow(elem_mid))
                               annotate!(mid[1], mid[2], text("$iel", 8, :center))
                           end
    return h
end

"""
    plot(m::HighOrderMesh, u; nsub=nothing, mesh_edges=false, contours=0)

Plot the FEM solution `u` on mesh `m`.

- 1D: line plot of `u` vs. `x`.
- 2D: filled triangle plot (color-mapped via `tripcolor`). Pass `contours=n` to
  overlay `n` contour lines (CG fields only). Pass `mesh_edges=true` to overlay
  interior and boundary edges in black.

`nsub` controls sub-element refinement (default: `1` for `p=1`, `3p` otherwise).
"""
function Plots.plot(m::HighOrderMesh{D,G,P,T}, u::Array{T};
                    nsub=nothing, mesh_edges=false, contours=0) where {D,G,P,T}
    allx, allu, allel = viz_solution(m, u; nsub)

    if D == 1
        return Plots.plot(allx[allel], allu[allel], color=:black, legend=false)
    elseif D == 2
        h = TriplotRecipes.tripcolor(allx[:,1], allx[:,2], allu, allel,
                                     color=:viridis, aspect_ratio=:equal)
        if contours != 0
            mesh_function_type(m, u) != :cg && error("Contours only supported for CG fields")
            h = TriplotRecipes.tricontour!(allx[:,1], allx[:,2], allu, allel,
                                           contours, colorbar_entry=false)
        end
        if mesh_edges
            _, int_lines, bnd_lines, _ = viz_mesh(m)
            Plots.plot!(vcat(int_lines, bnd_lines), linewidth=0.5, color=:black, legend=false)
        end
        return h
    else
        error("plot not implemented for D=$D")
    end
end

end
