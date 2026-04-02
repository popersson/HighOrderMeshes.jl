module HighOrderMeshesMakieExt

using HighOrderMeshes
using Makie

const meshgreen = "#CCFFCC"

# Reuse the current axis/figure if one exists, otherwise create a new Figure.
function init_figure(; aspect=DataAspect())
    ax = current_axis()
    if isnothing(ax)
        f  = Figure()
        ax = Axis(f[1,1])
    else
        f = current_figure()
        empty!(ax)
    end
    ax.aspect = aspect
    f, ax
end

"""
    plot(m::HighOrderMesh{2}; labels=(), reltol=1e-3, abstol=Inf, maxref=6,
         colors=(meshgreen, :black, :blue, :darkgray, :darkblue))

Plot a 2D mesh. Element interiors, interior edges, and boundary edges are drawn
in `colors[1..3]`. Optional `labels` may include `:nodes` (scatter the global
nodes) and/or `:elements` (annotate each element with its index).
`reltol`, `abstol`, `maxref` are forwarded to `viz_mesh` for curved-edge refinement.
"""
function Makie.plot(m::HighOrderMesh{2,G,P,T};
    labels=(),
    reltol=1e-3, abstol=Inf, maxref=6,
    colors=(meshgreen, :black, :blue, :darkgray, :darkblue)
    # colors: (element fill, interior edges, boundary edges, nodes, vertices)
) where {G,P,T}
    elem_lines, int_lines, bnd_lines, elem_mid =
        viz_mesh(m; reltol, abstol, maxref)

    # Makie's poly! expects a vector of polygons; split on NaN separators
    function split_nans(l)
        ix = findall(isnan.(l[:,1]))
        [ [NTuple{2,Float64}(r) for r in eachrow(l[ix[i]+1:ix[i+1]-1,:])] for i = 1:length(ix)-1 ]
    end

    f, ax = init_figure()
    poly!(ax,   split_nans(elem_lines), color=colors[1])
    lines!(ax,  int_lines, color=colors[2], linewidth=1)
    lines!(ax,  bnd_lines, color=colors[3], linewidth=2)

    labels isa Symbol && (labels = (labels,))
    :nodes    in labels && scatter!(ax, m.x[:,1], m.x[:,2],
                                    color=colors[4], markersize=8,
                                    marker=:circle, strokewidth=1.0)
    :elements in labels && text!(ax, elem_mid[:,1], elem_mid[:,2],
                                  text=string.(1:size(elem_mid,1)),
                                  align=(:center,:center), fontsize=16)
    return f
end

"""
    plot(m::HighOrderMesh, u; nsub=nothing, mesh_edges=false)

Plot the FEM solution `u` on mesh `m`.

- 1D: line plot of `u` vs. `x`.
- 2D: filled triangle plot (color-mapped). Pass `mesh_edges=true` to overlay
  interior and boundary edges in black.

`nsub` controls sub-element refinement (default: `1` for `p=1`, `3p` otherwise).
Note: contour lines are not supported in Makie.
"""
function Makie.plot(m::HighOrderMesh{D,G,P,T}, u::Array{T};
                    nsub=nothing, mesh_edges=false, contours=0) where {D,G,P,T}
    allx, allu, allel = viz_solution(m, u; nsub)

    if D == 1
        f, ax = init_figure(aspect=nothing)
        lines!(ax, allx[allel[:],1], allu[allel[:],1])
    elseif D == 2
        f, ax = init_figure()
        mesh!(ax, allx, allel', color=allu, shading=NoShading)
        contours != 0 && error("Contours not supported in Makie")
        if mesh_edges
            _, int_lines, bnd_lines, _ = viz_mesh(m)
            lines!(ax, vcat(int_lines, bnd_lines), linewidth=1.0, color=:black)
        end
        lims = (minimum(allu), maximum(allu))
        if typeof(f.content[end]) == Colorbar
            f.content[end].limits = lims
        else
            Colorbar(f[1,2], limits=lims)
        end
    else
        error("plot not implemented for D=$D")
    end
    return f
end

end
