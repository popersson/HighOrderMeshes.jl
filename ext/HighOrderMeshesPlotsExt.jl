module HighOrderMeshesPlotsExt

using HighOrderMeshes
using Plots, TriplotRecipes

const meshgreen = Plots.RGBX(0.8, 1, 0.8)

"""
    Plots.plot(m::HighOrderMesh{2,G,P,T};
               labels=(), reltol=1e-3, abstol=Inf, maxref=6,
               colors=(meshgreen, :black, :blue, :darkgray, :darkblue)
               # colors: elements, int edges, bnd edges, ho-nodes, vertices)

TBW
"""
function Plots.plot(m::HighOrderMesh{2,G,P,T};
    labels=(), reltol=1e-3, abstol=Inf, maxref=6,
    colors=(meshgreen, :black, :blue, :darkgray, :darkblue)
    # colors: elements, int edges, bnd edges, ho-nodes, vertices
) where {G,P,T}

    elem_lines, int_lines, bnd_lines, elem_mid =
        viz_mesh(m, reltol=reltol, abstol=abstol, maxref=maxref)

    h = Plots.plot(aspect_ratio=:equal, leg=false)
    Plots.plot!(elem_lines[:, 1], elem_lines[:, 2], seriestype=:shape,
                linecolor=nothing, fillcolor=colors[1])
    Plots.plot!(int_lines[:, 1], int_lines[:, 2], linewidth=1, color=colors[2])
    Plots.plot!(bnd_lines[:, 1], bnd_lines[:, 2], linewidth=2, color=colors[3])

    if isa(labels, Symbol)
        labels = (labels,)
    end

    if :nodes in labels
        scatter!(m.x[:,1], m.x[:,2], marker=(:circle, 2, 1.0, colors[4]))
    end

    if :elements in labels
        for (iel,mid) in enumerate(eachrow(elem_mid))
            annotate!(mid[1], mid[2], text("$iel", 8, :center))
        end
    end

    return h
end

"""
    Plots.plot(m::HighOrderMesh{D,G,P,T}, u::Array{T}; nsub=nothing, mesh_edges=false)

TBW
"""
function Plots.plot(m::HighOrderMesh{D,G,P,T}, u::Array{T}; nsub=nothing, mesh_edges=false, contours=0) where {D,G,P,T}
    allx, allu, allel = viz_solution(m, u, nsub=nsub)
    
    if D == 1
        return Plots.plot(allx[allel], allu[allel], color=:black, legend=false)
    elseif D == 2
        h = TriplotRecipes.tripcolor(allx[:,1], allx[:,2], allu, allel, color=:viridis, aspect_ratio=:equal)
        if contours != 0
            if mesh_function_type(m,u) != :cg
                throw("Contours only supported for CG functions")
            end
            h = TriplotRecipes.tricontour!(allx[:,1], allx[:,2], allu, allel,
                                           contours, colorbar_entry=false)
        end
        if mesh_edges
            _, int_lines, bnd_lines, _ =viz_mesh(m)
            all_lines = vcat(int_lines, bnd_lines)
            Plots.plot!(all_lines[:, 1], all_lines[:, 2], linewidth=0.5, color=:black, legend=false)
        end
        return h
    else
        throw("Dimension not implemented")
    end
end

end
