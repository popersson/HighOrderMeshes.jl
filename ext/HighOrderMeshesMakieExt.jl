module HighOrderMeshesMakieExt

using HighOrderMeshes
using Makie

const meshgreen = "#CCFFCC"

function init_figure(; aspect=DataAspect())
    ax = current_axis()
    if isnothing(ax)
        f = Figure()
        ax = Axis(f[1,1])
    else
        f = current_figure()
        empty!(ax)
    end
    ax.aspect = aspect
    f,ax
end


"""
    Makie.plot(m::HighOrderMesh{2,G,P,T};
               labels=(), reltol=1e-3, abstol=Inf, maxref=6,
               colors=(meshgreen, :black, :blue, :darkgray, :darkblue)
               # colors: elements, int edges, bnd edges, ho-nodes, vertices)

TBW
"""
function Makie.plot(m::HighOrderMesh{2,G,P,T};
    labels=(), reltol=1e-3, abstol=Inf, maxref=6,
    colors=(meshgreen, :black, :blue, :darkgray, :darkblue)
    # colors: elements, int edges, bnd edges, ho-nodes, vertices
) where {G,P,T}

    elem_lines, int_lines, bnd_lines, elem_mid =
        viz_mesh(m, reltol=reltol, abstol=abstol, maxref=maxref)

    function split_nans(l)
        ix = findall(isnan.(l[:,1]))
        el = [l[ix[i]+1:ix[i+1]-1,:] for i = 1:length(ix)-1]
        ps = [ [NTuple{2,Float64}(ee) for ee in eachrow(e)] for e in el]
    end
    elem_lines = split_nans(elem_lines)
   
    f,ax = init_figure()
    poly!(ax, elem_lines, color=colors[1])
    lines!(ax, int_lines, color=colors[2], linewidth=1)
    lines!(ax, bnd_lines, color=colors[3], linewidth=2)
    
    if isa(labels, Symbol)
        labels = (labels,)
    end

    if :nodes in labels
        scatter!(ax, m.x[:,1], m.x[:,2], color=colors[4], markersize=8, marker=:circle, strokewidth=1.0)
    end

    if :elements in labels
        text!(ax, elem_mid[:,1], elem_mid[:,2], text=string.(1:size(elem_mid,1)), align=(:center,:center), fontsize=16)
    end

    return f
end

"""
    Makie.plot(m::HighOrderMesh{D,G,P,T}, u::Array{T}; nsub=nothing, mesh_edges=false)

TBW
"""
function Makie.plot(m::HighOrderMesh{D,G,P,T}, u::Array{T}; nsub=nothing, mesh_edges=false, contours=0) where {D,G,P,T}
    allx, allu, allel = viz_solution(m, u, nsub=nsub)
  
    if D == 1
        f,ax = init_figure(aspect=nothing)
        lines!(ax, allx[allel[:],1], allu[allel[:],1])
    elseif D == 2
        f,ax = init_figure()
        mesh!(ax, allx, allel', color=allu, shading=NoShading)
        if contours != 0
            error("Contours not supported in Makie")
        end
        if mesh_edges
            _, int_lines, bnd_lines, _ =viz_mesh(m)
            all_lines = vcat(int_lines, bnd_lines)
            lines!(ax, all_lines[:, 1], all_lines[:, 2], linewidth=1.0, color=:black)
        end
        lims = (minimum(allu), maximum(allu))
        if typeof(f.content[end]) == Colorbar
            f.content[end].limits = lims
        else
            Colorbar(f[1,2], limits=lims)
        end
    else
        error("Dimension not implemented")
    end

    return f
end

end
