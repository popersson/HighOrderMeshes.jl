module HighOrderMeshesMakieExt

using HighOrderMeshes
using Makie

const meshgreen = "#CCFFCC"

"""
    MeshPlot(msh; labels=(), reltol=1e-3, abstol=Inf, maxref=6,
             colors=(meshgreen, :black, :blue, :darkgray, :darkblue))

Recipe drawing a 2D mesh. Element interiors, interior edges, and boundary
edges are drawn in `colors[1..3]`. Optional `labels` may include `:nodes`
(scatter the global nodes) and/or `:elements` (annotate each element with
its index). `reltol`, `abstol`, `maxref` are forwarded to `viz_mesh` for
curved-edge refinement.
"""
@recipe(MeshPlot, msh) do scene
    Attributes(
        labels = (),
        reltol = 1e-3,
        abstol = Inf,
        maxref = 6,
        colors = (meshgreen, :black, :blue, :darkgray, :darkblue),
        # colors: (element fill, interior edges, boundary edges, nodes, vertices)
    )
end

function Makie.plot!(mp::MeshPlot)
    m = mp.msh[]
    elem_lines, int_lines, bnd_lines, elem_mid =
        viz_mesh(m; reltol=mp.reltol[], abstol=mp.abstol[], maxref=mp.maxref[])

    # Makie's poly! expects a vector of polygons; split on NaN separators
    function split_nans(l)
        ix = findall(isnan.(l[:,1]))
        [ [NTuple{2,Float64}(r) for r in eachrow(l[ix[i]+1:ix[i+1]-1,:])] for i = 1:length(ix)-1 ]
    end

    colors = mp.colors[]
    poly!(mp,  split_nans(elem_lines), color=colors[1])
    lines!(mp, int_lines, color=colors[2], linewidth=1)
    lines!(mp, bnd_lines, color=colors[3], linewidth=2)

    labels = mp.labels[]
    labels isa Symbol && (labels = (labels,))
    :nodes    in labels && scatter!(mp, m.x[:,1], m.x[:,2],
                                    color=colors[4], markersize=8,
                                    marker=:circle, strokewidth=1.0)
    :elements in labels && text!(mp, elem_mid[:,1], elem_mid[:,2],
                                  text=string.(1:size(elem_mid,1)),
                                  align=(:center,:center), fontsize=16)
    mp
end

"""
    SolutionPlot(msh, u; nsub=nothing, mesh_edges=false, fe=nothing)

Recipe drawing the FEM solution `u` on mesh `msh`.

- 1D: line plot of `u` vs. `x`.
- 2D: filled triangle plot (color-mapped). Pass `mesh_edges=true` to overlay
  interior and boundary edges in black.

`nsub` controls sub-element refinement (default: `1` for `p=1`, `3p` otherwise).
`fe` is the reference element whose nodes `u` lives on; it defaults to the
mesh element, and a solver element (for example on Gauss-Legendre nodes)
plots a field straight from the solver's node layout, see `viz_solution`.
"""
@recipe(SolutionPlot, msh, u) do scene
    Attributes(
        nsub = nothing,
        mesh_edges = false,
        fe = nothing,
    )
end

function Makie.plot!(sp::SolutionPlot)
    m, u = sp.msh[], sp.u[]
    fe = something(sp.fe[], m.fe)
    allx, allu, allel = viz_solution(m, u; nsub=sp.nsub[], fe)
    D = dim(m)

    if D == 1
        lines!(sp, allx[allel[:],1], allu[allel[:],1])
    elseif D == 2
        mesh!(sp, allx, allel', color=allu, shading=NoShading)
        if sp.mesh_edges[]
            _, int_lines, bnd_lines, _ = viz_mesh(m)
            lines!(sp, vcat(int_lines, bnd_lines), linewidth=1.0, color=:black)
        end
    else
        error("plot not implemented for D=$D")
    end
    sp
end

Makie.plottype(::HighOrderMesh) = MeshPlot
Makie.plottype(::HighOrderMesh, ::AbstractArray) = SolutionPlot

"""
    plot(m::HighOrderMesh; kw...)

Plot a mesh in a new `Figure`, with `DataAspect()` in 2D. Keyword arguments
are forwarded to the [`MeshPlot`](@ref) recipe.
"""
function Makie.plot(m::HighOrderMesh{D}; kw...) where {D}
    f  = Figure()
    ax = Axis(f[1,1]; aspect = D == 2 ? DataAspect() : nothing)
    meshplot!(ax, m; kw...)
    f
end

"""
    plot(m::HighOrderMesh, u::AbstractArray; kw...)

Plot the FEM solution `u` on `m` in a new `Figure`, with `DataAspect()` and a
`Colorbar` in 2D. Keyword arguments are forwarded to the [`SolutionPlot`](@ref)
recipe.
"""
function Makie.plot(m::HighOrderMesh{D}, u::AbstractArray; kw...) where {D}
    f  = Figure()
    ax = Axis(f[1,1]; aspect = D == 2 ? DataAspect() : nothing)
    p  = solutionplot!(ax, m, u; kw...)
    # Attach the colorbar to the mesh plot so that its limits and colormap
    # match the plotted (interpolated) values, which can overshoot the nodal ones.
    D == 2 && Colorbar(f[1,2], p.plots[1])
    f
end

# Exercise the plotting paths once at precompile time so a fresh session
# does not pay for their first-use latency. Never allowed to break loading.
if ccall(:jl_generating_output, Cint, ()) == 1
    try
        let m = ex1mesh(nref=1), u = ex1solution(m)
            plot(m); plot(m, u); plot(set_degree(mshline(3), 2), zeros(3, 3))
        end
    catch
    end
end

end
