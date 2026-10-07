module HighOrderMeshesDistMeshExt

using HighOrderMeshes
using DistMesh: DMesh

# Straight-sided mesh from a DistMesh mesh; documented with HighOrderMesh
function HighOrderMeshes.HighOrderMesh(dm::DMesh{D,T,N}; kwargs...) where {D,T,N}
    x  = T[ p[k] for p in dm.p, k in 1:D ]
    el = Int[ t[i] for i in 1:N, t in dm.t ]
    HighOrderMesh(x, el; kwargs...)
end

end
