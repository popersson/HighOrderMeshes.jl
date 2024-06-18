abstract type ElementGeometry{D} end
struct Simplex{D} <: ElementGeometry{D} end
struct Block{D} <: ElementGeometry{D} end

dim(::ElementGeometry{D}) where {D} = D::Int

name(::Simplex{D}) where {D} = D > 3 ? "$(D)D simplex" :
    ("point", "line", "triangle", "tetrahedron")[D+1]
name(::Block{D}) where {D} = D > 3 ? "$(D)D simplex" :
    ("point", "line", "quadrilateral", "hexahedron")[D+1]

function Base.show(io::IO, eg::ElementGeometry)
    print(io, "ElementGeometry: $(dim(eg))D $(name(eg))")
end

nvertices(::ElementGeometry) = throw("Not implemented")
nvertices(::Simplex{D}) where {D} = D + 1
nvertices(::Block{D}) where {D} = 2^D

nfaces(::ElementGeometry) = throw("Not implemented")
nfaces(::Simplex{D}) where {D} = D + 1
nfaces(::Block{D}) where {D} = 2*D

nedges(::ElementGeometry) = throw("Not implemented")
nedges(::Simplex{D}) where {D} = binomial(D + 1, 2)
nedges(::Block{D}) where {D} = D * 2^(D - 1)

facemap(::ElementGeometry) = throw("Not implemented")
facemap(::Simplex{1}) = [2 1]
facemap(::Simplex{2}) = [[2,3] [3,1] [1,2]]
facemap(::Simplex{3}) = [[2,3,4] [1,4,3] [4,1,2] [3,2,1]]
facemap(::Block{1}) = [1 2]
facemap(::Block{2}) = [[3,1] [2,4] [1,2] [4,3]]
facemap(::Block{3}) = [[3,1,7,5] [2,4,6,8] [1,2,5,6] #=
                    =# [4,3,8,7] [3,4,1,2] [5,6,7,8]]

plot_face_order(::ElementGeometry) = throw("Not imlemented")
plot_face_order(::Simplex{2}) = [1, 2, 3]
plot_face_order(::Block{2}) = [3, 2, 4, 1]

edgemap(::ElementGeometry) = throw("Not implemented")
edgemap(::Simplex{1}) = [2,1]
edgemap(::Simplex{2}) = [[2,3] [3,1] [1,2]]
edgemap(::Simplex{3}) = [[1,2] [1,3] [1,4] [2,3] [2,4] [3,4]]
edgemap(::Block{1}) = [2,1]
edgemap(::Block{2}) = [[1,2] [2,4] [4,3] [3,1]]
edgemap(::Block{3}) = [[1,2] [2,4] [4,3] [3,1] [1,5] [2,6] #=
                    =# [4,8] [3,7] [5,6] [6,8] [8,7] [7,5]]

linegeom(::Simplex{D}) where {D} = Simplex{1}()
linegeom(::Block{D}) where {D} = Block{1}()
facegeom(::Simplex{D}) where {D} = Simplex{D-1}()
facegeom(::Block{D}) where {D} = Block{D-1}()

midpoint(::Simplex{D}) where D = fill(1/(D+1), D)
midpoint(::Block{D}) where D = fill(1/2, D)

function find_elgeom(D, nv)
    if nv == nvertices(Block{D}())
        return Block{D}()
    elseif nv == nvertices(Simplex{D}())
        return Simplex{D}()
    else
        error("Cannot determine element shape")
    end
end
