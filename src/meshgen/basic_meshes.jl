###########################################################################
## Structured hypercube meshes

"""
    blockmesh_hypercube(dims::NTuple{D,Int}; T=Float64)

Build a structured block mesh on `[0,1]^D` with `dims[d]` elements in each
direction. Returns a raw `(x, el)` pair (coordinates and connectivity) without
constructing a `HighOrderMesh`. Suitable as input to `mshhypercube`.
"""
function blockmesh_hypercube(dims::NTuple{D,Int}; T=Float64) where D
    x = zeros(T, 1, 0)
    for d in dims
        xx = (0:d) / T(d)
        x = hcat(repeat(x,length(xx),1),
                 reshape(repeat(xx',size(x,1),1),:,1))
    end

    el = [1]
    off = cumprod((0,dims...) .+ 1)
    for i = 1:D
        el2 = repeat(vcat(el, el .+ off[i]),1,1,dims[i]) .+ 
              off[i]*reshape(0:(dims[i]-1),1,1,dims[i])
        el2 = reshape(el2, 2*size(el,1), size(el,2)*dims[i])
        el = el2
    end

    x,el
end

"""
    mshhypercube(dims::NTuple{D,Int}; eg=Block{D}(), p=1, T=Float64, periodic_dirs=())
    mshcube(m, n, o; kwargs...)
    mshsquare(m, n; kwargs...)
    mshline(m; kwargs...)

Construct a `HighOrderMesh` on `[0,1]^D` with `dims[d]` cells per direction
and polynomial order `p`. Boundary regions are numbered `1..2D` in the order
`xmin, xmax, ymin, ymax, zmin, zmax`.

The element geometry `eg` is `Block{D}()` (quads, hexahedra) by default.
With `eg=Simplex{2}()` every cell is split into 2 triangles along the
diagonal from its lower left to its upper right corner, and with
`eg=Simplex{3}()` every cell is split into 6 tetrahedra around its main
diagonal (the Kuhn subdivision), so the faces of neighboring cells match.

Pass `periodic_dirs` to identify opposite face pairs periodically; e.g.
`periodic_dirs=(1,)` makes the mesh periodic in x.

`mshcube`, `mshsquare`, and `mshline` are convenience wrappers for 3D, 2D,
and 1D with default grid sizes of 5.

```julia
msh = mshsquare(8)                                # 8x8 quads, degree 1
msh = mshsquare(8, 4; eg=Simplex{2}(), p=3)       # 64 triangles, degree 3
msh = mshcube(4; eg=Simplex{3}(), p=2)            # 384 tetrahedra, degree 2
msh = mshsquare(8; periodic_dirs=(1,))            # periodic in x
msh = mshline(10; p=4)                            # 10 line elements on [0,1]
```
"""
function mshhypercube(dims::NTuple{D,Int}; eg::ElementGeometry=Block{D}(), p=1, T=Float64,
                      periodic_dirs=()) where D
    eg isa ElementGeometry{D} ||
        throw(ArgumentError("element geometry $eg does not match the dimension $D"))
    bndexpr(x) = [x'; 1 .- x'][:]
    x, el = blockmesh_hypercube(dims; T=T)
    m = HighOrderMesh(x, split_blocks(el, eg), bndexpr=bndexpr)
    for d in periodic_dirs
        set_bnd_periodic!(m, (2d-1, 2d), d)
    end
    set_degree(m, p)
end

# Split the cells of a block mesh (corners in lexicographic order) into
# positively oriented simplices that match across shared cell faces.
split_blocks(el, ::Block) = el
split_blocks(el, eg::Simplex) = reshape(el[block_split_pattern(eg), :], nvertices(eg), :)

# Two triangles sharing the diagonal from corner 1 to corner 4
block_split_pattern(::Simplex{2}) = [1 1; 2 4; 4 3]

# Kuhn subdivision: one tetrahedron per path from corner 1 to corner 8 along
# the coordinate directions, with two vertices swapped for odd permutations.
block_split_pattern(::Simplex{3}) = [1 1 1 1 1 1; 2 6 4 3 5 7; 4 2 3 7 6 5; 8 8 8 8 8 8]

mshcube(m=5, n=m, o=n; kwargs...)   = mshhypercube((m,n,o); kwargs...)
mshsquare(m=5, n=m; kwargs...)      = mshhypercube((m,n); kwargs...)
mshline(m=5; kwargs...)             = mshhypercube((m,); kwargs...)

###########################################################################
## Circle mesh

"""
    mshcircle(n=1; eg=Block{2}(), p=1, shape=:full)

Construct a high-order mesh of degree `p` on the unit disk, with all nodes
on the boundary exactly on the circle and the arcs split into equal angles.

- `eg=Block{2}()`: quads by transfinite interpolation. The coarsest mesh
  (`n=1`) consists of 3 quads (`:quarter`), 6 quads (`:half`), or 12 quads
  (`:full`); `n` subdivides each quad uniformly.
- `eg=Simplex{2}()`: triangles from a polygon of 2 (`:quarter`), 3 (`:half`)
  or 6 (`:full`) sectors, each split uniformly into `n^2` triangles. A
  smooth map takes the polygon onto the disk; it moves the nodes in
  proportion to the squared distance from the center, so elements near the
  center stay almost straight.
- `shape`: `:quarter` (first quadrant), `:half` (upper half), or `:full` (whole disk).

Boundary regions: `:quarter` → `[y=0, x=0, arc]`; `:half` → `[y=0, arc]`;
`:full` → `[arc]`.

```julia
msh = mshcircle(2; p=3)                       # 48 curved quads
msh = mshcircle(4; eg=Simplex{2}(), p=3)      # 96 curved triangles
msh = mshcircle(3; shape=:quarter, p=2)       # quarter disk, 3 boundary regions
```
"""
function mshcircle(n=1; eg::ElementGeometry{2}=Block{2}(), p=1, shape=:full)
    shape in (:quarter, :half, :full) ||
        throw(ArgumentError("`shape` must be :quarter, :half, or :full"))
    _mshcircle(eg, n, p, shape)
end

_circle_bndexpr(shape) = shape == :quarter ? (x -> [x[2], x[1], 0]) :
                         shape == :half    ? (x -> [x[2], 0]) : (x -> [0])

function _mshcircle(::Block{2}, n, p, shape)
    n1 = n * p
    s = (0:n1) ./ n1
    e = ones(n1+1)
    z = zeros(n1+1)
    a = (2 + π/8) / (4 + √2) # All coarsest quads same average (curved) side length

    # High-order node grid indices for the reference quad patch
    ix = reshape(1:(n1+1)^2, n1+1, n1+1)
    q0 = fill(0, p+1, p+1, n, n)
    for j = 1:n, i = 1:n
        ix0 = 0:p
        i0, j0 = (i-1)*p + 1, (j-1)*p + 1
        q0[:,:,i,j] .= ix[i0 .+ ix0, j0 .+ ix0]
    end
    q0 = reshape(q0, (p+1)^2, n^2)

    # Four boundary curves of the first quad patch (inner box → arc)
    phi = π/4 * (0:n1)/n1
    x1 = a*e;  y1 = a*s          # left edge (inner box)
    x2 = cos.(phi); y2 = sin.(phi)  # right edge (arc)
    x3 = range(a, stop=1,    length=n1+1); y3 = z   # bottom edge
    x4 = range(a, stop=1/√2, length=n1+1); y4 = x4  # top edge (diagonal)

    # Transfinite interpolation: blend the four boundary curves
    X1  = @. (1-s)*x1' + s*x2'
    Y1  = @. (1-s)*y1' + s*y2'
    X2  = @. x3*(1-s') + x4*s'
    Y2  = @. y3*(1-s') + y4*s'
    X12 = @. (1-s)*(1-s')*X1[1,1] + s*(1-s')*X1[end,1] + (1-s)*s'*X1[1,end] + s*s'*X1[end,end]
    Y12 = @. (1-s)*(1-s')*Y1[1,1] + s*(1-s')*Y1[end,1] + (1-s)*s'*Y1[1,end] + s*s'*Y1[end,end]
    X = @. X1 + X2 - X12
    Y = @. Y1 + Y2 - Y12
    p1 = [X[:] Y[:]]                      # first patch (lower-right sector)
    X2, Y2 = X[:,end:-1:1], Y[:,end:-1:1]
    p2 = [Y2[:] X2[:]]                    # second patch (upper-left sector, reflected)

    # Central square patch
    xx = a*(s.*e')
    yy = a*(e.*s')
    p0 = [xx[:] yy[:]]

    # Assemble the three quarter-circle patches and deduplicate shared nodes
    pp = [p1; p2; p0]
    el = [q0  q0 .+ (n1+1)^2  q0 .+ 2(n1+1)^2]
    pp, el = unique_mesh_nodes(pp, el)

    fe = FiniteElement(Block{2}(), p)
    if shape == :quarter
        return HighOrderMesh(fe, pp, el, bndexpr=_circle_bndexpr(shape))
    else
        # Rotate 90 degrees
        p1 = pp
        p2 = [-pp[:,2] pp[:,1]]
        el = hcat(el, el .+ size(p1,1))
        pp = vcat(p1,p2)
        pp,el = unique_mesh_nodes(pp,el)
        if shape == :half
            return HighOrderMesh(fe, pp, el, bndexpr=_circle_bndexpr(shape))
        else
            # Rotate 180 degrees
            p1 = pp
            p2 = [-pp[:,1] -pp[:,2]]
            el = hcat(el, el .+ size(p1,1))
            pp = vcat(p1,p2)
            pp,el = unique_mesh_nodes(pp,el)
            return HighOrderMesh(fe, pp, el, bndexpr=_circle_bndexpr(shape))
        end
    end
end

function _mshcircle(::Simplex{2}, n, p, shape)
    nsec = shape == :quarter ? 2 : shape == :half ? 3 : 6
    β    = shape == :quarter ? 1/4 : 1/3      # sector angle in units of π

    # Uniform subdivision of the reference sector (origin, v1, v2), lattice (i,j) ↦ (i v1 + j v2)/n
    ijs = [ (i, j) for j in 0:n for i in 0:n-j ]
    ix  = Dict(ij => k for (k, ij) in enumerate(ijs))
    t0  = Int[]
    for (i, j) in ijs
        i + j <= n-1 && append!(t0, (ix[(i,j)], ix[(i+1,j)], ix[(i,j+1)]))
        i + j <= n-2 && append!(t0, (ix[(i+1,j)], ix[(i+1,j+1)], ix[(i,j+1)]))
    end
    t0 = reshape(t0, 3, :)

    x  = zeros(0, 2)
    el = zeros(Int, 3, 0)
    for k in 0:nsec-1
        v1 = (cospi(k*β), sinpi(k*β))
        v2 = (cospi((k+1)*β), sinpi((k+1)*β))
        el = [el  t0 .+ size(x,1)]
        x  = [x; [ (i*v1[c] + j*v2[c]) / n for (i, j) in ijs, c in 1:2 ]]
    end
    x, el = unique_mesh_nodes(x, el)

    m = set_degree(HighOrderMesh(x, el; bndexpr=_circle_bndexpr(shape)), p)
    for r in eachrow(m.x)
        r .= _polygon_to_disk(r, β, nsec)
    end
    m
end

# Map a point of the polygon with vertices (cos kβπ, sin kβπ), k = 0..nsec, onto
# the unit disk. In the sector from v1 to v2, x = ρ ((1-t) v1 + t v2) moves by
# ρ² (u(t) - (1-t) v1 - t v2), with u(t) the point at the fraction t of the arc:
# the polygon edges map to the arcs with equal angles, the sector edges stay.
function _polygon_to_disk(x, β, nsec)
    iszero(x[1]) && iszero(x[2]) && return x
    θ  = mod(atan(x[2], x[1]) / π, 2)
    k  = min(floor(Int, θ / β), nsec - 1)
    c1, s1 = cospi(k*β), sinpi(k*β)
    c2, s2 = cospi((k+1)*β), sinpi((k+1)*β)
    det = c1*s2 - s1*c2
    a  = (x[1]*s2 - x[2]*c2) / det             # x = a v1 + b v2
    b  = (c1*x[2] - s1*x[1]) / det
    ρ  = a + b
    t  = b / ρ
    cu, su = cospi((k+t)*β), sinpi((k+t)*β)
    (x[1] + ρ^2 * (cu - (1-t)*c1 - t*c2),
     x[2] + ρ^2 * (su - (1-t)*s1 - t*s2))
end
