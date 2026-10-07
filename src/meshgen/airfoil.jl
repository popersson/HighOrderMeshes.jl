###########################################################################
## Airfoil meshes
#
# mshairfoil meshes the region around an airfoil with curved quads, using
# gmsh for the unstructured part. A thin structured band of quads follows
# the airfoil and a structured wedge follows the wake line from the
# trailing edge to the right far field. The band keeps the wall layer free
# of irregular nodes, so the isoparametric layer refinement gives clean,
# curved wall-normal layers, and the wedge gives quads stretched along the
# wake. Notes on the gmsh setup:
# - gmsh places transfinite nodes by the spline parameter, not by arc
#   length, so every wall cell is its own short gmsh spline through points
#   of the airfoil spline. That puts the wall nodes, and the band edge nodes
#   along the normals, exactly where they are wanted. The band edge is a
#   polyline, so only the band cells have a curved side.
# - Blossom recombination leaves a few triangles next to the band, so the
#   mesh is made all-quad with Mesh.SubdivisionAlgorithm = 1: gmsh meshes
#   with twice the requested sizes and splits every element.
# - In gmsh 4.14 a MathEval field that refers to other MathEval fields
#   deadlocks, so the MathEval fields only refer to Distance fields.

const _AIRFOIL_DIR = joinpath(@__DIR__, "airfoils")

"""
    airfoil_coordinates(name::Symbol)
    airfoil_coordinates(fname::AbstractString)

Airfoil coordinates as an `n × 2` matrix, ordered from the trailing edge along
the upper surface to the leading edge and back along the lower surface to the
trailing edge (Selig order). `name` is one of the sample airfoils:

- `:naca0012`: NACA 0012, the point set of 3DG's NACA 0012 mesh.
- `:rae2822`: RAE 2822, from the UIUC airfoil coordinate database.

`fname` is a coordinate file in Selig or Lednicer format, as in the UIUC
airfoil coordinate database. Lines that are not two numbers are skipped.

```julia
X = airfoil_coordinates(:rae2822)
```
"""
function airfoil_coordinates(name::Symbol)
    fname = joinpath(_AIRFOIL_DIR, "$name.dat")
    isfile(fname) || error("unknown sample airfoil :$name, the samples are " *
                           join((":" * first(splitext(f)) for f in readdir(_AIRFOIL_DIR)), ", "))
    airfoil_coordinates(fname)
end

function airfoil_coordinates(fname::AbstractString)
    rows = Vector{Float64}[]
    for line in eachline(fname)
        v = tryparse.(Float64, split(line))
        length(v) == 2 && !any(isnothing, v) && push!(rows, v)
    end
    X = reduce(vcat, permutedims.(rows))
    if X[1,1] > 1.5 && X[1,2] > 1.5
        # Lednicer: the point counts, then both surfaces from the LE to the TE
        nu, nl = Int(X[1,1]), Int(X[1,2])
        up, lo = X[2:1+nu, :], X[2+nu:1+nu+nl, :]
        X = [up[end:-1:1, :]; lo[(lo[1,:] == up[1,:] ? 2 : 1):end, :]]
    end
    X
end

"""
    naca4(code; n=201)

Coordinates of the NACA 4-digit airfoil `code` (for example `"2412"`) with a
closed trailing edge and chord 1, in Selig order (see
[`airfoil_coordinates`](@ref)). Each surface has `n` points with cosine
spacing.

```julia
msh = mshairfoil(naca4("2412"), aoa=4)
```
"""
function naca4(code::AbstractString; n=201)
    length(code) == 4 && all(isdigit, code) || error("NACA 4-digit code expected, got \"$code\"")
    m, p, t = parse(Int, code[1]) / 100, parse(Int, code[2]) / 10, parse(Int, code[3:4]) / 100
    x  = @. (1 - cos(π * (0:n-1) / (n-1))) / 2
    yt = @. 5t * (0.2969sqrt(x) - 0.1260x - 0.3516x^2 + 0.2843x^3 - 0.1036x^4)
    yt[end] = 0                                  # closed TE, without rounding
    if m == 0 || p == 0
        yc, θ = zero(x), zero(x)
    else
        yc = @. ifelse(x < p, m / p^2 * (2p*x - x^2), m / (1-p)^2 * (1 - 2p + 2p*x - x^2))
        θ  = @. atan(ifelse(x < p, 2m / p^2 * (p - x), 2m / (1-p)^2 * (p - x)))
    end
    up = [x .- yt .* sin.(θ)  yc .+ yt .* cos.(θ)]
    lo = [x .+ yt .* sin.(θ)  yc .- yt .* cos.(θ)]
    [up[end:-1:1, :]; lo[2:end, :]]
end

###########################################################################
## Cubic splines

# Not-a-knot cubic spline through the rows of y at the increasing parameters
# t, stored as the second derivatives M at the knots.
struct _CubicSpline
    t::Vector{Float64}
    y::Matrix{Float64}
    M::Matrix{Float64}
end

function _CubicSpline(t::AbstractVector, y::AbstractMatrix)
    n = length(t)
    n >= 4 || error("a spline needs at least 4 points")
    h = diff(t)
    I, J, V = Int[], Int[], Float64[]
    add!(i, j, v) = (push!(I, i); push!(J, j); push!(V, v))
    rhs = zeros(n, size(y, 2))
    # not-a-knot: the third derivative is continuous at t[2] and t[n-1]
    add!(1, 1, h[2]); add!(1, 2, -(h[1] + h[2])); add!(1, 3, h[1])
    for i in 2:n-1
        add!(i, i-1, h[i-1]); add!(i, i, 2(h[i-1] + h[i])); add!(i, i+1, h[i])
        rhs[i,:] = 6 * ((y[i+1,:] - y[i,:]) / h[i] - (y[i,:] - y[i-1,:]) / h[i-1])
    end
    add!(n, n-2, h[n-1]); add!(n, n-1, -(h[n-2] + h[n-1])); add!(n, n, h[n-2])
    _CubicSpline(collect(float(t)), Matrix{Float64}(y), sparse(I, J, V, n, n) \ rhs)
end

# The spline interval of t, its length h, the offset u of t in it, and the
# linear coefficient b of the cubic on it.
function _spline_piece(sp::_CubicSpline, t)
    i  = clamp(searchsortedlast(sp.t, t), 1, length(sp.t) - 1)
    h  = sp.t[i+1] - sp.t[i]
    M0, M1 = sp.M[i,:], sp.M[i+1,:]
    b  = (sp.y[i+1,:] - sp.y[i,:]) / h - h * (2M0 + M1) / 6
    i, h, t - sp.t[i], b, M0, M1
end

function _spline(sp::_CubicSpline, t)
    i, h, u, b, M0, M1 = _spline_piece(sp, t)
    sp.y[i,:] + b*u + M0*u^2/2 + (M1 - M0)*u^3/(6h)
end

function _dspline(sp::_CubicSpline, t)
    _, h, u, b, M0, M1 = _spline_piece(sp, t)
    b + M0*u + (M1 - M0)*u^2/(2h)
end

###########################################################################
## Airfoil geometry

# The airfoil as a spline in the chord-length parameter, oriented
# counterclockwise from the TE, with a table of arc length against the
# parameter. An open trailing edge is closed by shearing each surface
# linearly (in x from the LE) towards the midpoint of the gap.
function _airfoil_spline(X::AbstractMatrix)
    X = Float64.(X)
    keep = [true; [X[i,:] != X[i-1,:] for i in 2:size(X,1)]]
    X = X[keep, :]
    n = size(X, 1)
    area = sum(X[i,1]*X[mod1(i+1,n),2] - X[mod1(i+1,n),1]*X[i,2] for i in 1:n) / 2
    area < 0 && (X = X[end:-1:1, :])
    ile = argmin(X[:,1])
    if X[1,:] != X[end,:]
        mid = (X[1,:] + X[end,:]) / 2
        for (rng, te) in ((1:ile, X[1,:]), (ile:n, X[end,:]))
            for i in rng
                λ = clamp((X[i,1] - X[ile,1]) / (te[1] - X[ile,1]), 0, 1)
                X[i,:] += λ * (mid - te)
            end
        end
    end
    t  = [0; cumsum(hypot.(diff(X[:,1]), diff(X[:,2])))]
    sp = _CubicSpline(t, X)
    # arc length table, with the knots among the samples
    nsamp = 20
    T = [ t[i] + (t[i+1] - t[i]) * k / nsamp for k in 0:nsamp-1, i in 1:n-1 ]
    T = [vec(T); t[end]]
    P = reduce(vcat, (_spline(sp, τ)' for τ in T))
    S = [0; cumsum(hypot.(diff(P[:,1]), diff(P[:,2])))]
    (; sp, T, S, tle = t[ile], sle = S[(ile-1)*nsamp + 1])
end

# The spline parameter at the arc lengths s from the start.
function _airfoil_param(g, s)
    map(s) do si
        k = clamp(searchsortedlast(g.S, si), 1, length(g.S) - 1)
        g.T[k] + (g.T[k+1] - g.T[k]) * (si - g.S[k]) / (g.S[k+1] - g.S[k])
    end
end

_unitvec(v) = v ./ hypot(v...)

# Outward unit normal (the curve is counterclockwise).
_airfoil_normal(g, t) = (d = _dspline(g.sp, t); _unitvec([d[2], -d[1]]))

# Solve f(y) = 0 by bisection, for f increasing on [lo, hi].
function _bisection(f, lo, hi)
    for _ = 1:200
        mid = (lo + hi) / 2
        f(mid) < 0 ? (lo = mid) : (hi = mid)
    end
    (lo + hi) / 2
end

# Vinokur's two-sided stretching: n+1 points on [0, 1] with first spacing
# about Δ0 and last spacing about Δ1.
function _vinokur(n, Δ0, Δ1)
    A = sqrt(Δ1 / Δ0)
    B = 1 / (n * sqrt(Δ0 * Δ1))
    ξ = (0:n) ./ n
    u = if B > 1 + 1e-6
        y = _bisection(y -> sinh(y)/y - B, 1e-9, 50.0)
        @. 0.5 + tanh(y*(ξ - 0.5)) / (2tanh(y/2))
    elseif B < 1 - 1e-6
        y = _bisection(y -> B - sin(y)/y, 1e-9, π - 1e-9)
        @. 0.5 + tan(y*(ξ - 0.5)) / (2tan(y/2))
    else
        collect(ξ)
    end
    @. u / (A + (1 - A) * u)
end

# Vinokur's stretching with the end spacings corrected to be exactly h0 and h1.
function _stretching(n, h0, h1)
    Δ0, Δ1 = h0, h1
    for _ = 1:50
        u = _vinokur(n, Δ0, Δ1)
        Δ0 *= h0 / (u[2] - u[1])
        Δ1 *= h1 / (u[end] - u[end-1])
    end
    _vinokur(n, Δ0, Δ1)
end

# Number of cells n and last spacing of a geometric progression with ratio r
# and first spacing h0 that covers the length L.
function _progression(L, h0, r)
    n = max(1, round(Int, log(1 + L*(r-1)/h0) / log(r)))
    n, L * (r-1) / (r^n - 1) * r^(n-1)
end

# The trailing edge and the unit wake direction.
_airfoil_wake(X, aoa) = (g = _airfoil_spline(X); (_spline(g.sp, 0.0), [cosd(aoa), sind(aoa)]))

###########################################################################
## Gmsh geometry

"""
    airfoil_geo(foil=:naca0012; aoa=0, nfoil=48, hle=0.004, hte=0.004,
                tband=0.02, nband=2, rband=1.5, hwake=hte, gwake=0.1, awake=4,
                Lwake=10, dwake=0, hnear=Inf, dnear=0, growth=0.2, hmax=20,
                R=100, gmshopts="")

Gmsh `.geo` string for the airfoil mesh of [`mshairfoil`](@ref), which
describes the parameters. Useful to inspect the mesh in the gmsh GUI.
"""
function airfoil_geo(foil=:naca0012; aoa=0.0, nfoil=48, hle=0.004, hte=0.004,
                     tband=0.02, nband=2, rband=1.5, hwake=hte, gwake=0.1, awake=4.0,
                     Lwake=10.0, dwake=0.0, hnear=Inf, dnear=0.0, growth=0.2, hmax=20.0,
                     R=100.0, gmshopts="")
    iseven(nfoil) && iseven(nband) || error("nfoil and nband must be even")
    X = foil isa AbstractMatrix ? foil : airfoil_coordinates(foil)
    g = _airfoil_spline(X)
    # gmsh works with half the element counts and twice the sizes
    K, nn, δ = nfoil ÷ 2, nband ÷ 2, tband
    nsub = 8                                   # spline points inside each wall cell

    # Wall nodes by arc length σ from the TE: K cells on each surface, with
    # the lengths hte and hle at the ends. Each cell gets nsub extra points.
    Lup, Ltot = g.sle, g.S[end]
    σ = [ g.sle .* _stretching(K, 2hte / Lup, 2hle / Lup);
          g.sle .+ (Ltot - Lup) .* _stretching(K, 2hle / (Ltot - Lup), 2hte / (Ltot - Lup))[2:end] ]
    Σ = [ σ[j] + (σ[j+1] - σ[j]) * k / (nsub + 1) for k in 0:nsub+1, j in 1:2K ]
    θ = reshape(_airfoil_param(g, vec(Σ)), size(Σ))
    θ[1,1], θ[end,K], θ[1,K+1], θ[end,end] = 0.0, g.tle, g.tle, g.T[end]  # TE, LE exactly
    wall = [ [_spline(g.sp, τ) for τ in θ[:,j]] for j in 1:2K ]

    # The band edge is at the distance tband along the normals. Near the TE
    # it turns to the bisector of the airfoil and wake normals, so that the
    # band and the wake wedge share the edge from the TE.
    TE = wall[1][1]
    w  = [cosd(aoa), sind(aoa)]
    nw = [-w[2], w[1]]
    bisect(n1, n2) = (b = _unitvec(n1 + n2); b / (b' * n1))
    bup = bisect(_airfoil_normal(g, 0.0), nw)
    blo = bisect(_airfoil_normal(g, g.T[end]), -nw)
    function band(j, k)
        dTE = min(Σ[k,j], Ltot - Σ[k,j])
        λ   = (1 - min(dTE / max(4δ, σ[2], Ltot - σ[end-1]), 1))^2
        b   = j <= K ? bup : blo
        wall[j][k] + δ * ((1 - λ) * _airfoil_normal(g, θ[k,j]) + λ * b)
    end
    offs = [ [band(j, 1), band(j, nsub+2)] for j in 1:2K ]

    # Wake wedge: the progression ratio 1+2gwake gives gmsh lengths
    # ≈ 2hwake + 2gwake*s. The thickness grows so that the aspect ratio tends
    # to awake at Lwake, and reaches isotropic elements at the far field.
    LE  = wall[K+1][1]
    mid = (LE + TE) / 2
    xL, xR, yB, yT = mid[1] - R, mid[1] + 2R, mid[2] - R, mid[2] + R
    rw   = 1 + 2gwake
    Lw_tot = (xR - TE[1]) / w[1]
    Lw   = min(Lwake, Lw_tot / 2)
    n1, h1 = _progression(Lw, 2hwake, rw)
    n2, h2 = _progression(Lw_tot - Lw, h1 * rw, rw)
    δmid = δ + Lw * nband * gwake / awake
    δfar = max(nn * h2, δmid)
    M    = TE + Lw * w
    W    = TE + Lw_tot * w
    abs(W[2] - mid[2]) + δfar / w[1] < 0.9R ||
        error("the wake does not fit the far field, reduce aoa or increase R")

    # Points are numbered as they are written; wall cells are the curves
    # 1001:1000+2K (splines) and band edge cells 2001:2000+2K (lines), both
    # counterclockwise from the TE.
    io  = IOBuffer()
    npt = Ref(0)
    function pt(x)
        println(io, "Point($(npt[] += 1)) = { $(x[1]), $(x[2]), 0 };")
        npt[]
    end
    function curves(cells, tag0; closed)
        ids = [ pt(cells[1][1]) ]
        for (j, c) in enumerate(cells)
            inner = [ pt(p) for p in c[2:end-1] ]
            last  = closed && j == length(cells) ? ids[1] : pt(c[end])
            println(io, isempty(inner) ? "Line" : "Spline", "($(tag0+j)) = { ",
                    join([ids[end]; inner; last], ", "), " };")
            push!(ids, last)
        end
        ids
    end
    wid = curves(wall, 1000; closed=true)
    oid = curves(offs, 2000; closed=false)
    iTE, iLE = wid[1], wid[K+1]
    iPup, iOLE, iPlo = oid[1], oid[K+1], oid[end]
    iM, iMup, iMlo = pt(M), pt(M + δmid * nw), pt(M - δmid * nw)
    iW, iWup, iWlo = pt(W), pt(W + [0, δfar / w[1]]), pt(W - [0, δfar / w[1]])
    iF = [ pt(p) for p in ([xL, yB], [xR, yB], [xR, yT], [xL, yT]) ]
    wup, wlo = "1001:$(1000+K)", "$(1001+K):$(1000+2K)"
    oup, olo = "2001:$(2000+K)", "$(2001+K):$(2000+2K)"
    rev(r) = join(-reverse(r), ", ")           # reversed curves of a range
    hband(j) = hypot((offs[j][2] - offs[j][1])...)
    hmid   = maximum(diff(σ))
    hwk    = "$(2hwake) + $(2gwake)*F3"                          # wake spacing
    fields = [ "$hmid + $(2growth)*F1",                          # band edge
               "$(hband(K)) + $(2growth)*F2",                    # LE
               "$(hband(1)) + $(2growth)*F3",                    # TE
               "$hwk + $(2growth)*F4" ]                          # wedge edge
    dwake > 0 && push!(fields, "Max($hwk, $(isfinite(hnear) ? 2hnear : hmid)) + " *
                               "$(2growth)*Max(F4 - $dwake, 0)")
    isfinite(hnear) && push!(fields, "$(2hnear) + $(2growth)*Max(F1 - $dnear, 0)")
    print(io, """
    Line(5) = { $iLE, $iOLE };                 // band normal at the LE
    Line(6) = { $iTE, $iPup };                 // band normals at the TE
    Line(7) = { $iTE, $iPlo };
    Line(8) = { $iTE, $iM };                   // near wake line
    Line(9) = { $iPup, $iMup };                // near wake wedge edges
    Line(10) = { $iPlo, $iMlo };
    Line(11) = { $iM, $iMup };                 // wedge normals at Lwake
    Line(12) = { $iM, $iMlo };
    Line(13) = { $iM, $iW };                   // far wake line
    Line(14) = { $iMup, $iWup };               // far wake wedge edges
    Line(15) = { $iMlo, $iWlo };
    Line(16) = { $iW, $iWup };                 // wedge ends on the far field
    Line(17) = { $iW, $iWlo };
    Line(21) = { $(iF[1]), $(iF[2]) };         // far field
    Line(22) = { $(iF[2]), $iWlo };
    Line(23) = { $iWup, $(iF[3]) };
    Line(24) = { $(iF[3]), $(iF[4]) };
    Line(25) = { $(iF[4]), $(iF[1]) };

    Curve Loop(1) = { $(rev(1001:1000+K)), 6, $oup, -5 };       // band, upper
    Curve Loop(2) = { $(rev(1001+K:1000+2K)), 5, $olo, -7 };    // band, lower
    Curve Loop(3) = { 8, 11, -9, -6 };         // near wake, upper
    Curve Loop(4) = { 7, 10, -12, -8 };        // near wake, lower
    Curve Loop(5) = { 13, 16, -14, -11 };      // far wake, upper
    Curve Loop(6) = { 12, 15, -17, -13 };      // far wake, lower
    Curve Loop(7) = { 21, 22, -15, -10, $(rev(2001:2000+2K)), 9, 14, 23, 24, 25 };
    For i In {1:7}
      Plane Surface(i) = { i };
    EndFor

    Transfinite Curve{ $wup, $wlo, $oup, $olo } = 2;
    Transfinite Curve{ 5, 6, 7, 11, 12, 16, 17 } = $(nn+1) Using Progression $rband;
    Transfinite Curve{ 8, 9, 10 } = $(n1+1) Using Progression $rw;
    Transfinite Curve{ 13, 14, 15 } = $(n2+1) Using Progression $rw;
    Transfinite Surface{ 1 } = { $iLE, $iTE, $iPup, $iOLE };
    Transfinite Surface{ 2 } = { $iTE, $iLE, $iOLE, $iPlo };
    Transfinite Surface{ 3 } = { $iTE, $iM, $iMup, $iPup };
    Transfinite Surface{ 4 } = { $iTE, $iPlo, $iMlo, $iM };
    Transfinite Surface{ 5 } = { $iM, $iW, $iWup, $iMup };
    Transfinite Surface{ 6 } = { $iM, $iMlo, $iWlo, $iW };

    Physical Curve("Airfoil", 1) = { $wup, $wlo };
    Physical Curve("Far field", 2) = { 21, 22, 17, 16, 23, 24, 25 };
    Physical Surface("Domain", 1) = { 1:7 };

    // Sizes in the unstructured region, growing at the rate growth
    Field[1] = Distance;
    Field[1].CurvesList = { $oup, $olo };
    Field[1].Sampling = 10;
    Field[2] = Distance;
    Field[2].PointsList = { $iOLE };
    Field[3] = Distance;
    Field[3].PointsList = { $iTE };
    Field[4] = Distance;
    Field[4].CurvesList = { 9, 10, 14, 15 };
    Field[4].Sampling = 1000;
    """)
    for (i, f) in enumerate(fields)
        println(io, "Field[$(10+i)] = MathEval;\nField[$(10+i)].F = \"$f\";")
    end
    print(io, """
    Field[20] = Min;
    Field[20].FieldsList = { $(join(10 .+ eachindex(fields), ", ")) };
    Background Field = 20;
    Mesh.MeshSizeMax = $(2hmax);
    Mesh.MeshSizeExtendFromBoundary = 0;
    Mesh.MeshSizeFromPoints = 0;
    Mesh.MeshSizeFromCurvature = 0;

    Mesh.Algorithm = 6;
    Mesh.RecombineAll = 1;
    Mesh.SubdivisionAlgorithm = 1;
    Mesh.HighOrderOptimize = 0;
    $gmshopts
    """)
    String(take!(io))
end

###########################################################################
## Airfoil mesh

# C-mesh style layers: in pass i, split every edge with exactly one vertex
# on the airfoil (boundary 1) or on the wake line within lengths[i] of the
# TE. The airfoil and the wake line count as one curve, so edges along them
# are never split and the layers continue through the TE; a layer that stops
# in the wake ends with one 3-way split on each side of the line. A layer
# always reaches the first wake node, so that none ends at the TE.
function _airfoil_layer_refine(m::HighOrderMesh{2,Block{2}}, TE, w, lengths)
    fmap = facemap(Block{2}())
    for L in lengths
        cel = m.el[corner_nodes(m.fe), :]
        oncurve = falses(size(m.x, 1))
        for iel in axes(m.nb,2), j in axes(m.nb,1)
            nb = m.nb[j,iel]
            isboundary(nb) && bndtag(nb) == 1 && (oncurve[cel[fmap[:,j],iel]] .= true)
        end
        d = m.x .- TE'
        s = d * w                                     # distance along the wake
        online = @. abs(d[:,1]*w[2] - d[:,2]*w[1]) < 1e-10 * (1 + abs(s)) && s > 1e-12
        oncurve .|= online .& (s .<= max(L, minimum(s[online])) + 1e-12)
        marked = [ count(oncurve[cel[fmap[:,j],iel]]) == 1 for j in axes(fmap,2), iel in axes(cel,2) ]
        m = refine(m, marked)
    end
    m
end

"""
    mshairfoil(foil=:naca0012; aoa=0, p=3, ref=0, nbndlayers=6,
               layerlength=10, layerratio=0.5, verbose=false, kwargs...)

Quad mesh of degree `p` around an airfoil, generated with gmsh (which
must be on the `PATH`). Boundary 1 is the airfoil and boundary 2 is the far
field. `foil` is a sample name or coordinate file for
[`airfoil_coordinates`](@ref), or a coordinate matrix in the same format,
for example from [`naca4`](@ref). The coordinates are in chord units with
the leading edge near the origin; an open trailing edge is closed.

A thin structured band of quads follows the airfoil, and a structured
wedge follows a straight wake line from the trailing edge at the angle `aoa`
(degrees) to the right far field. The rest is unstructured. The gmsh mesh is
refined uniformly `ref` times and then gets `nbndlayers` isoparametric
boundary layers, which continue into the wake as in a C-mesh: layer `i` ends
`layerlength*layerratio^(i-1)` chords behind the trailing edge, so the
layers thin out downstream and the far field keeps isotropic elements. Use
`layerlength=Inf, layerratio=1` for layers all the way to the far field, or
`layerlength=0` to end them right behind the trailing edge.

The gmsh mesh is set by the `kwargs` (sizes and counts before `ref` and the
layers):

- `nfoil=48`, `hle=0.004`, `hte=0.004`: elements along each surface (even)
  and their lengths at the leading and trailing edges.
- `tband=0.02`, `nband=2`, `rband=1.5`: thickness of the structured band, its
  number of layers (even), and the thickness ratio of successive pairs of
  layers.
- `hwake=hte`, `gwake=0.1`: the element length along the wake is about
  `hwake + gwake*s` at the distance `s` behind the trailing edge.
- `awake=4`, `Lwake=10`: the wake elements approach the aspect ratio `awake`
  (length over thickness) up to `Lwake` chords behind the trailing edge and
  become isotropic towards the far field. `awake=1` gives isotropic wake
  elements.
- `dwake=0`: the wake element length (but at least `hnear`, or the mid-chord
  wall spacing) is kept up to the distance `dwake` from the wedge.
- `hnear=Inf`, `dnear=0`: the element size is at most `hnear` up to the
  distance `dnear` from the band, as needed for wall-resolved LES.
- `growth=0.2`, `hmax=20`: away from the band and the wake the element size
  grows at the rate `growth`, up to `hmax`.
- `R=100`: the far field is the rectangle `[xm-R, xm+2R] × [ym-R, ym+R]`
  around the mid-chord point `(xm, ym)`.
- `gmshopts=""`: extra gmsh commands appended to the `.geo` file (see
  [`airfoil_geo`](@ref)), for example `"Mesh.Smoothing = 10;"`.

With `verbose`, the gmsh output and the boundary names are printed.

```julia
msh = mshairfoil(aoa=5)                                    # NACA 0012, RANS
msh = mshairfoil(:rae2822, aoa=2.79, nbndlayers=10)
msh = mshairfoil(naca4("2412"), aoa=4, layerlength=Inf, layerratio=1)
```
"""
function mshairfoil(foil=:naca0012; aoa=0.0, p=3, ref=0, nbndlayers=6,
                    layerlength=10.0, layerratio=0.5, verbose=false, kwargs...)
    X   = foil isa AbstractMatrix ? foil : airfoil_coordinates(foil)
    m   = gmshstr2msh(airfoil_geo(X; aoa, kwargs...); p, verbose)
    m   = uniref(set_lobatto_nodes(m), ref)
    TE, w = _airfoil_wake(X, aoa)
    _airfoil_layer_refine(m, TE, w, [ layerlength * layerratio^(i-1) for i in 1:nbndlayers ])
end
