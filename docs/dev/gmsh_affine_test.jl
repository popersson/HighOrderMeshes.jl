# Affine round-trip check of the Gmsh node-order tables.
#
# For a straight-sided mesh every element is the affine (tris, tets) or
# bilinear (transfinite quads, hexes) image of the reference element, so the
# high-order nodes must equal the p=1 interpolation of the corner nodes. Any
# mistake in a Gmsh node-order table shows up as a mismatch. Requires the
# `gmsh` executable on the PATH.
#
# Used by step 10 of docs/dev/redesign-plan.md; run directly with
#   julia --project docs/dev/gmsh_affine_test.jl

using HighOrderMeshes

const geo_square = """
Point(1)={0,0,0}; Point(2)={1,0,0}; Point(3)={1,1,0}; Point(4)={0,1,0};
Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
Curve Loop(1)={1,2,3,4}; Plane Surface(1)={1};
Physical Curve("b",1)={1,2,3,4}; Physical Surface("d",1)={1};
Mesh.MeshSizeMax = 0.4;
"""
const geo_square_quads = geo_square *
    "Transfinite Curve{1,2,3,4}=4; Transfinite Surface{1}; Recombine Surface{1};\n"
const geo_box = """
SetFactory("OpenCASCADE");
Box(1)={0,0,0,1,1,1};
Physical Surface("b",1)={1,2,3,4,5,6}; Physical Volume("d",1)={1};
Mesh.MeshSizeMax = 0.6;
"""
const geo_box_hexes = """
Point(1)={0,0,0}; Point(2)={1,0,0}; Point(3)={1,1,0}; Point(4)={0,1,0};
Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,4}; Line(4)={4,1};
Curve Loop(1)={1,2,3,4}; Plane Surface(1)={1};
Transfinite Curve{1,2,3,4}=3; Transfinite Surface{1}; Recombine Surface{1};
out[] = Extrude {0,0,1} { Surface{1}; Layers{2}; Recombine; };
Physical Volume("d",1)={out[1]}; Physical Surface("b",1)={1,out[0],out[2],out[3],out[4],out[5]};
"""

"""
    gmsh_affine_error(geo, p) -> (mesh, max_error, nbnd_tagged, nbnd_untagged)

Mesh `geo` with Gmsh at order `p` and return the largest distance between the
imported high-order nodes and the linear interpolation of the corner nodes.
"""
function gmsh_affine_error(geo, p)
    m  = gmshstr2msh(geo; porder=p, cmdadd="-v 0")
    fe = m.fe
    xdg     = dg_nodes(m)
    corners = xdg[corner_nodes(fe), :, :]
    N1      = shapefcns(elgeom(m), ref_nodes(fe))
    pred    = interpolate(N1, corners)
    err     = maximum(abs.(pred - xdg))
    ntag    = count(nb -> nb[1] < 0, m.nb)
    nuntag  = count(nb -> nb[1] == 0, m.nb)
    m, err, ntag, nuntag
end

if abspath(PROGRAM_FILE) == @__FILE__
    for (label, geo) in (("tri", geo_square), ("quad", geo_square_quads),
                         ("tet", geo_box), ("hex", geo_box_hexes)), p in 1:5
        m, err, ntag, nuntag = gmsh_affine_error(geo, p)
        flag = err < 1e-10 ? "OK " : "MISMATCH"
        println("$flag $label p=$p: $(nel(m)) elements, max error $err, boundary faces tagged $ntag, untagged $nuntag")
    end
end
