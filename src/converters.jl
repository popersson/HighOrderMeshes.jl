## ==============================================================================

## gmsh import

# Parse GMSH file into Dict #
function parse_gmsh(fname)
    fid = open(fname)
    if fid == -1
        error("Can't open file")
    end

    gmsh = Dict()
    
    readuntil(fid, "\$MeshFormat\n")
    str = readline(fid)
    if str != "2.2 0 8"
        close(fid)
        error("MeshFormat must be \"2.2 0 8\"")
    end
    readuntil(fid, "\$EndMeshFormat\n")

    next_section = readline(fid)
    if next_section == "\$PhysicalNames"
        nbrnames = parse(Int, readline(fid))
        names = []
        for ii = 1:nbrnames
            cname = split(readline(fid))
            push!(names, (parse(Int, cname[1]), parse(Int, cname[2]), cname[3]))
        end
        gmsh[:PhysicalNames] = names
        readuntil(fid, "\$EndPhysicalNames\n")
        next_section = readline(fid)
    else
        gmsh[:PhysicalNames] = []
        if !occursin(next_section, "\$Nodes")
          error("No \$Nodes section")
        end
    end

    np = parse(Int, readline(fid))
    p = zeros(3, np)
    for ii = 1:np
        row = split(readline(fid))
        nbr = parse(Int, row[1])
        if nbr != ii
            close(fid)
            error("Only supports contiguous node numbering")
        end
        p[:,ii] = parse.(Float64, row[2:4])
    end
    gmsh[:Nodes] = p

    readuntil(fid, "\$Elements\n")
    ne = parse(Int, readline(fid))
    e = Dict(:type=>Int[], :tags=>Vector{Int}[], :node_numbers=>zeros(Int,0,0))
    nn = []
    max_nn = 0
    for ii = 1:ne
        row = split(readline(fid))
        nbr = parse(Int, row[1])
        push!(e[:type], parse(Int, row[2]))
        nbrtags = parse(Int, row[3])
        if nbrtags != 2
            close(fid)
            error("Only supports exactly 2 tags")
        end
        push!(e[:tags], parse.(Int, row[4:5]))
        push!(nn, parse.(Int, row[6:end]))
        max_nn = max(max_nn, length(nn[end]))
    end
    e[:node_numbers] = zeros(Int, max_nn, ne)
    for ii = 1:ne
        e[:node_numbers][1:length(nn[ii]),ii] = nn[ii]
    end
    gmsh[:Elements] = e
    
    close(fid)
    return gmsh
end

## Create Mesh from GMSH file ##
"""
    gmsh2msh(gmsh_fname)

TBW
"""
function gmsh2msh(gmsh_fname)
  
    # GMsh element types
    gmshpnts = [15,15,15,15,15]
    gmshlins = [1,8,26,27,28]
    gmshtris = [2,9,21,23,25]
    gmshquads = [3,10,36,37,38]
    gmshtets = [4,11,29,30,31]
    gmshhexs = [5,12,92,93,94]
    
    # Map to FEMesh element types
    gmshtype2elgeom = Vector{ElementGeometry}(undef,100)
    gmshtype2elgeom[gmshpnts] .= [Block{0}()]
    gmshtype2elgeom[gmshlins] .= [Block{1}()]
    gmshtype2elgeom[gmshtris] .= [Simplex{2}()]
    gmshtype2elgeom[gmshquads] .= [Block{2}()]
    gmshtype2elgeom[gmshtets] .= [Simplex{3}()]
    gmshtype2elgeom[gmshhexs] .= [Block{3}()]

    # Map to porder
    gmshtype2porder = Vector{Int}(undef,100)
    gmshtype2porder[gmshpnts] .= 1
    gmshtype2porder[gmshlins] .= 1:5
    gmshtype2porder[gmshtris] .= 1:5
    gmshtype2porder[gmshquads] .= 1:5
    gmshtype2porder[gmshtets] .= 1:5
    gmshtype2porder[gmshhexs] .= 1:5
    
    # Node ordering
    gmsh_node_order_map(_) = throw("Unknown element")
    
    gmsh_node_order_map(::Block{1}) =
      [[1,2],
       [1,3,2],
       [1,3,4,2],
       [1,3,4,5,2],
       [1,3,4,5,6,2]
       ]

    gmsh_node_order_map(::Simplex{2}) =
      [[1,2,3],
       [1,4,2,6,5,3],
       [1,4,5,2,9,10,6,8,7,3],
       [1,4,5,6,2,12,13,14,7,11,15,8,10,9,3],
       [1,4,5,6,7,2,15,16,19,17,8,14,21,20,9,13,18,10,12,11,3]
       ]

    gmsh_node_order_map(::Simplex{3}) =
      [1:4,
       [1,5,2,7,6,3,8,10,9,4],
       [1,5,6,2,10,17,7,9,8,3,12,18,16,19,20,14,11,15,13,4],
       [1,5,6,7,2,13,23,25,8,12,24,9,11,10,3,16,26,27,22,29,35,33,31,34,19,15,28,21,30,32,18,14,20,17,4],
[1,5,6,7,8,2,16,29,34,31,9,15,32,33,10,14,30,11,13,12,3,20,35,38,36,28,41,53,54,48,46,55,51,43,49,24,19,40,39,27,44,56,50,45,52,23,18,37,26,42,47,22,17,25,21,4]
       ]

    gmsh_node_order_map(::Block{2}) = 
      [[1,2,4,3],
       [1,5,2,8,9,6,4,7,3],
       [1;5;6;2;12;13;14;7;11;16;15;8;4;10;9;3],
       [1;5;6;7;2;16;17;21;18;8;15;24;25;22;9;14;20;23;19;10;4;13;12;11;3],
       [1;5;6;7;8;2;20;21;25;26;22;9;19;32;33;34;27;10;18;31;36;35;28;11;17;24;30;29;23;12;4;16;15;14;13;3]
       ]
    
    gmsh_node_order_map(::Block{3}) = 
      [[1,2,4,3,5,6,8,7],
       [1;9;2;10;21;12;4;14;3;11;22;13;23;27;24;16;25;15;5;17;6;18;26;19;8;20;7],
       [1;9;10;2;11;33;36;15;12;34;35;16;4;20;19;3;13;37;38;17;41;57;58;45;44;60;59;46;23;50;49;21;14;40;39;18;42;61;62;48;43;64;63;47;24;51;52;22;5;25;26;6;27;53;54;29;28;56;55;30;8;32;31;7],
       [1;9;10;11;2;12;45;52;48;18;13;49;53;51;19;14;46;50;47;20;4;26;25;24;3;15;54;58;55;21;63;99;107;100;72;70;108;119;110;76;66;102;112;101;73;30;82;85;81;27;16;61;62;59;22;67;109;120;111;79;71;121;125;122;80;69;114;123;113;77;31;86;89;88;28;17;57;60;56;23;64;103;115;104;75;68;116;124;117;78;65;106;118;105;74;32;83;87;84;29;5;33;34;35;6;36;90;94;91;39;37;97;98;95;40;38;93;96;92;41;8;44;43;42;7],
       [1;9;10;11;12;2;13;57;68;67;60;21;14;61;69;72;66;22;15;62;70;71;65;23;16;58;63;64;59;24;4;32;31;30;29;3;17;73;77;78;74;25;89;153;161;162;154;105;100;163;185;188;167;109;99;164;186;187;168;110;92;156;172;171;155;106;37;122;126;125;121;33;18;84;85;86;79;26;93;165;189;190;169;116;101;193;209;210;197;117;104;196;212;211;198;118;98;175;202;201;173;111;38;127;134;133;132;34;19;83;88;87;80;27;94;166;192;191;170;115;102;194;213;214;200;120;103;195;216;215;199;119;97;176;203;204;174;112;39;128;135;136;131;35;20;76;82;81;75;28;90;157;177;178;158;108;95;179;205;206;181;114;96;180;208;207;182;113;91;160;184;183;159;107;40;123;129;130;124;36;5;41;42;43;44;6;45;137;141;142;138;49;46;148;149;150;143;50;47;147;152;151;144;51;48;140;146;145;139;52;8;56;55;54;53;7]
       ]

    gmsh = parse_gmsh(gmsh_fname)

    node_numbers = gmsh[:Elements][:node_numbers]
    element_tags = gmsh[:Elements][:tags]
    gmsh_element_types = gmsh[:Elements][:type]
    unique_gmsh_element_types = unique(gmsh_element_types)

    xcg = gmsh[:Nodes]
    elem_dims = dim.(gmshtype2elgeom[unique_gmsh_element_types])

    ndims = 3
    if maximum(elem_dims) ≤ 2 && maximum(abs.(xcg[3:end,:])) == 0
        ndims = 2
        xcg = xcg[1:2,:]
    end

    gmsh_eltype_vol = unique_gmsh_element_types[elem_dims.==ndims]
    gmsh_eltype_surf = unique_gmsh_element_types[elem_dims.==ndims - 1]

    length(gmsh_eltype_vol) > 1 && throw("Multiple element types not supported")
    length(gmsh_eltype_surf) > 1 && throw("Multiple element types not supported")

    elgeom = gmshtype2elgeom[gmsh_eltype_vol[1]]
    porder = gmshtype2porder[gmsh_eltype_vol[1]]

    t1 = node_numbers[:, gmsh_element_types.==gmsh_eltype_vol[1]]
    el = t1[gmsh_node_order_map(elgeom)[porder],:]

    m = HighOrderMesh(FiniteElement(elgeom,porder), xcg', el)


    ### Boundary tags
    tagcol = any(first.(element_tags) .== 0) ? 2 : 1
    surf_elgeom = gmshtype2elgeom[gmsh_eltype_surf[1]]
    surf_porder = gmshtype2porder[gmsh_eltype_surf[1]]
    if !isempty(surf_porder)
        length(surf_porder) > 1 && throw("Must have unique order surface elements")
        porder != surf_porder && throw("Must have same porder for volume and surface elements")

        surf_loc = findall(gmsh[:Elements][:type] .== gmsh_eltype_surf[1])
        nbr_surf_nodes = length(gmsh_node_order_map(surf_elgeom)[porder])
        surf_nodes = node_numbers[1:nbr_surf_nodes,surf_loc]
        sort!(surf_nodes, dims=1)

        f2n = mkface2nodes(m)
        nf,nel = size(m.nb)
        m_surf_nodes = [ m.el[f2n[:,j],iel] for j = 1:nf, iel = 1:nel if m.nb[j,iel][1] < 1 ]
        m_surf_index = [ (j,iel) for j = 1:nf, iel = 1:nel if m.nb[j,iel][1] < 1 ]
        sort!.(m_surf_nodes)

        surf_map = indexin(m_surf_nodes, eachcol(surf_nodes))
        for i in eachindex(surf_map)
            m.nb[m_surf_index[i]...] = (-element_tags[surf_loc[surf_map[i]]][tagcol],0,0)
        end
    end
    
    
    return m
end

## Run GMSH, convert to Mesh ##
"""
    rungmsh2msh(gmshfname; porder=1, cmdadd="")

TBW
"""
function rungmsh2msh(gmshfname; porder=1, cmdadd="")
    fout = tempname() * ".msh"
    if isa(cmdadd, String)
        cmdadd = split(cmdadd)
    end
    gmshcmd = `gmsh -3 -format msh2 -o $fout $gmshfname -order $porder $cmdadd`
    gmsh_output = read(gmshcmd, String)
    msh = gmsh2msh(fout)
    rm(fout)
    msh
end

function gmshstr2msh(geostr; porder=1, cmdadd="")
    fbase = tempname()
    fin = fbase * ".geo"
    fout = fbase * ".msh"

    open(fin,"w") do io
        println(io, geostr)
    end
    
    if isa(cmdadd, String)
        cmdadd = split(cmdadd)
    end
    gmshcmd = `gmsh -3 -format msh2 -o $fout $fin -order $porder $cmdadd`
    gmsh_output = read(gmshcmd, String)
    msh = gmsh2msh(fout)

    rm(fin)
    rm(fout)
    msh
end


## ==============================================================================

## VTK export

vtk_node_order_map(_) = throw("Unknown element")

vtk_node_order_map(::Simplex{2}) = [
        [1,2,3],
        [1,4,2,6,5,3],
        [1,4,5,2,9,10,6,8,7,3],
        [1,4,5,6,2,12,13,14,7,11,15,8,10,9,3],
        [1,4,5,6,7,2,15,16,19,17,8,14,21,20,9,13,18,10,12,11,3],
        [1,4,5,6,7,8,2,18,19,22,23,20,9,17,27,28,24,10,16,26,25,11,15,21,12,14,13,3]
]   
vtk_node_order_map(::Simplex{3}) = [
        [1,2,3,4],
        [1,5,2,7,6,3,8,9,10,4],
        [1,5,6,2,10,20,7,9,8,3,11,17,13,19,18,15,12,14,16,4],
        [1,5,6,7,2,13,32,34,8,12,33,9,11,10,3,14,23,24,17,29,35,28,31,26,20,15,25,18,30,27,21,16,19,22,4],
        [1,5,6,7,8,2,16,47,52,49,9,15,50,51,10,14,48,11,13,12,3,17,29,32,30,21,41,53,54,37,46,55,40,43,35,25,18,34,33,22,44,56,39,45,38,26,19,31,23,42,36,27,20,24,28,4],
        [1,5,6,7,8,9,2,19,65,73,72,67,10,18,68,74,71,11,17,69,70,12,16,66,13,15,14,3,20,35,38,39,36,25,55,75,79,76,47,63,81,80,52,62,77,53,57,45,30,21,43,44,40,26,58,82,83,51,64,84,54,61,48,31,22,42,41,27,59,78,50,60,49,32,23,37,28,56,46,33,24,29,34,4]
]

vtk_node_order_map(::Block{2}) = [
        [1,2,4,3],
        [1,5,2,8,9,6,4,7,3],
        [1,5,6,2,11,13,14,7,12,15,16,8,4,9,10,3],
        [1,5,6,7,2,14,17,18,19,8,15,20,21,22,9,16,23,24,25,10,4,11,12,13,3],
        [1,5,6,7,8,2,17,21,22,23,24,9,18,25,26,27,28,10,19,29,30,31,32,11,20,33,34,35,36,12,4,13,14,15,16,3],
        [1,5,6,7,8,9,2,20,25,26,27,28,29,10,21,30,31,32,33,34,11,22,35,36,37,38,39,12,23,40,41,42,43,44,13,24,45,46,47,48,49,14,4,15,16,17,18,19,3],
        [1,5,6,7,8,9,10,2,23,29,30,31,32,33,34,11,24,35,36,37,38,39,40,12,25,41,42,43,44,45,46,13,26,47,48,49,50,51,52,14,27,53,54,55,56,57,58,15,28,59,60,61,62,63,64,16,4,17,18,19,20,21,22,3]
]

vtk_node_order_map(::Block{3}) = [
        [1,2,4,3,5,6,8,7],
        [1,9,2,12,25,10,4,11,3,17,23,18,21,27,22,19,24,20,5,13,6,16,26,14,8,15,7],
        [1,9,10,2,15,49,50,11,16,51,52,12,4,13,14,3,25,41,42,27,33,57,58,37,34,59,60,38,29,45,46,31,26,43,44,28,35,61,62,39,36,63,64,40,30,47,48,32,5,17,18,6,23,53,54,19,24,55,56,20,8,21,22,7],
        [1,9,10,11,2,18,81,82,83,12,19,84,85,86,13,20,87,88,89,14,4,15,16,17,3,33,63,64,65,36,45,99,100,101,54,46,102,103,104,55,47,105,106,107,56,39,72,73,74,42,34,66,67,68,37,48,108,109,110,57,49,111,112,113,58,50,114,115,116,59,40,75,76,77,43,35,69,70,71,38,51,117,118,119,60,52,120,121,122,61,53,123,124,125,62,41,78,79,80,44,5,21,22,23,6,30,90,91,92,24,31,93,94,95,25,32,96,97,98,26,8,27,28,29,7],
        [1,9,10,11,12,2,21,121,122,123,124,13,22,125,126,127,128,14,23,129,130,131,132,15,24,133,134,135,136,16,4,17,18,19,20,3,41,89,90,91,92,45,57,153,154,155,156,73,58,157,158,159,160,74,59,161,162,163,164,75,60,165,166,167,168,76,49,105,106,107,108,53,42,93,94,95,96,46,61,169,170,171,172,77,62,173,174,175,176,78,63,177,178,179,180,79,64,181,182,183,184,80,50,109,110,111,112,54,43,97,98,99,100,47,65,185,186,187,188,81,66,189,190,191,192,82,67,193,194,195,196,83,68,197,198,199,200,84,51,113,114,115,116,55,44,101,102,103,104,48,69,201,202,203,204,85,70,205,206,207,208,86,71,209,210,211,212,87,72,213,214,215,216,88,52,117,118,119,120,56,5,25,26,27,28,6,37,137,138,139,140,29,38,141,142,143,144,30,39,145,146,147,148,31,40,149,150,151,152,32,8,33,34,35,36,7],
        [1,9,10,11,12,13,2,24,169,170,171,172,173,14,25,174,175,176,177,178,15,26,179,180,181,182,183,16,27,184,185,186,187,188,17,28,189,190,191,192,193,18,4,19,20,21,22,23,3,49,119,120,121,122,123,54,69,219,220,221,222,223,94,70,224,225,226,227,228,95,71,229,230,231,232,233,96,72,234,235,236,237,238,97,73,239,240,241,242,243,98,59,144,145,146,147,148,64,50,124,125,126,127,128,55,74,244,245,246,247,248,99,75,249,250,251,252,253,100,76,254,255,256,257,258,101,77,259,260,261,262,263,102,78,264,265,266,267,268,103,60,149,150,151,152,153,65,51,129,130,131,132,133,56,79,269,270,271,272,273,104,80,274,275,276,277,278,105,81,279,280,281,282,283,106,82,284,285,286,287,288,107,83,289,290,291,292,293,108,61,154,155,156,157,158,66,52,134,135,136,137,138,57,84,294,295,296,297,298,109,85,299,300,301,302,303,110,86,304,305,306,307,308,111,87,309,310,311,312,313,112,88,314,315,316,317,318,113,62,159,160,161,162,163,67,53,139,140,141,142,143,58,89,319,320,321,322,323,114,90,324,325,326,327,328,115,91,329,330,331,332,333,116,92,334,335,336,337,338,117,93,339,340,341,342,343,118,63,164,165,166,167,168,68,5,29,30,31,32,33,6,44,194,195,196,197,198,34,45,199,200,201,202,203,35,46,204,205,206,207,208,36,47,209,210,211,212,213,37,48,214,215,216,217,218,38,8,39,40,41,42,43,7],
        [1,9,10,11,12,13,14,2,27,225,226,227,228,229,230,15,28,231,232,233,234,235,236,16,29,237,238,239,240,241,242,17,30,243,244,245,246,247,248,18,31,249,250,251,252,253,254,19,32,255,256,257,258,259,260,20,4,21,22,23,24,25,26,3,57,153,154,155,156,157,158,63,81,297,298,299,300,301,302,117,82,303,304,305,306,307,308,118,83,309,310,311,312,313,314,119,84,315,316,317,318,319,320,120,85,321,322,323,324,325,326,121,86,327,328,329,330,331,332,122,69,189,190,191,192,193,194,75,58,159,160,161,162,163,164,64,87,333,334,335,336,337,338,123,88,339,340,341,342,343,344,124,89,345,346,347,348,349,350,125,90,351,352,353,354,355,356,126,91,357,358,359,360,361,362,127,92,363,364,365,366,367,368,128,70,195,196,197,198,199,200,76,59,165,166,167,168,169,170,65,93,369,370,371,372,373,374,129,94,375,376,377,378,379,380,130,95,381,382,383,384,385,386,131,96,387,388,389,390,391,392,132,97,393,394,395,396,397,398,133,98,399,400,401,402,403,404,134,71,201,202,203,204,205,206,77,60,171,172,173,174,175,176,66,99,405,406,407,408,409,410,135,100,411,412,413,414,415,416,136,101,417,418,419,420,421,422,137,102,423,424,425,426,427,428,138,103,429,430,431,432,433,434,139,104,435,436,437,438,439,440,140,72,207,208,209,210,211,212,78,61,177,178,179,180,181,182,67,105,441,442,443,444,445,446,141,106,447,448,449,450,451,452,142,107,453,454,455,456,457,458,143,108,459,460,461,462,463,464,144,109,465,466,467,468,469,470,145,110,471,472,473,474,475,476,146,73,213,214,215,216,217,218,79,62,183,184,185,186,187,188,68,111,477,478,479,480,481,482,147,112,483,484,485,486,487,488,148,113,489,490,491,492,493,494,149,114,495,496,497,498,499,500,150,115,501,502,503,504,505,506,151,116,507,508,509,510,511,512,152,74,219,220,221,222,223,224,80,5,33,34,35,36,37,38,6,51,261,262,263,264,265,266,39,52,267,268,269,270,271,272,40,53,273,274,275,276,277,278,41,54,279,280,281,282,283,284,42,55,285,286,287,288,289,290,43,56,291,292,293,294,295,296,44,8,45,46,47,48,49,50,7]
]

vtk_node_order_map(::FiniteElement{D,G,P,T}) where {D,G,P,T} =
    vtk_node_order_map(G())[P]

vtk_celltype(_) = throw("Unknown element")
vtk_celltype(::Simplex{2}) = 69
vtk_celltype(::Simplex{3}) = 71
vtk_celltype(::Block{2}) = 70
vtk_celltype(::Block{3}) = 72
vtk_celltype(::FiniteElement{D,G,P,T}) where {D,G,P,T} =
    vtk_celltype(G())

"""
    mk_vtk_data(fid, typename, uname, u)

TBW
"""
function mk_vtk_data(fid, typename, uname, u)
    if typename == "VECTORS"
        println(fid, "$(typename) $(uname) float")
        if size(u,2) == 2
            u = hcat(u, 0*u[:,1])
        end
    else
        println(fid, "$(typename) $(uname) float $(size(u,2))")
        println(fid, "LOOKUP_TABLE default")
    end
    write_matrix(fid, u)
end

"""
    vtkwrite(fname, m::HighOrderMesh{D,G,P,T}, u::Array{T}=T[]; umap=nothing) where {D,G,P,T}

    umap: Vector with indices into u for each VTK solution field

          Ex:
                m = ex1mesh()
                u = hcat(m.x[:,1].^2, m.x[:,2], -m.x[:,1]) # 3 components
                vtkwrite("ex1.vtk", m, u, umap=[1, 2:3])   # 2 VTK solution fields:
                                                           #   component 1   (scalar)
                                                           #   component 2-3 (vector)


TBW
"""
function vtkwrite(fname, m::HighOrderMesh{D,G,P,T}, u::Array{T}=T[]; umap=nothing) where {D,G,P,T}
    if size(u,1) == size(m.x,1)
        # Assume CG solution
        x = copy(m.x)
        el = copy(m.el)
    elseif size(u,1) == size(m.el,1)
        # Assume DG solution
        if size(u,3) == size(m.el,2) # Assume 3DG solution format
            u = convert_3dg_solution(m, u)
        end
        ns,nel = size(m.el,1),size(m.el,2)
        x = reshape(m.x[m.el,:], ns*nel, D)
        el = collect(reshape(1:ns*nel, ns, nel))
        if !isempty(u)
            u = reshape(u, ns*nel, :)
        end
    else
        throw("u: Unsupported dimensions")
    end
    if size(x,2) == 2
        x = hcat(x, 0*x[:,1])
    end
    map = vtk_node_order_map(m.fe)
    celltype = vtk_celltype(m.fe)
    el[map,:] = el

    nx = size(x,1)
    ns,nel = size(el,1),size(el,2)

    fid = open(fname, "w")
    print(fid,
        """
        # vtk DataFile Version 4.2
        From HighOrderMeshes.jl
        ASCII
        DATASET UNSTRUCTURED_GRID
        POINTS $nx float
        """)
    write_matrix(fid, x)
    println(fid, "CELLS $nel $(nel*(1+ns))")
    write_matrix(fid, hcat(fill(ns, nel), el' .- 1))
    println(fid, "CELL_TYPES $nel")
    write_matrix(fid, fill(celltype, nel))
    println(fid, "POINT_DATA $nx")

    if isnothing(umap)
        mk_vtk_data(fid, "SCALARS", "scalar_data", u)
    else
        for i in 1:length(umap)
            comp = umap[i]
            if isa(comp, Number)
                mk_vtk_data(fid, "SCALARS", "u$comp", u[:,comp])
            else
                nbr = string(( "$(c)" for c in comp )...)
                mk_vtk_data(fid, "VECTORS", "u$nbr", u[:,comp])
            end
        end
    end

    close(fid)
end

## ==============================================================================

eltype3dg(::ElementGeometry) = throw("Unsupported element type")
eltype3dg(::Simplex) = 0
eltype3dg(::Block) = 1

node_order_3dg(m::HighOrderMesh) = 1:nbr_ho_nodes(m.fe)

function node_order_3dg(m::HighOrderMesh{D,Simplex{D},P}) where {D,P}
    s3dg = ([ (i...,P-sum(i)) for i in Iterators.product(fill(0:P,D)...) if sum(i) <= P ])
    shom = ([ (P-sum(i),i...) for i in Iterators.product(fill(0:P,D)...) if sum(i) <= P ])
    ix = indexin(s3dg, shom)
end

# Convert to the mesh format in the DG-FEM package 3DG
# Returns 3DG mesh fields as named tuple
function mshto3dg(m::HighOrderMesh{D,G,P,T}) where {D,G,P,T}
    eltype,nv,nf = eltype3dg(G()), nvertices(G()), nfaces(G())

    dim = D
    porder = P

    if eltype == 0 # Simplex
        s = ref_nodes(Simplex{D}(), equispaced(P))
        if !isapprox(s, ref_nodes(m.fe,D))
            throw("3DG conversion only supported for standard simplex node order. Consider using set_degree(m, porder(m))")
        end
    elseif eltype == 1 # Block
        s0 = gauss_lobatto_nodes(P+1, T=T)
        if !isapprox(s0, ref_nodes(m.fe,1))
            throw("3DG only supports quad meshes with Lobatto nodes. Consider using set_lobatto_nodes")
        end
    end
    
    p1 = permutedims(dg_nodes(m), (1,3,2))
    
    m1 = set_degree(m, 1)
    s = eval_shapefcns(m1.fe, ref_nodes(m.fe, D))
    sbnd = eval_shapefcns(m1.fe, ref_nodes(m.fe, D-1))
    # s0 = eval_shapefcns(m1.fe, ref_nodes(m.fe, 1)) # Should be this for consistency
    s0 = ref_nodes(m.fe, 1)  # ... but 3DG seems to use coordinates
    p = Cdouble.(transpose(m1.x))
    t = Cint.(m1.el .- 1)

    zixmap(i) = i>0 ? i-1 : i
    t2t = [ Cint(zixmap(nb[1])) for nb in m1.nb ]
    t2n = [ Cint(nb[2] - 1) for nb in m1.nb ]
    # TODO: Add neighbor face permutation in t2n

    if eltype == 0 # Simplex
        # Simplex node order in 3DG different than HOM
        s = s[:,vcat(2:end,1)]
        sbnd = sbnd[:,vcat(2:end,1)]
        p1 = p1[node_order_3dg(m),:,:]
    end
    
    np = size(p,2)
    nt = size(t,2)
    ns = size(s,1)
    nsbnd = size(sbnd,1)
    ns0 = size(s0,1)

    ncurved = fill(Cuchar(1), 1, np);
    ecurved = fill(Cuchar(1), 4, nt);
    tcurved = fill(Cuchar(1), 1, nt);

    return (; dim, np, nt, ns, nsbnd, ns0, nv, nf, eltype, porder,
            p, p1, s, sbnd, s0, t, t2t, t2n, ncurved, ecurved, tcurved)
end

## ==============================================================================
