// NACA mesh

// Wing
Include "naca0012.geo";
Spline(1) = { 1:53 };
Spline(2) = { 53:104,1 };
Line Loop(1) = { 1,2 };

// Outer boundary rectangle
R = 100;
Point(111) = { .5-R,-R,0 };
Point(112) = { .5+2*R,-R,0 };
Point(113) = { .5+2*R, R,0 };
Point(114) = { .5-R, R,0 };

Line(5) = { 111, 112 };
Line(6) = { 112, 113 };
Line(7) = { 113, 114 };
Line(8) = { 114, 111 };

Line Loop(2) = { 5, 6, 7, 8 };

Point(201) = { 2.5,0,0 };
Line(11) = { 1, 201 };

Point(202) = { 5,0,0 };
Line(12) = { 201,202 };

// Final geometry
Plane Surface(1) = { 2,1 };
//Line{11} In Surface{1};
//Line{12} In Surface{1};

Physical Line("Airfoil", 1) = {1,2};
Physical Line("Far field", 2) = {5,6,7,8};
Physical Surface("Domain", 1) = {1};

// Size field

hLE = 0.010485;
hTE = 0.01;
hwing = 0.035;
hwake1 = 0.05;
hwake2 = 0.3;

delta_wing = 0.15;
delta_wake1 = 0.2;
delta_wake2 = 1;

hgrw = 31;
dgrw = 100;

Field[1] = Attractor;
Field[1].NNodesByEdge = 100;
Field[1].EdgesList = { 1,2 };

Field[2] = Threshold;
Field[2].IField = 1;
Field[2].LcMin = hwing;
Field[2].LcMax = hwing + hgrw;
Field[2].DistMin = delta_wing;
Field[2].DistMax = delta_wing + dgrw;
Field[2].StopAtDistMax = 0;

Field[3] = Attractor;
Field[3].NodesList = { 53 };

Field[4] = Threshold;
Field[4].IField = 3;
Field[4].LcMin = hLE;
Field[4].LcMax = hLE + hgrw;
Field[4].DistMin = 0;
Field[4].DistMax = dgrw;
Field[4].StopAtDistMax = 0;

Field[5] = Attractor;
Field[5].NodesList = { 1 };

Field[6] = Threshold;
Field[6].IField = 5;
Field[6].LcMin = hTE;
Field[6].LcMax = hTE + hgrw;
Field[6].DistMin = 0;
Field[6].DistMax = dgrw;
Field[6].StopAtDistMax = 0;

Field[7] = Attractor;
Field[7].NNodesByEdge = 100;
Field[7].EdgesList = { 11 };

Field[8] = Threshold;
Field[8].IField = 7;
Field[8].LcMin = hwake1;
Field[8].LcMax = hwake1 + hgrw;
Field[8].DistMin = delta_wake1;
Field[8].DistMax = delta_wake1 + dgrw;
Field[8].StopAtDistMax = 1;

Field[9] = Attractor;
Field[9].NNodesByEdge = 100;
Field[9].EdgesList = { 12 };

Field[10] = Threshold;
Field[10].IField = 9;
Field[10].LcMin = hwake2;
Field[10].LcMax = hwake2 + hgrw;
Field[10].DistMin = delta_wake2;
Field[10].DistMax = delta_wake2 + dgrw;
Field[10].StopAtDistMax = 1;

Field[11] = Min;
Field[11].FieldsList = { 2,4,6,8,10 };
Background Field = 11;

//Define Boundary Layer
//Field[21] = BoundaryLayer;
//Field[21].EdgesList = {5, 227, 208, -226, 206, 6};
//Field[21].NodesList = {1,1021, 1022, 160};
//Field[21].FacesList = {232, 228};
//Field[21].FanNodesList = {1,1021, 1022, 160};
//Field[21].FansList = {5, 6};
//Field[21].hfar = 0.05;
//Field[21].hwall_n = 0.0005;
//Field[21].thickness = 0.02;
//Field[21].ratio = 1.1;
//Field[21].AnisoMax = 10;
//Field[21].Quads = 1;
//Field[21].IntersectMetrics = 0;
//BoundaryLayer Field = 21;

Mesh.Algorithm = 8;
Mesh.HighOrderOptimize = 1;
Recombine Surface(1);
Mesh.CharacteristicLengthExtendFromBoundary = 0;

Mesh 2;
