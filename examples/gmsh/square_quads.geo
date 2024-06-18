SetFactory("OpenCASCADE");

Rectangle(1) = { 0,0,0, 1,1 };

Field[3] = Distance;
Field[3].NodesList = { 1,2,3,4 };
Field[4] = Threshold;
Field[4].IField = 3;
Field[4].LcMin = 0.02;
Field[4].LcMax = 0.3;
Field[4].DistMin = 0;
Field[4].DistMax = 0.5;
Background Field = 4;

// Mesh.Algorithm = 8;
Recombine Surface(1);
