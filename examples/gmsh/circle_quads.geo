R = 1;

Point(1) = { .5,0,0 };
Point(2) = { .5 - R, 0, 0 };
Point(3) = { .5 + R, 0, 0 };
Circle(4) = { 2, 1, 3 };
Circle(5) = { 3, 1, 2 };
Line Loop(1) = { 4, 5 };

Plane Surface(1) = { 1 };

Physical Line("Circle", 1) = {4,5};
Physical Surface("Domain", 1) = {1};

//Mesh.Algorithm = 8;
Recombine Surface(1);
