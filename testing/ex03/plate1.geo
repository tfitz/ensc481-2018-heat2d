// Gmsh project created on Fri Apr 27 14:35:05 2018
SetFactory("OpenCASCADE");


Geometry.PointNumbers = 1;
Geometry.LineNumbers = 1;
Geometry.SurfaceNumbers = 1;

Mesh.ElementOrder = 1;
Mesh.Smoothing = 2;
Mesh.SecondOrderLinear = 0;
Mesh.PointNumbers = 1;

//-----------------------------------------
// Geometry
//
Rectangle(1) = {0, 0, 0, 3, 1, 0};
Ellipse(5) = {1.5, 0.5, 0, 0.7, 0.3, 0, 2*Pi};
Line Loop(2) = {5};
Surface(2) = {2};
BooleanDifference(3) = { Surface{1}; Delete; }{ Surface{2}; Delete; };
// set relative size of elements
Characteristic Length{1,2,3,4,5} = 0.25;
// recombine to get quads
Recombine Surface {3};


//-----------------------------------------
// Generate Physical Outputs
Physical Surface("body") = {3};
Physical Line("BC_south") = {1};
Physical Line("BC_east") = {3};
Physical Line("BC_north") = {4};
Physical Line("BC_west") = {2};
