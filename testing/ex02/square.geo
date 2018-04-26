// Flexible Plate

//-----------------------------------------
// Define Parameters
L_x = 1.000000; // X direction
L_y = 1.000000; // Y direction

// Number of elements in each direction
Nx = 10;
Ny = 10;

//-----------------------------------------
// GMSH Options
Geometry.PointNumbers = 1;
Geometry.LineNumbers = 1;
Geometry.SurfaceNumbers = 1;

Mesh.ElementOrder = 1;
Mesh.Algorithm = 6;
//Mesh.Algorithm3D = 4;
Mesh.Smoothing = 2;
Mesh.SecondOrderLinear = 0;
Mesh.PointNumbers = 1;

//-----------------------------------------
// Geometry
//
// Define base points, lines and surface
Point(1) = {    0,    0,0};
Point(2) = {  L_x,    0,0};
Point(3) = {  L_x,  L_y,0};
Point(4) = {    0,  L_y,0};
Line(1) = {1, 2};
Transfinite Line{1} = Nx+1;
Line(2) = {2, 3};
Transfinite Line{2} = Ny+1;
Line(3) = {3, 4};
Transfinite Line{3} = Nx+1;
Line(4) = {4, 1};
Transfinite Line{4} = Ny+1;
Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Transfinite Surface{1} = {1,2,3,4};
Recombine Surface {1};


//-----------------------------------------
// Generate Physical Outputs
Physical Surface("body") = {1};
Physical Line("BC_south") = {1};
Physical Line("BC_east") = {2};
Physical Line("BC_north") = {3};
Physical Line("BC_west") = {4};

