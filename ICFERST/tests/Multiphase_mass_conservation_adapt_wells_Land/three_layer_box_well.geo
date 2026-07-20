// Three-layer box WITH EMBEDDED WELL-TRACE NODES for the well saturation
// conservation-across-adapt test. 
//
// Generate with (gmsh v2 format):
//   gmsh -3 -format msh2 -o three_layer_box_well.msh three_layer_box_well.geo

SetFactory("OpenCASCADE");

Box(1) = {0, 0,  0, 200, 200, 60};   // bottom
Box(2) = {0, 0, 60, 200, 200, 20};   // barrier
Box(3) = {0, 0, 80, 200, 200, 40};   // top

Coherence;

// ---- well trace vertices (must match co2_toy_well_welltrace.bdf) ----
lc = 10;
Point(101) = {100, 100,  20, lc};
Point(102) = {100, 100,  30, lc};
Point(103) = {100, 100,  40, lc};
Point(104) = {100, 100,  50, lc};
Point(105) = {100, 100,  60, lc};   // on the region-1/2 interface
Point(106) = {100, 100,  70, lc};
Point(107) = {100, 100,  80, lc};   // on the region-2/3 interface
Point(108) = {100, 100,  90, lc};
Point(109) = {100, 100, 100, lc};
Point(110) = {100, 100, 110, lc};
Point(111) = {100, 100, 120, lc};   // on the top surface

// interior points into their volumes
Point{101, 102, 103, 104} In Volume{1};
Point{106}                In Volume{2};
Point{108, 109, 110}      In Volume{3};

// interface / boundary points into their surfaces
s12() = Surface In BoundingBox{-1, -1, 59.9, 201, 201, 60.1};
s23() = Surface In BoundingBox{-1, -1, 79.9, 201, 201, 80.1};
stp() = Surface In BoundingBox{-1, -1, 119.9, 201, 201, 120.1};
Point{105} In Surface{s12(0)};
Point{107} In Surface{s23(0)};
Point{111} In Surface{stp(0)};

Physical Volume(1) = {1};
Physical Volume(2) = {2};
Physical Volume(3) = {3};

xmin() = Surface In BoundingBox{ -0.1, -1, -1,    0.1, 201, 121};
xmax() = Surface In BoundingBox{199.9, -1, -1,  200.1, 201, 121};
ymin() = Surface In BoundingBox{ -1,  -0.1, -1,  201,   0.1, 121};
ymax() = Surface In BoundingBox{ -1, 199.9, -1,  201, 200.1, 121};
top()  = Surface In BoundingBox{ -1, -1, 119.9,  201, 201, 120.1};
bot()  = Surface In BoundingBox{ -1, -1,  -0.1,  201, 201,   0.1};

Physical Surface(1) = {xmin()};
Physical Surface(2) = {xmax()};
Physical Surface(3) = {ymin()};
Physical Surface(4) = {ymax()};
Physical Surface(5) = {top()};
Physical Surface(6) = {bot()};

Mesh.CharacteristicLengthMin = 16;
Mesh.CharacteristicLengthMax = 24;
