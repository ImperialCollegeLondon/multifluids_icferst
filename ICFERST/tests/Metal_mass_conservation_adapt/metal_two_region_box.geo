// Two-region box for the metal (solid+fluid) mass conservation-across-adapt test.
//
// Generate with (gmsh v2 format):
//   gmsh -3 -format msh2 -o metal_two_region_box.msh metal_two_region_box.geo

SetFactory("OpenCASCADE");

Box(1) = {-50, -50, -50,  100, 100, 50};   // top half
Box(2) = {-50, -50, -100, 100, 100, 50};   // bottom half

// Make the two boxes conforming and share the interface surface
Coherence;

// Volumes
Physical Volume(1) = {1};
Physical Volume(2) = {2};

// Outer surfaces
top()  = Surface In BoundingBox{-51, -51,   -0.1, 51, 51,   0.1};
bot()  = Surface In BoundingBox{-51, -51, -100.1, 51, 51, -99.9};
xmin() = Surface In BoundingBox{-50.1, -51, -101, -49.9, 51, 1};
xmax() = Surface In BoundingBox{ 49.9, -51, -101,  50.1, 51, 1};
ymin() = Surface In BoundingBox{-51, -50.1, -101, 51, -49.9, 1};
ymax() = Surface In BoundingBox{-51,  49.9, -101, 51,  50.1, 1};

Physical Surface(1) = {top()};   // top     (z =    0)
Physical Surface(2) = {bot()};   // bottom  (z = -100)
Physical Surface(3) = {xmin()};  // side x = -50
Physical Surface(4) = {xmax()};  // side x = +50
Physical Surface(5) = {ymin()};  // side y = -50
Physical Surface(6) = {ymax()};  // side y = +50

Mesh.CharacteristicLengthMin = 8;
Mesh.CharacteristicLengthMax = 15;
