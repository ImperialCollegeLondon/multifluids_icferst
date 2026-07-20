// Three-layer box for the saturation conservation-across-adapt test.
//
// Generate with (gmsh v2 format):
//   gmsh -3 -format msh2 -o three_layer_box.msh three_layer_box.geo

SetFactory("OpenCASCADE");

Box(1) = {0, 0,  0, 200, 200, 60};   // bottom
Box(2) = {0, 0, 60, 200, 200, 20};   // barrier
Box(3) = {0, 0, 80, 200, 200, 40};   // top

Coherence;

Physical Volume(1) = {1};
Physical Volume(2) = {2};
Physical Volume(3) = {3};

xmin() = Surface In BoundingBox{ -0.1, -1, -1,    0.1, 201, 121};
xmax() = Surface In BoundingBox{199.9, -1, -1,  200.1, 201, 121};
ymin() = Surface In BoundingBox{ -1,  -0.1, -1,  201,   0.1, 121};
ymax() = Surface In BoundingBox{ -1, 199.9, -1,  201, 200.1, 121};
top()  = Surface In BoundingBox{ -1, -1, 119.9,  201, 201, 120.1};
bot()  = Surface In BoundingBox{ -1, -1,  -0.1,  201, 201,   0.1};

Physical Surface(1) = {xmin()};  // side x = 0
Physical Surface(2) = {xmax()};  // side x = 200
Physical Surface(3) = {ymin()};  // side y = 0
Physical Surface(4) = {ymax()};  // side y = 200
Physical Surface(5) = {top()};   // top    z = 120
Physical Surface(6) = {bot()};   // bottom z = 0

Mesh.CharacteristicLengthMin = 20;
Mesh.CharacteristicLengthMax = 30;
