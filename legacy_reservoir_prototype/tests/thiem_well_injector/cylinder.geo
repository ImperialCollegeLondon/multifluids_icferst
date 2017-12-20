
s = 200;

Point(1) = {500 , 500 , 0, s/25};
Point(2) = {1000, 500 , 0, s};
Point(3) = {500 , 1000, 0, s};
Point(4) = {500 , 0   , 0, s};
Point(5) = {0   , 500 , 0, s};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 5};
Circle(3) = {5, 1, 4};
Circle(4) = {4, 1, 2};

Line(5) = {1, 2};
Line(6) = {1, 3};
Line(7) = {1, 5};
Line(8) = {1, 4};
Line Loop(9) = {4, -5, 8};
Plane Surface(10) = {9};
Line Loop(11) = {1, -6, 5};
Plane Surface(12) = {11};
Line Loop(13) = {2, -7, 6};
Plane Surface(14) = {13};
Line Loop(15) = {7, 3, -8};
Plane Surface(16) = {15};

Extrude {0, 0, 20} {
  Surface{16, 10, 12, 14}; Layers{8};
}

Physical Surface(1) = {33, 50, 67, 84}; // top
Physical Surface(2) = {16, 10, 12, 14}; // bottom
Physical Surface(3) = {41, 58, 75, 28}; // sides
Physical Volume(4)  = {1, 2, 3, 4};     // volume

Mesh.Algorithm = 4;
Mesh.Remove4Triangles = 1;
Mesh.RecombinationAlgorithm = 0;
Mesh.CharacteristicLengthMin = 1;
Mesh.CharacteristicLengthMax = 30;
