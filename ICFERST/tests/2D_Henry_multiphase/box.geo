Point(1) = {0, 0, 0, 1.0};
Point(2) = {2, 0, 0, 1.0};
Point(3) = {2, -1, 0, 1.0};
Point(4) = {0, -1, 0, 1.0};
Line(1) = {1, 4};
Line(2) = {4, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
Physical Line(7) = {1};
Physical Line(8) = {3};
Physical Line(9) = {2};
Physical Line(10) = {4};
Physical Surface(11) = {6};


Mesh.Algorithm = 1;
Mesh.MshFileVersion = 2;
Mesh.RecombinationAlgorithm = 0;
Mesh.CharacteristicLengthMin = 0.01;
Mesh.CharacteristicLengthMax = 0.05;
