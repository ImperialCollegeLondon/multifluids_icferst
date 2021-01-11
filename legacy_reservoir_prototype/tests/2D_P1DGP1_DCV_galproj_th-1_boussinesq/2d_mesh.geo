// Gmsh project created on Tue Jul  4 10:40:54 2017
nxlayers=22;
nylayers=20;
Point(1) = {0., 0., 0.};
Extrude {2.22, 0, 0} {
  Point{1}; Layers{nxlayers};
}
Extrude {0, 2.0, 0} {
  Line{1}; Layers{nylayers};
}

Transfinite Surface{5}={1,2,3,4} Alternate;

// left
Physical Line(6) = {3};
// right
Physical Line(7) = {4};
// bottom
Physical Line(8) = {1};
// top
Physical Line(9) = {2};
// surface
Physical Surface(10) = {5};
