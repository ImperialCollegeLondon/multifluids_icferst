// Gmsh project created on Tue Jul  4 10:40:54 2017
nxlayers=11;
nzlayers=10;
nylayers=1;
Point(1) = {0., 0., 0.};
Extrude {2.22, 0, 0} {
  Point{1}; Layers{nxlayers};
}
Extrude {0, 0, 2.0} {
  Line{1}; Layers{nzlayers};
}
Extrude {0, 0.5, 0} {
  Surface{5}; Layers{nylayers};
}

Physical Volume(101) = {1};

// left
Physical Surface(6) = {26};
// bottom 
Physical Surface(8) = {14};
// top
Physical Surface(9) = {22};
// right
Physical Surface(7) = {18};
// side 1
Physical Surface(11) = {5};
// side 2
Physical Surface(12) = {27};




