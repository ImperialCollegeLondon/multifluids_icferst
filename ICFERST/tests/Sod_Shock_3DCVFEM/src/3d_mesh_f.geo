Point(1) = {0., 0., 0.};
Extrude {1, 0, 0} {
  Point{1};
Layers{150};
}
Extrude {0, 0.025, 0} {
  Line{1};
Layers{1};
}
Extrude {0, 0, 0.025} {
  Surface{5};
Layers{1};
}
// Inflow
Physical Surface(1) = {26}; //left
// Sides
Physical Surface(2) = {14, 22}; //y-top
// Top-Bottom
Physical Surface(3) = {5, 27}; //z-top
// Outflow
Physical Surface(4) = {18}; //right
// Volume
Physical Volume(5) = {1};
