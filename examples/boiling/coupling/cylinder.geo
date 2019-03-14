n=2;
l=0.4;
r=0.15;
nr=5*n;
nl=5*n;
Point(1) = {0,0,0,0};

Point(2) = {r,0,0,0};
Point(3) = {0,r,0,0};
Point(4) = {-r,0,0,0};
Point(5) = {0,-r,0,0};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

Line Loop(5) = {1,2,3,4};
Plane Surface(6) = {5};

Extrude {0,0,l} {
  Surface{6};Layers{nl};
}
Physical Surface(1) = {6};         // bottom
Physical Surface(2) = {15,19,23,27}; //sides
Physical Surface(3) = {28};         // top
Physical Volume (4) = {1};          // volume
