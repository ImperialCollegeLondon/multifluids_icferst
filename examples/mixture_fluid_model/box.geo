l=5.0;
w=0.25;
h=1.25;

s=10;
nlx=4*s;
nly=1;
nlz=1*s;

Point(1) = {0,0,0,1};

Extrude {0,w,0} {
  Point{1};Layers{nly};
}

Extrude {0,0,h} {
  Line{1};Layers{nlz};
}

Line Loop(6) = {2,-4,-1,3};

Plane Surface(7) = {6};

Extrude {l,0,0} {
  Surface{5};Layers{nlx};
}

Physical Surface(3) = {16};         // bottom
Physical Surface(1) = {5,20,28,29}; //sides
Physical Surface(4) = {24};         // top
Physical Volume (2) = {1};          // volume
