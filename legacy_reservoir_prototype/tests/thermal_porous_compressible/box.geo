l=82;
w=95;
h=25;
nlx=2;
nly=2;
nlz=1;

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

Physical Surface(1) = {16};         // bottom
Physical Surface(2) = {5,20,28,29}; //sides
Physical Surface(3) = {24};         // top
Physical Volume (4) = {1};          // volume


