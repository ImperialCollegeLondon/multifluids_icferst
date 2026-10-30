// ATES doublet. Regions: 1 overburden, 2 aquifer, 3 fracture zone; well sleeves 101-103 (well 1) and 201-203 (well 2).
// Surfaces: 1 bottom, 2 top, 3 y=-150, 4 x=200, 5 y=150, 6 x=-200, 22 well 1 head, 21 well 2 head.
SetFactory("OpenCASCADE");
S = 3; LC = 50; LW = 4; e = 1e-3;
Box(1) = {-200, -150, -84, 400, 300, 84};
Box(2) = {-200, -150, -110, 400, 300, 26};
Box(3) = {-200, -150, -125, 400, 300, 15};
Box(4) = {-200, -150, -200, 400, 300, 75};
Box(5) = {-30-S, -S, -160, 2*S, 2*S, 160};
Box(6) = { 30-S, -S, -160, 2*S, 2*S, 160};
Point(101) = {-30, 0, 0}; Point(102) = {-30, 0, -160};
Point(103) = { 30, 0, 0}; Point(104) = { 30, 0, -160};
Line(101) = {101, 102}; Line(102) = {103, 104};
BooleanFragments{ Volume{1:6}; Curve{101, 102}; Delete; }{}
For w In {0:1}
  xw = -30 + 60*w;
  Physical Volume(101 + 100*w) = Volume In BoundingBox{xw-S-e, -S-e, -84-e, xw+S+e, S+e, e};
  Physical Volume(102 + 100*w) = {Volume In BoundingBox{xw-S-e, -S-e, -110-e, xw+S+e, S+e, -84+e},
                                  Volume In BoundingBox{xw-S-e, -S-e, -160-e, xw+S+e, S+e, -125+e}};
  Physical Volume(103 + 100*w) = Volume In BoundingBox{xw-S-e, -S-e, -125-e, xw+S+e, S+e, -110+e};
EndFor
s1() = Volume In BoundingBox{-30-S-e, -S-e, -200-e, -30+S+e, S+e, e};
s2() = Volume In BoundingBox{ 30-S-e, -S-e, -200-e,  30+S+e, S+e, e};
r1() = Volume In BoundingBox{-200-e, -150-e, -84-e, 200+e, 150+e, e};
r2() = Volume In BoundingBox{-200-e, -150-e, -110-e, 200+e, 150+e, -84+e};
r2() += Volume In BoundingBox{-200-e, -150-e, -200-e, 200+e, 150+e, -125+e};
r3() = Volume In BoundingBox{-200-e, -150-e, -125-e, 200+e, 150+e, -110+e};
r1() -= {s1(), s2()}; r2() -= {s1(), s2()}; r3() -= {s1(), s2()};
Physical Volume(1) = {r1()};
Physical Volume(2) = {r2()};
Physical Volume(3) = {r3()};
h1() = Surface In BoundingBox{-30-S-e, -S-e, -e, -30+S+e, S+e, e};
h2() = Surface In BoundingBox{ 30-S-e, -S-e, -e,  30+S+e, S+e, e};
top() = Surface In BoundingBox{-200-e, -150-e, -e, 200+e, 150+e, e};
top() -= {h1(), h2()};
Physical Surface(1) = Surface In BoundingBox{-200-e, -150-e, -200-e, 200+e, 150+e, -200+e};
Physical Surface(2) = {top()};
Physical Surface(3) = Surface In BoundingBox{-200-e, -150-e, -200-e, 200+e, -150+e, e};
Physical Surface(4) = Surface In BoundingBox{200-e, -150-e, -200-e, 200+e, 150+e, e};
Physical Surface(5) = Surface In BoundingBox{-200-e, 150-e, -200-e, 200+e, 150+e, e};
Physical Surface(6) = Surface In BoundingBox{-200-e, -150-e, -200-e, -200+e, 150+e, e};
Physical Surface(22) = {h1()};
Physical Surface(21) = {h2()};
Field[1] = Distance; Field[1].CurvesList = {Curve In BoundingBox{-30-e, -e, -160-e, -30+e, e, e}, Curve In BoundingBox{30-e, -e, -160-e, 30+e, e, e}}; Field[1].Sampling = 100;
Field[2] = Threshold; Field[2].InField = 1; Field[2].SizeMin = LW; Field[2].SizeMax = LC; Field[2].DistMin = S; Field[2].DistMax = 100;
Background Field = 2;
Mesh.MeshSizeExtendFromBoundary = 0; Mesh.MeshSizeFromPoints = 0; Mesh.MeshSizeFromCurvature = 0;
Mesh.MshFileVersion = 2.2;
