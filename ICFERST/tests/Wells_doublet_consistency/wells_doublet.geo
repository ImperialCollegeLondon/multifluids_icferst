// 3-layer box with a well doublet. Regions: 1 overburden, 2 aquifer, 3 underburden,
// well sleeves 11/12 (wells 1/2 in region 1) and 21/22 (in region 2).
// Surfaces: 1 bottom, 2 top, 3 y=-200, 4 x=300, 5 y=200, 6 x=-300, 101/102 well heads.
SetFactory("OpenCASCADE");
S = 4; LC = 60; LW = 5;
Box(1) = {-300, -200, -150, 600, 400, 150};
Box(2) = {-300, -200, -200, 600, 400, 50};
Box(3) = {-300, -200, -250, 600, 400, 50};
Box(4) = {-100-S, -S, -150, 2*S, 2*S, 150};
Box(5) = { 100-S, -S, -150, 2*S, 2*S, 150};
Box(6) = {-100-S, -S, -200, 2*S, 2*S, 50};
Box(7) = { 100-S, -S, -200, 2*S, 2*S, 50};
Point(101) = {-100, 0, 0}; Point(102) = {-100, 0, -200};
Point(103) = { 100, 0, 0}; Point(104) = { 100, 0, -200};
Line(101) = {101, 102}; Line(102) = {103, 104};
BooleanFragments{ Volume{1:7}; Curve{101, 102}; Delete; }{}
e = 1e-3;
Physical Volume(11) = Volume In BoundingBox{-100-S-e, -S-e, -150-e, -100+S+e, S+e, e};
Physical Volume(12) = Volume In BoundingBox{ 100-S-e, -S-e, -150-e,  100+S+e, S+e, e};
Physical Volume(21) = Volume In BoundingBox{-100-S-e, -S-e, -200-e, -100+S+e, S+e, -150+e};
Physical Volume(22) = Volume In BoundingBox{ 100-S-e, -S-e, -200-e,  100+S+e, S+e, -150+e};
v1() = Volume In BoundingBox{-300-e, -200-e, -150-e, 300+e, 200+e, e};
v2() = Volume In BoundingBox{-300-e, -200-e, -200-e, 300+e, 200+e, -150+e};
sl() = Volume In BoundingBox{-100-S-e, -S-e, -200-e, -100+S+e, S+e, e};
sl() += Volume In BoundingBox{ 100-S-e, -S-e, -200-e,  100+S+e, S+e, e};
v1() -= {sl()}; v2() -= {sl()};
Physical Volume(1) = {v1()};
Physical Volume(2) = {v2()};
Physical Volume(3) = Volume In BoundingBox{-300-e, -200-e, -250-e, 300+e, 200+e, -200+e};
Physical Surface(1) = Surface In BoundingBox{-300-e, -200-e, -250-e, 300+e, 200+e, -250+e};
Physical Surface(101) = Surface In BoundingBox{-100-S-e, -S-e, -e, -100+S+e, S+e, e};
Physical Surface(102) = Surface In BoundingBox{ 100-S-e, -S-e, -e,  100+S+e, S+e, e};
top() = Surface In BoundingBox{-300-e, -200-e, -e, 300+e, 200+e, e};
top() -= Surface In BoundingBox{-100-S-e, -S-e, -e, -100+S+e, S+e, e};
top() -= Surface In BoundingBox{ 100-S-e, -S-e, -e,  100+S+e, S+e, e};
Physical Surface(2) = {top()};
Physical Surface(3) = Surface In BoundingBox{-300-e, -200-e, -250-e, 300+e, -200+e, e};
Physical Surface(4) = Surface In BoundingBox{300-e, -200-e, -250-e, 300+e, 200+e, e};
Physical Surface(5) = Surface In BoundingBox{-300-e, 200-e, -250-e, 300+e, 200+e, e};
Physical Surface(6) = Surface In BoundingBox{-300-e, -200-e, -250-e, -300+e, 200+e, e};
Field[1] = Distance; Field[1].CurvesList = {Curve In BoundingBox{-100-e, -e, -200-e, -100+e, e, e}, Curve In BoundingBox{100-e, -e, -200-e, 100+e, e, e}}; Field[1].Sampling = 100;
Field[2] = Threshold; Field[2].InField = 1; Field[2].SizeMin = LW; Field[2].SizeMax = LC; Field[2].DistMin = S; Field[2].DistMax = 120;
Background Field = 2;
Mesh.MeshSizeExtendFromBoundary = 0; Mesh.MeshSizeFromPoints = 0; Mesh.MeshSizeFromCurvature = 0;
Mesh.MshFileVersion = 2.2;
