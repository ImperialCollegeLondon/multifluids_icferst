// Box geometry for the priming toy z_primetoy_box_compressible_rhoC.
// 100 m cube, coarse initial mesh, physical surface ids 1..6 and volume id 1.
// Generate the mesh with: gmsh -3 -format msh2 box.geo -o box.msh
lc = 12.0;
Point(1) = {0,   0,   0,   lc};
Point(2) = {100, 0,   0,   lc};
Point(3) = {100, 100, 0,   lc};
Point(4) = {0,   100, 0,   lc};
Point(5) = {0,   0,   100, lc};
Point(6) = {100, 0,   100, lc};
Point(7) = {100, 100, 100, lc};
Point(8) = {0,   100, 100, lc};
Line(1) = {1,2}; Line(2) = {2,3}; Line(3) = {3,4}; Line(4) = {4,1};
Line(5) = {5,6}; Line(6) = {6,7}; Line(7) = {7,8}; Line(8) = {8,5};
Line(9) = {1,5}; Line(10) = {2,6}; Line(11) = {3,7}; Line(12) = {4,8};
Line Loop(1) = {1,2,3,4};      Plane Surface(1) = {1};   // bottom z=0
Line Loop(2) = {5,6,7,8};      Plane Surface(2) = {2};   // top z=100
Line Loop(3) = {1,10,-5,-9};   Plane Surface(3) = {3};   // y=0
Line Loop(4) = {2,11,-6,-10};  Plane Surface(4) = {4};   // x=100
Line Loop(5) = {3,12,-7,-11};  Plane Surface(5) = {5};   // y=100
Line Loop(6) = {4,9,-8,-12};   Plane Surface(6) = {6};   // x=0
Surface Loop(1) = {1,2,3,4,5,6};
Volume(1) = {1};
Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};
Physical Surface(5) = {5};
Physical Surface(6) = {6};
Physical Volume(1) = {1};
