R=1;
r_0=0.1;



Point(1) = { 0, 0, 0.};
Point(2) = { R, 0, 0.};
Point(3) = { 0, R, 0.};
Point(4) = {-R, 0, 0.};
Point(5) = {0 ,-R, 0.};
Point(6) = {r_0, 0, 0.};
Point(7) = {0, r_0, 0.};
Point(8) = {-r_0,0, 0.};
Point(9) = {0, -r_0, 0.};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};
Circle(5) = {6,1,7};
Circle(6) = {7,1,8};
Circle(7) = {8,1,9};
Circle(8) = {9,1,6};


Line(11)={2,6};
Line(12)={3,7};
Line(13)={4,8};
Line(14)={5,9};

Line Loop(5) = {1,12,-5,-11};
Line Loop(6) = {2,13,-6,-12};
Line Loop(7) = {3,14,-7,-13};
Line Loop(8) = {4,11,-8,-14};


Plane Surface(6) = {5};
Plane Surface(7) = {6};
Plane Surface(8) = {7};
Plane Surface(9) = {8};

//Transfinite Line {1,2,3,4}=8;
//Transfinite Line {5,6,7,8}=12;

Field[1]=Attractor;
Field[1].NodesList={1};

Field[2]=MathEval;
Field[2].F="0.1+0.1*(F1-0.1)";


Background Field=2;


Physical Line(7) = {5,6,7,8};//inner
Physical Line(8) = {1,2,3,4};//outer

Physical Surface(11) = {6};
Physical Surface(12) = {7};
Physical Surface(13) = {8};
Physical Surface(14) = {9};
