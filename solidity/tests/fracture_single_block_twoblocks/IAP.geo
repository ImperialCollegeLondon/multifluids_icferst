cl1 = 0.5;
a=3;

Point(1) = {-a, -a, -a, cl1};
Point(2) = {a, -a, -a, cl1};
Point(3) = {a, a, -a, cl1};
Point(4) = {-a, a, -a, cl1};
Point(5) = {-a, a, a, cl1};
Point(6) = {-a, -a, a, cl1};
Point(7) = {a, -a, a, cl1};
Point(8) = {a, a, a, cl1};

Line(1) = {7, 2};
Line(2) = {2, 3};
Line(3) = {3, 8};
Line(4) = {7, 8};
Line(5) = {7, 6};
Line(6) = {6, 5};
Line(7) = {5, 8};
Line(8) = {5, 4};
Line(9) = {4, 3};
Line(10) = {1, 2};
Line(11) = {1, 4};
Line(12) = {6, 1};

Line Loop(14) = {5, 12, 10, -1};
Plane Surface(14) = {14};
Line Loop(16) = {4, -3, -2, -1};
Plane Surface(16) = {16};
Line Loop(18) = {10, 2, -9, -11};
Plane Surface(18) = {18};
Line Loop(20) = {11, -8, -6, 12};
Plane Surface(20) = {20};
Line Loop(22) = {8, 9, 3, -7};
Plane Surface(22) = {22};
Line Loop(24) = {5, 6, 7, -4};
Plane Surface(24) = {24};


Surface Loop(26) = {24, 14, 20, 18, 16, 22};

Volume(26) = {26};


Physical Surface(1) = {20};
Physical Surface(2) = {14};
Physical Surface(3) = {18};
Physical Surface(4) = {16};
Physical Surface(5) = {24};
Physical Surface(6) = {22};
Physical Volume(27) = {26};
