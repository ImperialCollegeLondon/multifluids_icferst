
s=0.2;

d=0.5;

Point(1) = {0  , 0  , 0, s};
Point(2) = {1  , 0  , 0, s/3};
Point(3) = {1  , 0.5, 0, s/6};
Point(4) = {1.5, 0.5, 0, s/6};
Point(5) = {1.5, 0  , 0, s/3};
Point(6) = {5.0, 0  , 0, s};

Point(7) = {0  , 2-d  , 0, s};
Point(8) = {1  , 2-d  , 0, s/3};
Point(9) = {1  , 1.5-d, 0, s/6};
Point(10)= {1.5, 1.5-d, 0, s/6};
Point(11)= {1.5, 2-d  , 0, s/3};
Point(12)= {5.0, 2-d  , 0, s};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 12};
Line(7) = {12, 11};
Line(8) = {11, 10};
Line(9) = {10, 9};
Line(10) = {9, 8};
Line(11) = {8, 7};
Line(12) = {7, 1};

Line Loop(13) = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
Plane Surface(14) = {13};




Physical Line(1) = {12}; // left
Physical Line(2) = {1, 2, 3, 4, 5, 7, 8, 9, 10, 11}; // top-bottom
Physical Line(3) = {6}; // right
Physical Surface(4) = {14}; // volume
