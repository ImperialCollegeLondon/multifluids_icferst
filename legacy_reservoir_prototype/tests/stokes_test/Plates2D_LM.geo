resolution_low  = 0.2;
resolution_high = 0.05;
resolution_med  = 0.05;

xmax = 6*Pi;
ymax = (2.0/3.0)*Pi;
zmax = 1.0;

//Computational Domain 

Point(1) = {  0.0,  0.0,    0.0, resolution_low};
Point(2) = {  xmax, 0.0,    0.0, resolution_low};
Point(3) = {  xmax, zmax/2, 0.0, resolution_low};
Point(4) = {  xmax, zmax,   0.0, resolution_low};
Point(5) = {  0.0,  zmax,   0.0, resolution_low};
Point(6) = {  0.0,  zmax/2, 0.0, resolution_low};


//Horizontal Lines // 
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Line(7) = {6, 3};

// Line Loop
Line Loop(1) = {1, 2, -7, 6};
Line Loop(2) = {7, 3,  4, 5};

// Surface
Plane Surface(1) = {1};
Plane Surface(2) = {2};

Surface Loop(1) = {1 ,2};

Coherence;

Physical Line(1) = {5,6};   // INLET
Physical Line(3) = {2,3};   // OUTLET
Physical Line(4) = {4};     // TOP
Physical Line(5) = {1};     // BOTTOM

Physical Surface(1) = {1, 2};
