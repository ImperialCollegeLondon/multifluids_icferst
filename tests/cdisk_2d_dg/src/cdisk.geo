Point (1) = {-1, -1, 0, 0.05};
Point (2) = {1, -1, 0, 0.05};
Point (3) = {1, 1, 0, 0.05};
Point (4) = {-1, 1, 0, 0.05};
Line (1) = {4, 1};
Line (2) = {1, 2};
Line (3) = {2, 3};
Line (4) = {3, 4};
Line Loop (1) = {1, 2, 3, 4};
Plane Surface (1) = {1};

// Volume number for whole domain.
Physical Surface (1) = {1};
// Top of the box.
Physical Line(1) = {4};
// Rest of the walls.
Physical Line(2) = {1,2,3};
