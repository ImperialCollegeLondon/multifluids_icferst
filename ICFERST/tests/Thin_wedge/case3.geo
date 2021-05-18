H_X=3.0;
H_Y=0.3;
H_Z=0.0;

s=0.02;

Point(1)={0,-0.1,0,s};
Point(2)={0,H_Y,0,s};
Point(3)={H_X,H_Y,H_Z,s};
Point(4)={H_X,-0.1,H_Z,s};

Point(5)={0,0.25*H_Y,0,s};
Point(6)={0,0.3*H_Y,0,s};
Point(7)={0.75*H_X,0.25*H_Y,0,s};

Line(1)={1,4};
Line(2)={4,3};
Line(3)={3,2};
Line(4)={2,6};
Line(5)={6,7};
Line(6)={7,5};
Line(7)={5,1};
Line(8)={5,6};

Line Loop(1)={1,2,3,4,5,6,7};
Line Loop(2)={5,6,8};

Plane Surface(1)={1};
Plane Surface(2)={2};

Physical Surface(2)={1};
Physical Surface(1)={2};

Physical Line(6)={4,7,8}; // inflow
Physical Line(7)={2}; // outflow
Physical Line(8)={1,3}; // sides

Mesh.Algorithm = 4;
Mesh.CharacteristicLengthMin = 0.1;
Mesh.CharacteristicLengthMax = 40;
