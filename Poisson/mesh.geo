// Gmsh project created on Sun Sep 28 11:06:04 2025
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Transfinite Curve {1} = 11 Using Progression 1;
//+
Physical Point("Left", 2) = {1};
//+
Physical Point("Right", 3) = {2};
//+
Physical Curve("Domain", 4) = {1};
