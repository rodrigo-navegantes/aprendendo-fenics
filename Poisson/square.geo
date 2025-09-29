// square.geo
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {0, 1, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};

// Physical groups
Physical Surface("domain", 1) = {6};
Physical Line("right", 2) = {2};  // boundary x=1
Physical Line("top",   3) = {3};
Physical Line("left",  4) = {4};  // boundary x=0
Physical Line("bottom",5) = {1};//+
Transfinite Curve {4, 1, 2, 3} = 11 Using Progression 1;
//+
Transfinite Surface {6} Alternated;
