Point(1) = {0,0,0,1.0};
Point(2) = {1,0,0,1.0};
Point(3) = {0.5,-0.0001,0,1.0};

Line(1) = {1,2};
Circle(2) = {1,3,2};

Line Loop(3) = {1,2};

Plane Surface(10) = {3};