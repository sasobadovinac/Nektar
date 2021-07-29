// Cube containing all element types. Adapted from example from Gmsh
// mailing list: http://geuz.org/pipermail/gmsh/2012/007042.html
//
// Does not work at high order thanks to pyramids.
//
// Run: gmsh -3 CubeAllElements.geo

lc = 0.5;

Point(1) = {-2, -2, -2, lc};
Point(2) = {2, -2, -2, lc};
Point(3) = {2, 2, -2, lc};
Point(4) = {-2, 2, -2, lc};

Line(1) = {1,2} ;
Line(2) = {3,2} ;
Line(3) = {3,4} ;
Line(4) = {4,1} ;

Line Loop(5) = {1,-2,3,4} ;

Plane Surface(6) = {5} ;

Transfinite Line{1} = 5;
Transfinite Line{2} = 5;
Transfinite Line{3} = 5;
Transfinite Line{4} = 5;

Transfinite Surface{6} = {1,2,3,4};

Recombine Surface{6};

extrudetovol1[] = Extrude {0.0, 0.0, 0.8} {Surface {6}; Layers {1}; Recombine;};
surfOpposite = extrudetovol1[0];

extrudetovol2[] = Extrude {0.0, 0.0, 0.8} {Surface {surfOpposite}; Layers{2}; Recombine;};
surfOpposite2 = extrudetovol2[0];

extrudetovol3[] = Extrude {0.0, 0.0, 0.8} {Surface {surfOpposite2}; Layers{3}; QuadTriAddVerts; Recombine;};

//surfOpposite3 = extrudetovol3[0];
Extrude {0.0, 0.0, 0.8} {Surface {6}; Layers{2}; Recombine;};

//surfOpposite4 = extrudetovol4[0];


Extrude {0.0, 0.0, 0.8} {Surface {6}; Layers{1}; Recombine;}

Physical Volume(0) = {1,2,3};
Physical Surface(1) = {6};
Physical Surface(2) = {15, 37, 59};
Physical Surface(3) = {19, 41, 63};
Physical Surface(4) = {23, 45, 67};
Physical Surface(5) = {27, 49, 71};
Physical Surface(6) = {72};
