cl__1 = 1;
Point(1) = {-0.5, 0, 0, 1};
Point(2) = {0.5, 0, 0, 1};
Point(3) = {0.5, 0.1, 0, 1};
Point(5) = {-0.5, 0.1, 0, 1};
Line(1) = {1, 2};
Transfinite Line {1} = 50Using Progression 1;
Line(2) = {2, 3};
Transfinite Line {2} = 10Using Progression 1;
Line(3) = {3, 5};
Transfinite Line {3} = 50Using Progression 1;
Line(4) = {5, 1};
Transfinite Line {4} = 10Using Progression 1;
Line Loop(6) = {1, 2, 3, 4};
Plane Surface(6) = {6};
/*Transfinite Surface {6};*/
/*Recombine Surface {6};*/
Physical Point("forcepoint") = {3};
Physical Line("bottom") = {1};
Physical Line("rightside") = {2};
Physical Line("top") = {3};
Physical Line("leftside") = {4};
Physical Surface("Interior") = {6};
