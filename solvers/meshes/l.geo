n = 0.0625;

Point(1) = {0, 0, 0, n};
Point(2) = {2, 0, 0, n};
Point(3) = {2, 1, 0, n};
Point(4) = {1, 1, 0, n};
Point(5) = {1, 2, 0, n};
Point(6) = {0, 2, 0, n};
Line(1) = {6, 1};
Line(2) = {1, 2};
Line(3) = {2, 3};
Line(4) = {3, 4};
Line(5) = {4, 5};
Line(6) = {5, 6};
Line Loop(7) = {6, 1, 2, 3, 4, 5};
Plane Surface(8) = {7};
Physical Surface(9) = {8};
