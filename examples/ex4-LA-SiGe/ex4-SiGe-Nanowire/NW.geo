// Gmsh project
SetFactory("OpenCASCADE");

// Added substrate
Box(1) = {0, 0, 23.666, 27.15, 27.15, 19766.334};
Physical Volume(1)={1};

// Substrate
Box(2) = {0, 0, 0, 27.15, 27.15, 23.666};

// Nano-wire
Cylinder(3) = {13.575, 13.575, -10, 0, 0, 10, 6, 2*Pi};
Physical Volume(2)={2,3};

Box(4) = {0, 0, -10, 27.15, 27.15, 33.666};
v() = BooleanFragments{ Physical Volume{2}; Delete;}{ Volume{4}; Delete; };
Physical Volume(3)={4};

//+  Air
Box(0) = {0, 0, -210, 27.15, 27.15, 20000};
v() = BooleanFragments{ Physical Volume{1:3}; Delete;}{ Volume{0}; Delete; };
Physical Volume(0)={5};

//Meshing Size Air
out[] = PointsOf { Physical Volume{0}; };
Characteristic Length {out[]} = 6.7875 ;

//Meshing AddedSubstrate
out[] = PointsOf { Physical Volume{1}; };
Characteristic Length {out[]} = 27.15 ;

//Meshing Substrate
out[] = PointsOf { Physical Volume{2}; };
Characteristic Length {out[]} = 1.5 ;

//Meshing Oxide
out[] = PointsOf { Physical Volume{3}; };
Characteristic Length {out[]} = 1.5 ;
