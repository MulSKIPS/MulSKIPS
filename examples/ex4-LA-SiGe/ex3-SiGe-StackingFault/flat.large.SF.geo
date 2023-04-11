// Gmsh project
SetFactory("OpenCASCADE");

// Added substrate
Box(2) = {0, 0, 20.951, 27.15, 27.15, 19769.049};
Physical Volume(2)={2};

// Substrate
Box(1) = {0, 0, -10, 27.15, 27.15, 30.951};
Physical Volume(1)={1};

//+  Air
Box(0) = {0, 0, -210, 27.15, 27.15, 20000};
v() = BooleanFragments{ Physical Volume{1:2}; Delete;}{ Volume{0}; Delete; };
Physical Volume(0)={3};


//Meshing AddedSubstrate
out[] = PointsOf { Physical Volume{2}; };
Characteristic Length {out[]} = 27.15 ;

//Meshing Size Air
out[] = PointsOf { Physical Volume{0}; };
Characteristic Length {out[]} = 6.7875 ;

//Meshing Substrate
out[] = PointsOf { Physical Volume{1}; };
Characteristic Length {out[]} = 1.5 ;
