// Gmsh project
SetFactory("OpenCASCADE");

// Added substrate
Box(2) = {0, 0, 54.780, 10.860, 10.860, 19736};
Physical Volume(2)={2};

// Substrate
Box(1) = {0, 0, -10, 10.860, 10.860, 64.780};
Physical Volume(1)={1};

//+  Air
Box(0) = {0, 0, -210, 10.860, 10.860, 20000.780};
v() = BooleanFragments{ Physical Volume{1:2}; Delete;}{ Volume{0}; Delete; };
Physical Volume(0)={3};


//Meshing AddedSubstrate
out[] = PointsOf { Physical Volume{2}; };
Characteristic Length {out[]} = 10.860 ;

//Meshing Size Air
out[] = PointsOf { Physical Volume{0}; };
Characteristic Length {out[]} = 10.860 ;

//Meshing Substrate
out[] = PointsOf { Physical Volume{1}; };
Characteristic Length {out[]} = 1 ;
