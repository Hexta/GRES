#ifndef GEOMETRY
#    define GEOMETRY

#    include "functions.h"
#    include "dataTypes.h"

void createAtomsAndBondes( surface3D* surface, vector<atomType>*, atomsCoords* cellAts, float xs_, float ys_, float zs_, int z_min, float scaling, vector<atomName>* atNames_, Bonds* outBonds );
void createSphere( GLdouble radius, GLint slices, GLint stacks, int* vSize1, int* vSize2, int* vSize3 );
void normalize( float v[3] );
void normalize( float v[3], coords3D* );
coords3D normalize( coords3D );
void norm( coords3D &in );
coords3D normcrossprod( coords3D in1, coords3D in2 );

#endif
