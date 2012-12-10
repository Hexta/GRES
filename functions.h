#ifndef FUNCTIONS
#    define FUNCTIONS

#    include <cmath>
#    include <vector>
#    include <algorithm>
#    include <cstdlib>

using namespace std;

struct atomType
{
    int x; //сдвиг ячейки по OX (-1;0;+1)
    int y;
    int z;
    char type; //тип атома (1 -- 8)
    bool toDel;
} ;
typedef vector<atomType> soseds; //соседи одного атома
typedef vector<soseds> allSoseds; //соседи всех атомов

struct atomInfo
{
    soseds neighbours;
    char fNbCount;
    bool deleted;
} ;

struct coords3D
{
    float x;
    float y;
    float z;
} ;

struct length
{
    float x1;
    float y1;
    float z1;
    float x2;
    float y2;
    float z2;
    float length;
} ;

struct rectangle
{
    float x1;
    float y1;
    float z1;
    float x2;
    float y2;
    float z2;
    float x3;
    float y3;
    float z3;
    float x4;
    float y4;
    float z4;
    float S;
} ;
typedef vector<coords3D> cellAtom;
typedef vector<coords3D> atomsCoords;
typedef vector<atomsCoords> cells;
typedef vector<length>lengthes;
typedef vector<rectangle>rectangles;
typedef vector<atomInfo> cell;
typedef vector<cell> surface1D;
typedef vector<surface1D> surface2D;
typedef vector<surface2D> surface3D;
void findSoseds ( allSoseds*, atomsCoords*, float* xs, float* ys, float* zs );
bool rect_comp( const rectangle  &r1, const rectangle &r2 );
bool cells_comp( const atomsCoords  &c1, const atomsCoords &c2 );
double distance( const double& x1, const double& y1, const double& z1, const double& x2, const double& y2, const double& z2 );
double distance( coords3D V ); //Вычисляет длину вектора
coords3D pointShifting( coords3D P1, coords3D P2 ); //считает сдвиг между точками (вектор)
coords3D pointShift( coords3D A, coords3D V ); //сдвигает точку на вектор
bool compareTranslCell( atomsCoords cellAtoms, atomsCoords* allAtoms, coords3D V ); //Совпадение атомов после трансляции
void coordsMove( atomsCoords *ca, coords3D O, coords3D Vx, coords3D Vy, coords3D Vz ); // смещение координат в ячейку
void cellOptim( atomsCoords *ca, float*, float*, float* ); //чистка от общих атомов
double ScalarMult( coords3D V1, coords3D V2 );
void recallNeighbours( surface3D* surface, vector<atomType> *surfAtoms, int x, int y, int z, int type );
bool selAtom( surface3D* surface, vector<atomType> *surfAtoms, allSoseds*, int z_min, atomsCoords *tA, vector<bool>* mask, float* rates );
bool selAtomCA( surface3D* surface, vector<atomType> *surfAtoms, int z_min, atomsCoords *tA, vector<bool>* mask, float* rates );

atomsCoords findCell ( int h, int k, int l, float* xs, float* ys, float* zs, coords3D* Vx, coords3D* Vy, coords3D* Vz ); //Поиск эл ячейки
void addLayer( surface3D* surface, allSoseds*, int sX, int sY, int sZ );
double VectorQuad( coords3D V );
bool coords3Dcompare( coords3D, coords3D );
atomsCoords atomsInBox( atomsCoords* atoms, coords3D Vx, coords3D Vy, coords3D Vz, coords3D P1 );
void findZmin( surface3D* surface, int* z_min );
void optimizeSurface( surface3D* surface, int z_min );
void delAtom( surface3D* surface, vector<atomType> *surfAtoms, int x, int y, int z, int type, int surfAtN );
bool operator==(const  atomType &a1, const  atomType &a2 );
coords3D operator +(const coords3D& v1, const coords3D& v2 );
coords3D operator *(const int& n, const coords3D& v );

#endif
