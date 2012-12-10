#include "functions.h"
#include <QTime>

float P1 = 1.0;
float P2 = 0.030;
float P3 = 0.00010;

const int NUMBER_OF_ATOMS_IN_CELL = 8;

coords3D atomTypes[] = {
	{1.0, 1.0, 1.0 },
	{1.0, 0.5, 0.5 },
	{0.5, 1.0, 0.5 },
	{0.5, 0.5, 1.0 },
	{0.25, 0.25, 0.25 },
	{0.75, 0.75, 0.25 },
	{0.75, 0.25, 0.75 },
	{0.25, 0.75, 0.75 } };

void
addLayer( surface3D* surface, allSoseds* sosedi, int sX, int sY, int sZ )
{
	surface2D surfaceXY;
	for ( int y = 0; y < sY; ++y )
	{
		surface1D surfaceX;
		for ( int x = 0; x < sX; ++x )
		{
			cell cell;
			for ( unsigned int a = 0; a < sosedi->size( ); ++a )
			{
				soseds neighbs;
				char numberNeighbs = 0; //Число первых соседей
				for ( int nb = 0; nb < 4; ++nb )
				{
					if ( x + ( *sosedi )[a][nb].x >= 0 && y + ( *sosedi )[a][nb].y >= 0 && x + ( *sosedi )[a][nb].x < sX && y + ( *sosedi )[a][nb].y < sY )
					{
						++numberNeighbs;
						atomType neighb = { x + ( *sosedi )[a][nb].x, y + ( *sosedi )[a][nb].y, sZ + ( *sosedi )[a][nb].z, ( *sosedi )[a][nb].type, false };
						neighbs.push_back( neighb );
					}
				}
				atomInfo atom = { neighbs, numberNeighbs, !numberNeighbs };
				cell.push_back( atom );
			}
			surfaceX.push_back( cell );
		}
		surfaceXY.push_back( surfaceX );
	}
	surface->push_back( surfaceXY );
}

atomsCoords
atomsInBox( atomsCoords* atoms, coords3D Vx, coords3D Vy, coords3D Vz, coords3D P1 )
/*
определяет атомы лежащие внутри прямоугольного параллелипипеда
 */
{
	atomsCoords cellAts;
	for ( unsigned int v = 0; v < atoms->size( ); ++v )
	{
		const float x = ( *atoms )[v].x;
		const float y = ( *atoms )[v].y;
		const float z = ( *atoms )[v].z;
		coords3D X = { x, y, z };
		coords3D V = pointShifting( P1, X );

		double k = ScalarMult( V, Vz ) / VectorQuad( Vz );
		if ( k >= 0.0 && k <= 1.0 )
		{
			k = ScalarMult( V, Vy ) / VectorQuad( Vy );
			if ( k >= 0.0 && k <= 1.0 )
			{
				k = ScalarMult( V, Vx ) / VectorQuad( Vx );
				if ( k >= 0.0 && k <= 1.0 )
					cellAts.push_back( X ); //Записываем атомы в ячейке
			}
		}
	}
	return cellAts;
}

void
cellOptim( atomsCoords *ca, float* xs, float* ys, float* zs )
{
	atomsCoords newCell;
	for ( unsigned int i = 0; i < ca->size( ) - 4; ++i )
	{
		const float x = ( *ca )[i].x;
		const float y = ( *ca )[i].y;
		const float z = ( *ca )[i].z;

		if ( x > 0.01 && y > 0.01 && z > 0.01 )
		{
			coords3D atom = { x, y, z };
			newCell.push_back( atom );
		}
	}
	*xs = distance( ( *ca )[ ca->size( ) - 3 ] );
	*ys = distance( ( *ca )[ ca->size( ) - 2 ] );
	*zs = distance( ( *ca )[ ca->size( ) - 1 ] );
	*ca = newCell;
}

bool
cells_comp( const atomsCoords &c1, const atomsCoords &c2 )
{
	return c1.size( ) < c2.size( );
}

bool
compareTranslCell( atomsCoords cellAtoms, atomsCoords* allAtoms, coords3D V )
{
	bool happy = false;
	for ( unsigned int j = 0; j < cellAtoms.size( ) - 4; ++j )
	{
		happy = false;
		coords3D A0 = cellAtoms[j]; //исходный атом
		coords3D At = pointShift( A0, V );
		for ( unsigned int v = 0; v < allAtoms->size( ); ++v )
		{
			if ( coords3Dcompare( At, ( *allAtoms )[v] ) )
			{
				happy = true;
				break;
			}
		}
		if ( !happy )
			break;
	}
	return happy;
}

bool
coords3Dcompare( coords3D coords1, coords3D coords2 )
{
	//сравнение координат
	return coords1.x == coords2.x && coords1.y == coords2.y && coords1.z == coords2.z;
}

void
coordsMove( atomsCoords *ca, coords3D O, coords3D Vx, coords3D Vy, coords3D Vz )
{
	for ( unsigned int i = 0; i < ca->size( ) - 4; ++i )
	{
		const double x = ScalarMult( pointShifting( O, ( *ca )[i] ), Vx ) / distance( Vx );
		const double y = ScalarMult( pointShifting( O, ( *ca )[i] ), Vy ) / distance( Vy );
		const double z = ScalarMult( pointShifting( O, ( *ca )[i] ), Vz ) / distance( Vz );
		( *ca )[i].x = x;
		( *ca )[i].y = y;
		( *ca )[i].z = z;
	}
}

double
distance( const double& x1, const double& y1, const double& z1, const double& x2, const double& y2, const double& z2 )
{
	return sqrt( pow( ( x2 - x1 ), 2 ) + pow( ( y2 - y1 ), 2 ) + pow( ( z2 - z1 ), 2 ) );
}

double
distance( coords3D V )
{
	return sqrt( pow( V.x, 2 ) + pow( V.y, 2 ) + pow( V.z, 2 ) );
}

coords3D
pointShift( coords3D A, coords3D V )
{
	//сдвигает точку A на вектор V
	coords3D result = { A.x + V.x, A.y + V.y, A.z + V.z };
	return result;
}

coords3D
pointShifting( coords3D P1, coords3D P2 )
{
	//Разность координат двух точек
	coords3D result = { P2.x - P1.x, P2.y - P1.y, P2.z - P1.z };
	return result;
}

atomsCoords
findCell( int h, int k, int l, float* Xsize, float* Ysize, float* Zsize, coords3D* vX, coords3D* vY, coords3D* vZ )
{
	const int SIZE_X = 5;
	const int SIZE_Y = 5;
	const int SIZE_Z = 5;
	atomsCoords allAtoms, atomsP1; // атомы, лежащие на плоскостях
	lengthes ls;
	rectangles rectanglesP1; //прямоугольники
	//создадим кристалл из нескольких ячеек
	for ( int z = 0; z < 9; ++z )
		for ( int y = 0; y < 9; ++y )
			for ( int x = 0; x < 9; ++x )
				for ( int a = 0; a < NUMBER_OF_ATOMS_IN_CELL; ++a )
				{
					coords3D atom = { x + atomTypes[a].x, y + atomTypes[a].y, z + atomTypes[a].z };
					allAtoms.push_back( atom );
				}
	//Найдем свободный член в уравнении секущей плоскости hx+ky+lz-C=0
	int C = 0.5 * ( h * ( SIZE_X - 1 ) + k * ( SIZE_Y - 1 ) + l * ( SIZE_Z - 1 ) ) + 1;
	//Найдем атомы, лежащие на плоскости №1

	atomsP1.reserve( 10000 );
	for ( int z = 0; z < SIZE_Z; ++z )
		for ( int y = 0; y < SIZE_Y; ++y )
			for ( int x = 0; x < SIZE_X; ++x )
				for ( int a = 0; a < NUMBER_OF_ATOMS_IN_CELL; ++a )
					if ( h * ( x + atomTypes[a].x ) + k * ( y + atomTypes[a].y ) + l * ( z + atomTypes[a].z ) - C == 0.0 )
					{
						coords3D atom = { x + atomTypes[a].x, y + atomTypes[a].y, z + atomTypes[a].z };
						atomsP1.push_back( atom );
					}

	//Найдем прямоугольники, лежащие на плоскости №1
	ls.reserve( 50000 );
	for ( unsigned int i = 0; i < atomsP1.size( ); ++i )
		for ( unsigned int j = i + 1; j < atomsP1.size( ); ++j )
		{
			float l1 = distance( atomsP1[i].x, atomsP1[i].y, atomsP1[i].z, atomsP1[j].x, atomsP1[j].y, atomsP1[j].z );
			length L = { atomsP1[i].x, atomsP1[i].y, atomsP1[i].z, atomsP1[j].x, atomsP1[j].y, atomsP1[j].z, l1 };
			ls.push_back( L );
		}
	float kx = -1;
	float ky = -10;
	float kz = -100;

	for ( unsigned int i = 0; i < ls.size( ); ++i )
		for ( unsigned int j = i + 1; j < ls.size( ); ++j )
		{
			float x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, l1;
			x1 = ls[i].x1;
			y1 = ls[i].y1;
			z1 = ls[i].z1;
			x2 = ls[i].x2;
			y2 = ls[i].y2;
			z2 = ls[i].z2;
			x3 = ls[j].x1;
			y3 = ls[j].y1;
			z3 = ls[j].z1;
			x4 = ls[j].x2;
			y4 = ls[j].y2;
			z4 = ls[j].z2;
			l1 = ls[i].length;
			if ( fabs( ls[i].length - ls[j].length ) < 0.01 )
			{
				if ( fabs( x4 - x3 ) > 0.01 )
					kx = ( x2 - x1 ) / ( x4 - x3 );
				else if ( fabs( x2 - x1 ) < 0.01 )
					kx = 1;

				if ( fabs( y4 - y3 ) > 0.01 )
					ky = ( y2 - y1 ) / ( y4 - y3 );
				else if ( fabs( y2 - y1 ) < 0.01 )
					ky = 1;

				if ( fabs( z4 - z3 ) > 0.01 )
					kz = ( z2 - z1 ) / ( z4 - z3 );
				else if ( fabs( z2 - z1 ) < 0.01 )
					kz = 1;
				if ( ( kx == ky && kx == kz ) &&
					 !( ( x2 - x1 )*( x4 - x2 ) + ( y2 - y1 )*( y4 - y2 ) + ( z2 - z1 )*( z4 - z2 ) ) )
				{
					rectangle rect = { x1, y1, z1, x2, y2, z2,
									  x3, y3, z3, x4, y4, z4, ls.at( i ).length };
					rectanglesP1.push_back( rect );

					//удаляем дубли
					for ( unsigned int v = j + 1; v < ls.size( ); ++v )
						if ( ( ls[v].x1 == x1 && ls[v].x2 == x3 && ls[v].y1 == y1 && ls[v].y2 == y3 && ls[v].z1 == z1 && ls[v].z2 == z3 ) ||
							 ( ls[v].x1 == x2 && ls[v].x2 == x4 && ls[v].y1 == y2 && ls[v].y2 == y4 && ls[v].z1 == z2 && ls[v].z2 == z4 ) )
						{
							ls.erase( ls.begin( ) + v );
							--v;
						}
				}
			}
		}
	stable_sort( rectanglesP1.begin( ), rectanglesP1.end( ), rect_comp );

	cells allCells;

	for ( unsigned int i = 0; i < rectanglesP1.size( ); ++i )
	{
		float x1 = rectanglesP1[i].x1, y1 = rectanglesP1[i].y1, z1 = rectanglesP1[i].z1;
		float x2 = rectanglesP1[i].x2, y2 = rectanglesP1[i].y2, z2 = rectanglesP1[i].z2;
		float x3 = rectanglesP1[i].x3, y3 = rectanglesP1[i].y3, z3 = rectanglesP1[i].z3;
		float x4 = rectanglesP1[i].x4, y4 = rectanglesP1[i].y4, z4 = rectanglesP1[i].z4;

		for ( float n = 0.5; n < 5; n += 0.5 )
		{
			int atoms = 0;

			for ( unsigned int j = 0; j < allAtoms.size( ); ++j )
			{
				if ( ( allAtoms[j].x == x1 + n * h && allAtoms[j].y == y1 + n * k && allAtoms[j].z == z1 + n * l )
					 || ( allAtoms[j].x == x2 + n * h && allAtoms[j].y == y2 + n * k && allAtoms[j].z == z2 + n * l )
					 || ( allAtoms[j].x == x3 + n * h && allAtoms[j].y == y3 + n * k && allAtoms[j].z == z3 + n * l )
					 || ( allAtoms[j].x == x4 + n * h && allAtoms[j].y == y4 + n * k && allAtoms[j].z == z4 + n * l ) )
					++atoms;
			}

			if ( atoms == 4 )
			{
				atomsCoords cellAtoms;
				coords3D Vx, Vy, Vz, P1;
				P1.x = x3;
				P1.y = y3;
				P1.z = z3;
				coords3D P2 = { x4, y4, z4 }, P3 = { x1, y1, z1 }, P4 = { x3 + n*h, y3 + n*k, z3 + n * l };
				Vz = pointShifting( P1, P4 );
				Vy = pointShifting( P1, P3 );
				Vx = pointShifting( P1, P2 );
				cellAtoms = atomsInBox( &allAtoms, Vx, Vy, Vz, P1 );

				// а в конец списка атомов запишем векторы координат и координаты начала координат :)
				cellAtoms.push_back( P1 );
				cellAtoms.push_back( Vx );
				cellAtoms.push_back( Vy );
				cellAtoms.push_back( Vz );
				allCells.push_back( cellAtoms );
			}
		}
	}
	stable_sort( allCells.begin( ), allCells.end( ), cells_comp );

	for ( unsigned int i = 0; i < allCells.size( ); ++i )
	{
		bool happy = false;

		//считаем векторы координат
		coords3D P1 = allCells[i][allCells[i].size( ) - 4];
		coords3D Vx = allCells[i][allCells[i].size( ) - 3];
		coords3D Vy = allCells[i][allCells[i].size( ) - 2];
		coords3D Vz = allCells[i][allCells[i].size( ) - 1];
		//транслируем по OX
		happy = compareTranslCell( allCells[i], &allAtoms, Vx );

		//транслируем по OY
		if ( happy )
			happy = compareTranslCell( allCells[i], &allAtoms, Vy );
		//транслируем по OZ
		if ( happy )
			happy = compareTranslCell( allCells[i], &allAtoms, Vz );

		//транслируем по -OZ
		if ( happy )
			happy = compareTranslCell( allCells[i], &allAtoms, -1 * Vz );
		//транслируем по -OY
		if ( happy )
			happy = compareTranslCell( allCells[i], &allAtoms, -1 * Vy );
		//транслируем по -OX
		if ( happy )
			happy = compareTranslCell( allCells[i], &allAtoms, -1 * Vx );

		if ( happy )
		{
			atomsCoords translCell = atomsInBox( &allAtoms, Vx, Vy, Vz, pointShift( P1, Vx ) );
			happy = ( allCells[i].size( ) - 4 == translCell.size( ) );
		}

		if ( happy )
		{
			atomsCoords translCell = atomsInBox( &allAtoms, Vx, Vy, Vz, pointShift( P1, Vy ) );
			happy = ( allCells[i].size( ) - 4 == translCell.size( ) );
		}

		if ( happy )
		{
			atomsCoords translCell = atomsInBox( &allAtoms, Vx, Vy, Vz, pointShift( P1, Vz ) );
			happy = ( allCells[i].size( ) - 4 == translCell.size( ) );
		}

		if ( happy )
		{
			atomsCoords translCell = atomsInBox( &allAtoms, Vx, Vy, Vz, pointShift( P1, -1 * Vx ) );
			happy = ( allCells[i].size( ) - 4 == translCell.size( ) );
		}

		if ( happy )
		{
			atomsCoords translCell = atomsInBox( &allAtoms, Vx, Vy, Vz, pointShift( P1, -1 * Vy ) );
			happy = ( allCells[i].size( ) - 4 == translCell.size( ) );
		}

		if ( happy )
		{
			atomsCoords translCell = atomsInBox( &allAtoms, Vx, Vy, Vz, pointShift( P1, -1 * Vz ) );
			happy = ( allCells[i].size( ) - 4 == translCell.size( ) );
		}
		if ( happy )
		{
			atomsCoords tempCell = allCells[i];
			allCells.clear( );
			allCells.push_back( tempCell );
			*vX = Vx;
			*vY = Vy;
			*vZ = Vz;
			break;
		}
	}

	coordsMove( &( allCells[0] ), allCells[0][allCells[0].size( ) - 4], allCells[0][allCells[0].size( ) - 3], allCells[0][allCells[0].size( ) - 2], allCells[0][allCells[0].size( ) - 1] );

	cellOptim( &( allCells[0] ), Xsize, Ysize, Zsize );

	return allCells[0];
}

void
findSoseds( allSoseds* allSosedi, atomsCoords* atom_Types, float* xs, float* ys, float* zs )
{
	const int NUMBER_OF_ATOMS_IN_CELL = atom_Types->size( );
	allSosedi->clear( );
	for ( int type0 = 0; type0 < NUMBER_OF_ATOMS_IN_CELL; ++type0 )
	{
		soseds sosedi;

		double x0 = ( *atom_Types )[type0].x;
		double y0 = ( *atom_Types )[type0].y;
		double z0 = ( *atom_Types )[type0].z;

		for ( int xc1 = -1; xc1 < 2; ++xc1 )//смещение соседней ячейки №1
			for ( int yc1 = -1; yc1 < 2; ++yc1 )
				for ( int zc1 = -1; zc1 < 2; ++zc1 )
					for ( char type1 = 0; type1 < NUMBER_OF_ATOMS_IN_CELL; ++type1 )
						if ( !( type1 == type0 && !xc1 && !yc1 && !zc1 ) )
						{
							double x1 = ( *atom_Types )[type1].x + xc1**xs;
							double y1 = ( *atom_Types )[type1].y + yc1**ys;
							double z1 = ( *atom_Types )[type1].z + zc1**zs;
							if ( fabs( pow( x1 - x0, 2 ) + pow( y1 - y0, 2 ) + pow( z1 - z0, 2 ) - 3.0 / 16.0 ) <= 0.001 ) //квадрат длины связи #1
							{
								atomType s1 = { xc1, yc1, zc1, type1, false };
								sosedi.push_back( s1 );
							}
						}
		allSosedi->push_back( sosedi );
	}
}

void
findZmin( surface3D* surface, int* zm )
{
	for ( unsigned int z = *zm; z < surface->size( ) - 1; ++z )
		for ( int y = ( *surface )[z].size( ) - 2; --y >= 2; )
			for ( int x = ( *surface )[z][y].size( ) - 2; --x >= 2; )
				for ( int a = ( *surface )[z][y][x].size( ); --a >= 0; )
					if ( !( ( *surface )[z][y][x][a].deleted ) )
					{
						*zm = z;
						return;
					}
}

void
recallNeighbours( surface3D* surface, vector<atomType> *surfAtoms, int x, int y, int z, int type )
{
	for ( unsigned int nb = 0; nb < ( *surface )[z][y][x][type].neighbours.size( ); ++nb )
	{
		const int xNb = ( *surface )[z][y][x][type].neighbours[nb].x;
		const int yNb = ( *surface )[z][y][x][type].neighbours[nb].y;
		const int zNb = ( *surface )[z][y][x][type].neighbours[nb].z;
		const char typeNb = ( *surface )[z][y][x][type].neighbours[nb].type;
		for ( unsigned int i = 0; i < ( *surface )[zNb][yNb][xNb][typeNb].neighbours.size( ); ++i )
			if ( ( *surface )[zNb][yNb][xNb][typeNb].neighbours[i].x == x && ( *surface )[zNb][yNb][xNb][typeNb].neighbours[i].y == y
				 && ( *surface )[zNb][yNb][xNb][typeNb].neighbours[i].z == z
				 && ( *surface )[zNb][yNb][xNb][typeNb].neighbours[i].type == type )
			{
				( *surface )[zNb][yNb][xNb][typeNb].neighbours.erase( ( *surface )[zNb][yNb][xNb][typeNb].neighbours.begin( ) + i );
				( *surface )[zNb][yNb][xNb][typeNb].fNbCount -= 1;
				break;
			}
		if ( !( *surface )[zNb][yNb][xNb][typeNb].fNbCount )
		{
			( *surface )[zNb][yNb][xNb][typeNb].deleted = true;
			( *surface )[zNb][yNb][xNb][typeNb].neighbours.clear( );
		}

		if ( ( ( *surface )[zNb][yNb][xNb][typeNb].fNbCount == 3 ) &&
			 ( ( xNb > 1 && xNb < ( *surface )[zNb][yNb].size( ) - 2 ) && ( yNb > 1 && yNb < ( *surface )[zNb].size( ) - 2 ) ) )
		{
			atomType aT = { xNb, yNb, zNb, typeNb, false };
			surfAtoms->push_back( aT );
		}
	}
}

bool
rect_comp( const rectangle &r1, const rectangle &r2 )
{
	const double S1 = distance( r1.x1, r1.y1, r1.z1, r1.x2, r1.y2, r1.z2 ) *
		distance( r1.x3, r1.y3, r1.z3, r1.x4, r1.y4, r1.z4 );
	const double S2 = distance( r2.x1, r2.y1, r2.z1, r2.x2, r2.y2, r2.z2 ) *
		distance( r2.x3, r2.y3, r2.z3, r2.x4, r2.y4, r2.z4 );
	return S1 < S2;
}

double
ScalarMult( coords3D V1, coords3D V2 )
{
	//Скалярное произведение векторов
	return V1.x * V2.x + V1.y * V2.y + V1.z * V2.z;
}

bool
selAtom( surface3D* surface, vector<atomType> *surfAtoms, allSoseds* sosedi, int z_min, atomsCoords *tA, vector<bool>* mask, float* rates )
{
	P1 = rates[0];
	P2 = rates[1];
	P3 = rates[2];
	const int yC = ( *surface )[z_min].size( );
	const int xC = ( *surface )[z_min][0].size( );
	bool result = false;


	int i = rand( ) * surfAtoms->size( ) / RAND_MAX;
	{
		short x = ( *surfAtoms )[i].x;
		short y = ( *surfAtoms )[i].y;
		unsigned short z = ( *surfAtoms )[i].z;
		unsigned short a = ( *surfAtoms )[i].type;
		if ( ( ( !( *mask ).empty( ) ) && ( ( ( z == 0 ) && !( *mask ) [ ( y - 2 )*( ( *surface )[z][y].size( ) - 4 ) + x - 2 ] )
										 || ( z + ( *tA )[a].z >= 0.5 ) ) )
			 || ( ( *mask ).empty( ) ) )
		{
			float randN = ( float ) rand( ) / RAND_MAX;
			int bonds = ( *surface )[z][y][x][a].fNbCount;
			if ( !bonds )
			{
				surfAtoms->erase( surfAtoms->begin( ) + i );
				--i;
				//continue;
				return false;
			}

			if ( ( bonds == 1 && randN < P1 )\
						  || ( bonds == 2 && randN < P2 )\
						  || ( bonds == 3 && randN < P3 ) )
			{
				delAtom( surface, surfAtoms, x, y, z, a, i );
				/*
				удалим и соостветсвенный краевой атом
				 */
				if ( x == 4 || x == 5 )
				{
					delAtom( surface, surfAtoms, 5 - x, y, z, a, -1 );
				}
				if ( x == xC - 6 || x == xC - 5 )
				{
					delAtom( surface, surfAtoms, 2 * ( xC - 1 ) - 5 - x, y, z, a, -1 );
				}
				if ( y == 4 || y == 5 )
				{
					delAtom( surface, surfAtoms, x, 5 - y, z, a, -1 );
				}
				if ( y == yC - 6 || y == yC - 5 )
				{
					delAtom( surface, surfAtoms, x, 2 * yC - 7 - y, z, a, -1 );
				}
				if ( z >= surface->size( ) - 3 )
				{
					addLayer( surface, sosedi, xC, yC, ( *surface ).size( ) );
				}

			}
		}
	}
	return result;
}

bool
selAtomCA( surface3D* surface, vector<atomType> *surfAtoms, int z_min, atomsCoords *tA, vector<bool>* mask, float* rates )
{
	vector<atomType> atoms2del;
	atoms2del.reserve( surfAtoms->size( ) / 2 );
	P1 = rates[0];
	P2 = rates[1];
	P3 = rates[2];
	const int yC = ( *surface )[z_min].size( );
	const int xC = ( *surface )[z_min][0].size( );
	bool result = false;

	for ( unsigned int i = 0; i < surfAtoms->size( ); ++i )
	{
		short x = ( *surfAtoms )[i].x;
		short y = ( *surfAtoms )[i].y;
		unsigned short z = ( *surfAtoms )[i].z;
		unsigned short a = ( *surfAtoms )[i].type;
		if ( ( ( !( *mask ).empty( ) ) && ( ( ( z + ( *tA )[a].z < 0.5 ) && !( *mask ) [ ( y - 2 )*( ( *surface )[z][y].size( ) - 4 ) + x - 2 ] )
										 || ( z + ( *tA )[a].z >= 0.5 ) ) )
			 || ( ( *mask ).empty( ) ) )
		{
			float randN = ( float ) rand( ) / RAND_MAX;
			int bonds = ( *surface )[z][y][x][a].fNbCount;
			if ( !bonds )
			{
				surfAtoms->erase( surfAtoms->begin( ) + i-- );
				continue;
				return false;
			}

			if ( ( bonds == 1 && randN < P1 )\
						  || ( bonds == 2 && randN < P2 )\
						  || ( bonds == 3 && randN < P3 ) )
			{
				( *surfAtoms )[i].toDel = true;
				if ( z == surface->size( ) - 3 )
					result = true;
			}
		}
	}

	for ( unsigned int i = 0; i < surfAtoms->size( ); ++i )
	{
		if ( ( *surfAtoms )[i].toDel )
		{
			short x = ( *surfAtoms )[i].x;
			short y = ( *surfAtoms )[i].y;
			unsigned short z = ( *surfAtoms )[i].z;
			unsigned short a = ( *surfAtoms )[i].type;
			delAtom( surface, surfAtoms, x, y, z, a, i );

			if ( x == 4 || x == 5 )
			{
				delAtom( surface, surfAtoms, 5 - x, y, z, a, -1 );
			}
			if ( x == xC - 6 || x == xC - 5 )
			{
				delAtom( surface, surfAtoms, 2 * ( xC - 1 ) - 5 - x, y, z, a, -1 );
			}
			if ( y == 4 || y == 5 )
			{
				delAtom( surface, surfAtoms, x, 5 - y, z, a, -1 );
			}
			if ( y == yC - 6 || y == yC - 5 )
			{
				delAtom( surface, surfAtoms, x, 2 * yC - 7 - y, z, a, -1 );
			}
			--i;
		}
	}
	return result;
}

void
delAtom( surface3D* surface, vector<atomType> *surfAtoms, int x, int y, int z, int type, int i )
{
	if ( i != -1 )
		surfAtoms->erase( surfAtoms->begin( ) + i );
	recallNeighbours( surface, surfAtoms, x, y, z, type );
	( *surface )[z][y][x][type].fNbCount = 0;
	( *surface )[z][y][x][type].deleted = true;
	( *surface )[z][y][x][type].neighbours.clear( );
}

double
VectorQuad( coords3D V )
{
	//Возведение вектора в квадрат
	return pow( V.x, 2 ) + pow( V.y, 2 ) + pow( V.z, 2 );
}

coords3D operator +( const coords3D& v1, const coords3D& v2 )
{
	coords3D temp = { v1.x + v2.x, v1.y + v2.y, v1.z + v2.z };
	return temp;
}

coords3D operator *( const int& n, const coords3D& v )
{
	coords3D temp = { n * v.x, n * v.y, n * v.z };
	return temp;
}

void
optimizeSurface( surface3D* surface, int z_min )
{
	if ( z_min > 1 )
		( *surface )[z_min - 2].clear( );
	surface -> swap( *surface );

}

bool operator==( const atomType &a1, const atomType &a2 )
{
	return (( a1.x == a2.x ) && ( a1.y == a2.y ) && ( a1.z == a2.z ) && ( a1.type == a2.type ) );
}
