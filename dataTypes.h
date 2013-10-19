#ifndef DATA_TYPES
#define DATA_TYPES
#include <vector>
using namespace std;

struct atomName
{
	int name;
	int xC;
	int yC;
	int zC;
	float x;
	float y;
	float z;
	unsigned char type;
	int fNbCount;
};
	typedef vector<atomName>AtomsNames;
	
	struct Bond
{
	float x1; float y1; float z1;
	float x2; float y2; float z2;
};
	typedef vector<Bond> Bonds;
#endif