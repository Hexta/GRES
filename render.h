#include <QWidget>
#include <QGLWidget>
#include <QMouseEvent>
#include <QKeyEvent>
#include <QColorDialog>
#include <QVector>
#include <vector>
#include <list>
#include <cmath>
#include "functions.h"
#include "consts.h"
#include "dataTypes.h"
#ifdef _WIN32
#include <windows.h>
#define GL_GLEXT_PROTOTYPES
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glext.h>
#endif
using namespace std;

class SelectAtomMenu;

class Render : public QGLWidget
{
    Q_OBJECT
public:
    Render(QWidget *parent = 0);
     QSize minimumSizeHint() const;
     QSize sizeHint() const;
	public slots:
		void view(surface3D* surface, vector<atomType>* ,atomsCoords*, float* Xsize,float* Ysize,float* Zsize, int z_min, int z_center, int width, int height, coords3D* Vx, coords3D* Vy, coords3D* Vz, vizType vizualType);
		void changeVizType(vizType* type);
		void saveResult();
	signals:
		void etching(void);
protected:
    void initializeGL();
    void resizeGL(int width, int height);
    void paintGL();
    void mousePressEvent(QMouseEvent *event);
    void mouseMoveEvent(QMouseEvent *event);
	void mouseReleaseEvent(QMouseEvent *event);
    void keyPressEvent ( QKeyEvent * event );
private:
    void draw();
	GLfloat drotationX;
	GLfloat drotationY;
	GLfloat drotationZ;
	GLfloat movX;
	GLfloat sumMovX;
	GLfloat sumMovY;
	GLfloat movY;
    GLfloat scale;
    QPoint lastPos;
	SelectAtomMenu* selAtomMenu;

	struct SelAtomType
	{
		short xC;
		short yC;
		short zC;
		float x;
		float y;
		float z;
		char type;
		char fNbCount;
	}selAtomType;
	struct TypeAndCoordsOfAtom
	{
		float x;
		float y;
		float z;
		char fNbCount;
	};
	vector<TypeAndCoordsOfAtom> typeAndCoordsOfAtoms;
	AtomsNames atNames;
	float sR;
	void sphereTemplate(float);
	void cylinderTemplate(float);
	coords3D Vx, Vy, Vz;
	GLuint theSphere;
	QColor clearColor;
	atomsCoords cellAtoms;
	surface3D surfaceXYZ;
	float xs, ys, zs;
	int z_center, z_min;
	int SIZE_X, SIZE_Y;
	bool scribble;
	bool moving;
	bool rotating;
	bool dataChanged;
	float scaling;
	cellAtom coordsOfAtoms; //список координат атомов
	Bonds bonds;
	QAction *exitAction;
	void createActions();
	void processSelection(int x, int y);
	void processSelectionMenu();
	void processAtom (GLuint *pSelectBuff);
	bool dataInitialized;
	vizType visualType;//тип визуализации
	cells surfPoints;
	atomsCoords surfVertex;
	atomsCoords surfNormals;
	vector <GLuint> buffers;
	int sphereQual;
	int vSize1, vSize2, vSize3;
	vector<GLfloat> matrix;
	void createAtomsAndBonds(surface3D* surface,atomsCoords*, float Xsize,float Ysize,float Zsize, int z_min, AtomsNames*, Bonds*);
	void createSurfacePoints(surface3D* surface, float Xsize,float Ysize,float Zsize, int z_min, atomsCoords* surfV, atomsCoords* surfN);
	void initMatrix(vector<GLfloat>*);
	void drawAxis();
	void setGeometry(GLfloat zCenter=0);
	vector<atomType> surfAtoms;

#ifdef _WIN32
	PFNGLBINDBUFFERARBPROC pglBindBufferARB;
	PFNGLDELETEBUFFERSARBPROC pglDeleteBuffersARB;
	PFNGLBUFFERDATAARBPROC pglBufferDataARB;
	PFNGLBUFFERSUBDATAARBPROC pglBufferSubDataARB;
#endif
};
