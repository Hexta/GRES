#ifndef SHNJAGA
#define SHNJAGA

#include <QtGui>
#include <QMainWindow>
#    include "functions.h"
#    include "consts.h"

using namespace std;

class Render;
class Settings;
class QActionGroup;
class EtchingMenu;
class MaskMenu;
class MainW : public QMainWindow
{
	Q_OBJECT
	public:
		MainW(QWidget *parent=0 , int argc=NULL, char *const argv[]=NULL);
	public slots:
		void newDocument();
		void getSettings(int*, int*, int *, int*, int *, int*);//считываем значения виджетов
		void changeVizType(QAction* type);
		void etch(int, int* , float* rates); //травление

	private:
		QAction *aboutQtAction;
		QAction *exitAction;
		QAction *newAct;
		QAction *saveAct;
		QAction *viewAsAtsAndBonds_SurfaceAndBulkAct;
		QAction *viewAsAtsAndBonds_SurfaceAct;
		
		QAction *viewAsAts_SurfaceAndBulkAct;
		QAction *viewAsAts_SurfaceAct;
		QAction *viewAsCellsSurface;
		QAction *etchingAction;
		QAction *maskAction;
		
		QMenu *fileMenu;
		QMenu *viewMenu;
		QMenu *viewSurfaceMenu;
		QMenu *viewBulkMenu;
		QMenu *helpMenu;
		QToolBar *fileToolBar;
		QToolBar *etchingToolBar;
		QActionGroup *viewGroup;

		MaskMenu *maskMenu;
		
 		Render* result;
		Settings *settings;
		EtchingMenu* etchMenu;
		void createActions();
		void createMenus();
		void createToolBars();
		void setSettings();//устанавливаем значения виджетов с настройками

		int h,k,l;
		int SIZE_X, SIZE_Y, SIZE_Z;
		atomsCoords cellAtoms;
		surface3D surfaceXYZ;
		float Xsize, Ysize, Zsize;
		vector<atomType> surfAtoms;
		int z_center, z_min;
		coords3D Vx, Vy, Vz;
		bool perfect;
		allSoseds sosedi;
		cells surfacePoints;
		vizType vizualType;//тип визуализации
		void drawResult();
		vector<bool> mask;
	private slots:
		void showEtchMenu();
		void showMenuMask();
		void setMask(vector<bool>*);
};
#endif
