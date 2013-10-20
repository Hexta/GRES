#define _STDCALL_SUPPORTED
#define GLUT_NO_LIB_PRAGMA

#include <QApplication>
#include "render.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include "face.h"

using namespace std;

int
main(int argc, char *argv[]) {
    Q_INIT_RESOURCE(resources);
    QApplication app(argc, argv);
#ifdef _WIN32
    QGLFormat fmt;
    fmt.setAlpha(true);
    fmt.setDoubleBuffer(false);
    QGLFormat::setDefaultFormat(fmt);
#endif
    MainW *content;
    content = new MainW(0, argc, argv);
    content->setWindowTitle(QObject::tr("GRES"));
    content->show();
    return app.exec();
}
