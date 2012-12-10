#ifndef SETTINGS_DOCK
#    define SETTINGS_DOCK
#    include <QWidget>

class QGroupBox;
class QSpinBox;

class Settings : public QWidget
{
    Q_OBJECT
public:
    Settings ( QWidget *parent = 0 );
public slots:
    void set ( int* h, int* k , int* l, int xS, int yS, int zS );
signals:
    void settingsChanged( int*, int*, int*, int*, int*, int* );
private:
    QSpinBox*	hSpinBox;
    QSpinBox*	kSpinBox;
    QSpinBox*	lSpinBox;
    QSpinBox*	xSpinBox;
    QSpinBox*	ySpinBox;
    QSpinBox*	zSpinBox;

private slots:
    void get( );

} ;
#endif
