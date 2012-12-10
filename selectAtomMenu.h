#ifndef SELECT_ATOM_MENU
#    define SELECT_ATOM_MENU

#    include <QWidget>

class QLineEdit;
class QLabel;

class SelectAtomMenu : public QWidget
{
    Q_OBJECT
public:
    SelectAtomMenu ( QWidget *parent = 0 );
    QSize sizeHint( ) const;
public slots:
    void setInfo( int x, int y, int z, int type, int fNbCount );
private:
    QLabel *atomPosition;
    QLabel *firstNeighboursCount;
} ;
#endif