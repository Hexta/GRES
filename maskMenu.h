#ifndef MASK_MENU
#define MASK_MENU
#include <QWidget>
#include <vector>

using namespace std;

class QLabel;
class QPushButton;
class QImage;

class MaskMenu : public QWidget {
    Q_OBJECT
public:
    MaskMenu(QWidget *parent = 0, int xS = 10, int yS = 10);
signals:
    void maskChanged(vector<bool>*);
private:

    QLabel* maskPreview;
    QPushButton *loadMaskButton;
    QPushButton *setMaskButton;
    QImage* maskImage;
    int xS, yS;
private slots:
    void loadMask();
    void setMask();
};
#endif
