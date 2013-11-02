#ifndef ETCHING_MENU
#define ETCHING_MENU
#include <QWidget>

class QLineEdit;
class QPushButton;
class QComboBox;
class EtchingMenu : public QWidget
{
	Q_OBJECT
	public:
		EtchingMenu (QWidget *parent=0);
	signals:
		void startEtching(int, int iterations, float* rates);
		
	private:
		QComboBox *simTypeComboBox;
		QLineEdit *iterationsLine;
		QLineEdit *f1NbLine;
		QLineEdit *f2NbLine;
		QLineEdit *f3NbLine;
		QLineEdit *f4NbLine;
		QLineEdit *simTypeLine;
		QPushButton *etchButton;
		float removRates[3];
	private slots:
		void startEtching();
};
#endif
