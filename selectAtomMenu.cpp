#include "selectAtomMenu.h"
#include <QLabel>
#include <QGridLayout>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QLineEdit>
#include <QGroupBox>

SelectAtomMenu::SelectAtomMenu( QWidget *parent ) : QWidget( parent )
{
	QLabel *firstNeighboursCountLabel = new QLabel( tr( "first" ) );

	firstNeighboursCount = new QLabel( );
	firstNeighboursCount -> setTextInteractionFlags( Qt::TextSelectableByMouse );

	atomPosition = new QLabel( );

	QGroupBox *atomPositionBox = new QGroupBox( tr( "Atom`s position" ) );

	QGridLayout *atomPositionGrid = new QGridLayout;

	atomPositionGrid ->addWidget( atomPosition, 0, 0 );

	QSpacerItem *positionSpace = new QSpacerItem( 1, 1, QSizePolicy::Expanding, QSizePolicy::Minimum );
	atomPositionGrid->addItem( positionSpace, 0, 1 );
	atomPositionBox->setLayout( atomPositionGrid );

	QGridLayout *grid = new QGridLayout;
	grid->addWidget( atomPositionBox, 0, 0 );

	QGroupBox *neighboursBox = new QGroupBox( tr( "Atom`s neighbours" ) );
	QGridLayout *neighbourGrid = new QGridLayout;
	neighbourGrid -> addWidget( firstNeighboursCountLabel, 0, 0 );
	neighbourGrid -> addWidget( firstNeighboursCount, 0, 1 );
	QSpacerItem *neighbSpace = new QSpacerItem( 1, 1, QSizePolicy::Expanding, QSizePolicy::Minimum );
	neighbourGrid -> addItem( neighbSpace, 0, 2 );
	neighboursBox->setLayout( neighbourGrid );

	grid->addWidget( neighboursBox, 1, 0 );

	QSpacerItem *space = new QSpacerItem( 10, 10, QSizePolicy::Expanding, QSizePolicy::Expanding );
	grid->addItem( space, 2, 0 );

	setLayout( grid );
}

void
SelectAtomMenu::setInfo( int x, int y, int z, int type, int fNbCount )
{
	QString position;
	position = "(" + QString::number( x ) + "," + QString::number( y ) + "," + QString::number( z ) +
		"," + QString::number( type ) + ")";
	atomPosition -> setText( position );
	firstNeighboursCount->setText( QString::number( fNbCount ) );
}

QSize
SelectAtomMenu::sizeHint( ) const
{
	return QSize( 100, 300 );
}
