#ifndef MASK_PREVIEW
#    define MASK_PREVIEW
#    include <QLabel>

class MaskPreview : public QLabel
{
    Q_OBJECT
public:
    MaskPreview ( QWidget *parent = 0 );
} ;
#endif
