#include "mainwindow.h"
#include <QApplication>

#include "QVTKOpenGLNativeWidget.h"
//#include "/home/jeroen/Desktop/ProgramFiles/VTK-9.0.0/GUISupport/Qt/QVTKOpenGLWidget.h"
#include <QSurfaceFormat>

#include "vtkInteractorStyleImage.h"
#include "vtkRenderWindowInteractor.h"

#include "vtkJPEGReader.h"



int main(int argc, char *argv[])
{
    
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    
    return a.exec();
    
}
