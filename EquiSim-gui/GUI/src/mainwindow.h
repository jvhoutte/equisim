#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QSlider>
#include <QLCDNumber>
#include <QRadioButton>
#include <QButtonGroup>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkCompositePolyDataMapper2.h>
#include <vtkActor.h>
#include <vtkScalarBarActor.h>
#include <vtkRenderer.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkDICOMImageReader.h>
#include <vtkImageViewer2.h>
#include <vtkImageSlice.h>
#include <vtkImageSliceMapper.h>
#include <vtkWindowToImageFilter.h>
#include <vtkSphereSource.h>
#include <vtkImageMapToColors.h>
#include <vtkLookupTable.h>
#include <vtkImageActor.h>
#include <vtkSTLReader.h>
#include <vtkCompositeDataDisplayAttributes.h>
#include <vtkRendererCollection.h>
#include <vtkPNGWriter.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <random>
#include <vtkTIFFWriter.h>

#include "SKELETON/vtkStatModelReader.h"
#include "SKELETON/DeformRigidSurface.h"
#include "SKELETON/vtkBlockColorSelector.h"
#include "SKELETON/vtkSkeletonVisualiser.h"

#include "GUI/CorrelationFileReader.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    
private slots:
    void on_btnOpenFile_clicked();
    void on_btnSave_clicked();
    void on_btnScreenshot_clicked();
    
    void on_slider_angle1_sliderMoved();
    void on_slider_angle2_sliderMoved();
    void on_slider_angle3_sliderMoved();
    
    void on_slider_pc1_sliderMoved();
    void on_slider_pc2_sliderMoved();
    void on_slider_pc3_sliderMoved();
    void on_slider_pc4_sliderMoved();
    void on_slider_pc5_sliderMoved();
    void on_slider_pc6_sliderMoved();
    void on_slider_pc7_sliderMoved();
    void on_slider_pc8_sliderMoved();
    void on_slider_pc9_sliderMoved();
    void on_slider_pc10_sliderMoved();
    void on_slider_pc11_sliderMoved();
    void on_slider_pc12_sliderMoved();
    void on_slider_pc13_sliderMoved();
    void on_slider_pc14_sliderMoved();
    void on_slider_pc15_sliderMoved();
    void on_slider_pc16_sliderMoved();
    void on_slider_pc17_sliderMoved();
    void on_slider_pc18_sliderMoved();
    void on_slider_pc19_sliderMoved();
    void on_slider_pc20_sliderMoved();
    void on_slider_pc21_sliderMoved();
    void on_slider_pc22_sliderMoved();
    void on_slider_pc23_sliderMoved();
    void on_slider_pc24_sliderMoved();
    void on_slider_pc25_sliderMoved();
    void on_slider_pc26_sliderMoved();
    void on_slider_pc27_sliderMoved();
    void on_slider_pc28_sliderMoved();
    void on_slider_pc29_sliderMoved();
    void on_slider_pc30_sliderMoved();
    void on_slider_pc31_sliderMoved();
    void on_slider_pc32_sliderMoved();
    void on_slider_pc33_sliderMoved();
    void on_slider_pc34_sliderMoved();
    void on_slider_pc35_sliderMoved();
    void on_slider_pc36_sliderMoved();
    void on_slider_pc37_sliderMoved();
    void on_slider_pc38_sliderMoved();
    void on_slider_pc39_sliderMoved();
    void on_slider_pc40_sliderMoved();
    void on_slider_pc41_sliderMoved();
    void on_slider_pc42_sliderMoved();
    void on_slider_pc43_sliderMoved();
    void on_slider_pc44_sliderMoved();
    void on_slider_pc45_sliderMoved();
    void on_slider_pc46_sliderMoved();
    void on_slider_pc47_sliderMoved();
    void on_slider_pc48_sliderMoved();
    void on_slider_pc49_sliderMoved();
    void on_slider_pc50_sliderMoved();
    void on_slider_pc51_sliderMoved();
    void on_slider_pc52_sliderMoved();
    void on_slider_pc53_sliderMoved();
    void on_slider_pc54_sliderMoved();
    void on_slider_pc55_sliderMoved();
    void on_slider_pc56_sliderMoved();
    void on_slider_pc57_sliderMoved();
    void on_slider_pc58_sliderMoved();
    void on_slider_pc59_sliderMoved();
    void on_slider_pc60_sliderMoved();
    void on_slider_pc61_sliderMoved();
    void on_slider_pc62_sliderMoved();
    void on_slider_pc63_sliderMoved();
    void on_slider_pc64_sliderMoved();
    void on_slider_pc65_sliderMoved();
    void on_slider_pc66_sliderMoved();
    void on_slider_pc67_sliderMoved();
    void on_slider_pc68_sliderMoved();
    
    void on_biometric_slider_sliderMoved();
    
    void on_radioButton_1_clicked();
    void on_radioButton_2_clicked();
    void on_radioButton_3_clicked();
    void on_radioButton_4_clicked();
    void on_radioButton_5_clicked();
    void on_radioButton_6_clicked();
    void on_radioButton_7_clicked();
    void on_radioButton_8_clicked();
    void on_radioButton_9_clicked();
    void on_radioButton_10_clicked();
    void on_radioButton_11_clicked();
    void on_radioButton_12_clicked();
    void on_radioButton_13_clicked();
    void on_radioButton_14_clicked();
    void on_radioButton_15_clicked();
    void on_radioButton_16_clicked();
    void on_radioButton_17_clicked();
    void on_radioButton_18_clicked();
    void on_radioButton_19_clicked();
    void on_radioButton_20_clicked();
    void on_radioButton_21_clicked();
    void on_radioButton_22_clicked();
    void on_radioButton_23_clicked();
    void on_radioButton_24_clicked();
    void on_radioButton_25_clicked();
    void on_radioButton_26_clicked();
    void on_radioButton_27_clicked();
    void on_radioButton_28_clicked();
    void on_radioButton_29_clicked();
    void on_radioButton_30_clicked();
    void on_radioButton_31_clicked();
    void on_radioButton_32_clicked();
    void on_radioButton_33_clicked();
    void on_radioButton_34_clicked();
    void on_radioButton_35_clicked();
    void on_radioButton_36_clicked();
    void on_radioButton_37_clicked();
    void on_radioButton_38_clicked();
    void on_radioButton_39_clicked();
    void on_radioButton_40_clicked();
    void on_radioButton_41_clicked();
    void on_radioButton_42_clicked();
    void on_radioButton_43_clicked();
    void on_radioButton_44_clicked();
    void on_radioButton_45_clicked();
    void on_radioButton_46_clicked();
    void on_radioButton_47_clicked();
    void on_radioButton_48_clicked();
    void on_radioButton_49_clicked();
    void on_radioButton_50_clicked();
    void on_radioButton_51_clicked();
    void on_radioButton_52_clicked();
    void on_radioButton_53_clicked();
    void on_radioButton_54_clicked();
    void on_radioButton_55_clicked();
    void on_radioButton_56_clicked();
    void on_radioButton_57_clicked();
    void on_radioButton_58_clicked();
    void on_radioButton_59_clicked();
    void on_radioButton_60_clicked();
    void on_radioButton_61_clicked();
    void on_radioButton_62_clicked();
    void on_radioButton_63_clicked();
    void on_radioButton_64_clicked();
    void on_radioButton_65_clicked();
    void on_radioButton_66_clicked();
    void on_radioButton_67_clicked();
    void on_radioButton_68_clicked();
    
private:
    Ui::MainWindow *ui;
    
    vtkSmartPointer<vtkStatModelReader> statmodelreader;
    vtkSmartPointer<vtkRigidSurfaceTransformFilter> skinning;
    vtkSmartPointer<vtkSkeleton> skeleton;
        
    void BoneAngleSlider(QSlider* slider, int boneid, bool update=true);
    void PcSlider(int pc, bool update=true);
    void UpdateStatModelReader();
    void RenderModel(bool keepcamera=true);
    void DisplayAngles();
    void CheckButton();
    void ColorModel();
    
    std::vector<double> angles;
    std::vector<double> pcweights;
    vtkSmartPointer<vtkDoubleArray> anglerange; 
    vtkSmartPointer<vtkDoubleArray> boneanglearray; 
    vtkSmartPointer<vtkRenderer> renderer;
    vtkSmartPointer<vtkBlockColorSelector> colorselector;
    vtkSmartPointer<vtkSkeletonVisualiser> skeletonvis;
        
    QRadioButton *nocolorbutton;
    QButtonGroup *buttonGroup;
    int lastcheckedbutton;
    std::vector<QSlider*> pcslidervec;
    
    void ChangeBiometric();
    BiometricCorrelations::BiometricCorrMapType biomap;
    
    vtkSmartPointer<vtkMultiBlockDataSet> SplitIntoBlocks(bool linstructure = false);
    vtkSmartPointer<vtkStringArray> GetFieldNameArray(int segmidx);
    vtkSmartPointer<vtkPolyData> AddLandmarkArray(vtkSmartPointer<vtkPolyData> poly);
    
    int NUMBER_PCWEIGHTS = 68;
    
private slots:
    void indexChanged(int index);
    
};

#endif // MAINWINDOW_H
