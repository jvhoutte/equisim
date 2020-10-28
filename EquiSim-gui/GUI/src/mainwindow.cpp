#include "mainwindow.h"
#include "src/ui_mainwindow.h"
#include <QFileDialog>



MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    
    statmodelreader = vtkSmartPointer<vtkStatModelReader>::New();
    skinning = vtkSmartPointer<vtkRigidSurfaceTransformFilter>::New();
    colorselector = vtkSmartPointer<vtkBlockColorSelector>::New();
    
    angles = std::vector<double>(3,0);    
    pcweights = std::vector<double>(NUMBER_PCWEIGHTS,0);
    
    nocolorbutton = new QRadioButton;
    buttonGroup = new QButtonGroup;
    buttonGroup->addButton(ui->radioButton_1); buttonGroup->setId(ui->radioButton_1, 1);
    buttonGroup->addButton(ui->radioButton_2); buttonGroup->setId(ui->radioButton_2, 2);
    buttonGroup->addButton(ui->radioButton_3); buttonGroup->setId(ui->radioButton_3, 3);
    buttonGroup->addButton(ui->radioButton_4); buttonGroup->setId(ui->radioButton_4, 4);
    buttonGroup->addButton(ui->radioButton_5); buttonGroup->setId(ui->radioButton_5, 5);
    buttonGroup->addButton(ui->radioButton_6); buttonGroup->setId(ui->radioButton_6, 6);
    buttonGroup->addButton(ui->radioButton_7); buttonGroup->setId(ui->radioButton_7, 7);
    buttonGroup->addButton(ui->radioButton_8); buttonGroup->setId(ui->radioButton_8, 8);
    buttonGroup->addButton(ui->radioButton_9); buttonGroup->setId(ui->radioButton_9, 9);
    buttonGroup->addButton(ui->radioButton_10); buttonGroup->setId(ui->radioButton_10, 10);
    buttonGroup->addButton(ui->radioButton_11); buttonGroup->setId(ui->radioButton_11, 11);
    buttonGroup->addButton(ui->radioButton_12); buttonGroup->setId(ui->radioButton_12, 12);
    buttonGroup->addButton(ui->radioButton_13); buttonGroup->setId(ui->radioButton_13, 13);
    buttonGroup->addButton(ui->radioButton_14); buttonGroup->setId(ui->radioButton_14, 14);
    buttonGroup->addButton(ui->radioButton_15); buttonGroup->setId(ui->radioButton_15, 15);
    buttonGroup->addButton(ui->radioButton_16); buttonGroup->setId(ui->radioButton_16, 16);
    buttonGroup->addButton(ui->radioButton_17); buttonGroup->setId(ui->radioButton_17, 17);
    buttonGroup->addButton(ui->radioButton_18); buttonGroup->setId(ui->radioButton_18, 18);
    buttonGroup->addButton(ui->radioButton_19); buttonGroup->setId(ui->radioButton_19, 19);
    buttonGroup->addButton(ui->radioButton_20); buttonGroup->setId(ui->radioButton_20, 20);
    buttonGroup->addButton(ui->radioButton_21); buttonGroup->setId(ui->radioButton_21, 21);
    buttonGroup->addButton(ui->radioButton_22); buttonGroup->setId(ui->radioButton_22, 22);
    buttonGroup->addButton(ui->radioButton_23); buttonGroup->setId(ui->radioButton_23, 23);
    buttonGroup->addButton(ui->radioButton_24); buttonGroup->setId(ui->radioButton_24, 24);
    buttonGroup->addButton(ui->radioButton_25); buttonGroup->setId(ui->radioButton_25, 25);
    buttonGroup->addButton(ui->radioButton_26); buttonGroup->setId(ui->radioButton_26, 26);
    buttonGroup->addButton(ui->radioButton_27); buttonGroup->setId(ui->radioButton_27, 27);
    buttonGroup->addButton(ui->radioButton_28); buttonGroup->setId(ui->radioButton_28, 28);
    buttonGroup->addButton(ui->radioButton_29); buttonGroup->setId(ui->radioButton_29, 29);
    buttonGroup->addButton(ui->radioButton_30); buttonGroup->setId(ui->radioButton_30, 30);
    buttonGroup->addButton(ui->radioButton_31); buttonGroup->setId(ui->radioButton_31, 31);
    buttonGroup->addButton(ui->radioButton_32); buttonGroup->setId(ui->radioButton_32, 32);
    buttonGroup->addButton(ui->radioButton_33); buttonGroup->setId(ui->radioButton_33, 33);
    buttonGroup->addButton(ui->radioButton_34); buttonGroup->setId(ui->radioButton_34, 34);
    buttonGroup->addButton(ui->radioButton_35); buttonGroup->setId(ui->radioButton_35, 35);
    buttonGroup->addButton(ui->radioButton_36); buttonGroup->setId(ui->radioButton_36, 36);
    buttonGroup->addButton(ui->radioButton_37); buttonGroup->setId(ui->radioButton_37, 37);
    buttonGroup->addButton(ui->radioButton_38); buttonGroup->setId(ui->radioButton_38, 38);
    buttonGroup->addButton(ui->radioButton_39); buttonGroup->setId(ui->radioButton_39, 39);
    buttonGroup->addButton(ui->radioButton_40); buttonGroup->setId(ui->radioButton_40, 40);
    buttonGroup->addButton(ui->radioButton_41); buttonGroup->setId(ui->radioButton_41, 41);
    buttonGroup->addButton(ui->radioButton_42); buttonGroup->setId(ui->radioButton_42, 42);
    buttonGroup->addButton(ui->radioButton_43); buttonGroup->setId(ui->radioButton_43, 43);
    buttonGroup->addButton(ui->radioButton_44); buttonGroup->setId(ui->radioButton_44, 44);
    buttonGroup->addButton(ui->radioButton_45); buttonGroup->setId(ui->radioButton_45, 45);
    buttonGroup->addButton(ui->radioButton_46); buttonGroup->setId(ui->radioButton_46, 46);
    buttonGroup->addButton(ui->radioButton_47); buttonGroup->setId(ui->radioButton_47, 47);
    buttonGroup->addButton(ui->radioButton_48); buttonGroup->setId(ui->radioButton_48, 48);
    buttonGroup->addButton(ui->radioButton_49); buttonGroup->setId(ui->radioButton_49, 49);
    buttonGroup->addButton(ui->radioButton_50); buttonGroup->setId(ui->radioButton_50, 50);
    buttonGroup->addButton(ui->radioButton_51); buttonGroup->setId(ui->radioButton_51, 51);
    buttonGroup->addButton(ui->radioButton_52); buttonGroup->setId(ui->radioButton_52, 52);
    buttonGroup->addButton(ui->radioButton_53); buttonGroup->setId(ui->radioButton_53, 53);
    buttonGroup->addButton(ui->radioButton_54); buttonGroup->setId(ui->radioButton_54, 54);
    buttonGroup->addButton(ui->radioButton_55); buttonGroup->setId(ui->radioButton_55, 55);
    buttonGroup->addButton(ui->radioButton_56); buttonGroup->setId(ui->radioButton_56, 56);
    buttonGroup->addButton(ui->radioButton_57); buttonGroup->setId(ui->radioButton_57, 57);
    buttonGroup->addButton(ui->radioButton_58); buttonGroup->setId(ui->radioButton_58, 58);
    buttonGroup->addButton(ui->radioButton_59); buttonGroup->setId(ui->radioButton_59, 59);
    buttonGroup->addButton(ui->radioButton_60); buttonGroup->setId(ui->radioButton_60, 60);
    buttonGroup->addButton(ui->radioButton_61); buttonGroup->setId(ui->radioButton_61, 61);
    buttonGroup->addButton(ui->radioButton_62); buttonGroup->setId(ui->radioButton_62, 62);
    buttonGroup->addButton(ui->radioButton_63); buttonGroup->setId(ui->radioButton_63, 63);
    buttonGroup->addButton(ui->radioButton_64); buttonGroup->setId(ui->radioButton_64, 64);
    buttonGroup->addButton(ui->radioButton_65); buttonGroup->setId(ui->radioButton_65, 65);
    buttonGroup->addButton(ui->radioButton_66); buttonGroup->setId(ui->radioButton_66, 66);
    buttonGroup->addButton(ui->radioButton_67); buttonGroup->setId(ui->radioButton_67, 67);
    buttonGroup->addButton(ui->radioButton_68); buttonGroup->setId(ui->radioButton_68, 68);
    buttonGroup->addButton(nocolorbutton);     buttonGroup->setId(nocolorbutton, 69);
    nocolorbutton->setChecked(true);
    buttonGroup->setExclusive(true);
    
    // Combobox for biometrics
    connect(ui->biometric_combobox, SIGNAL(currentIndexChanged(int)), this, SLOT(indexChanged(int)));
    
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_btnOpenFile_clicked(){
    
    std::string skeleton_file = "armature_mean.bvh";
    std::string shape_file = "statmodel_v4.vtm";
    
    biomap = BiometricCorrelations::ReadCorrelationFile("regressioncoeff.txt");
    
    statmodelreader->SetBVHFileName(skeleton_file.c_str());
    statmodelreader->SetVTMFileName(shape_file.c_str());    
    statmodelreader->Update();
    
    skeleton = statmodelreader->GetArmature();
    
    skinning->SetInputData(0, statmodelreader->GetShape()); 
    skinning->SetInputData(1, skeleton);
    skinning->Update();
    
    anglerange = vtkDoubleArray::SafeDownCast(skeleton->GetVertexData()->GetArray("AngleRange"));
    boneanglearray = vtkDoubleArray::SafeDownCast(skeleton->GetVertexData()->GetArray("BoneAngles"));
    
    colorselector->SetInputData(skinning->GetOutput());
    colorselector->SetColorId(-1);
    colorselector->Update();
    
    skeletonvis = vtkSmartPointer<vtkSkeletonVisualiser>::New();
    skeletonvis->SetInputData(skeleton);
    skeletonvis->Update();
    
    // MAPPER FOR GEOMETRY
    
    vtkSmartPointer<vtkCompositePolyDataMapper2> mapper = vtkSmartPointer<vtkCompositePolyDataMapper2>::New();
    mapper->SetInputConnection( colorselector->GetOutputPort() );
    vtkSmartPointer<vtkCompositeDataDisplayAttributes> cdsa = vtkSmartPointer<vtkCompositeDataDisplayAttributes>::New();
    mapper->SetCompositeDataDisplayAttributes(cdsa.Get());
    //double color[] = {1, 1, 1};
    //mapper->SetBlockColor(0,color); // make all blocks white
    
    mapper->SetBlockVisibility(0,true);
    mapper->SetBlockOpacity(0,0.5);
    mapper->ScalarVisibilityOn();
    mapper->SetScalarRange(0,4);
    
    vtkSmartPointer<vtkScalarBarActor> scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
    scalarBar->SetLookupTable(mapper->GetLookupTable());
    scalarBar->SetTitle("PC variance");
    scalarBar->SetNumberOfLabels(3);
    scalarBar->SetBarRatio(0.2);
    scalarBar->SetHeight(5);
    scalarBar->SetMaximumHeightInPixels(400);
    
    // Create a lookup table to share between the mapper and the scalarbar
    vtkSmartPointer<vtkLookupTable> hueLut = vtkSmartPointer<vtkLookupTable>::New();
    hueLut->SetHueRange (0.667, 0.0);
    hueLut->Build();
    mapper->SetLookupTable( hueLut );
    scalarBar->SetLookupTable( hueLut );
    
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
    
    
    // MAPPER FOR SKELETON 
    vtkSmartPointer<vtkPolyDataMapper> polymapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    polymapper->SetInputConnection( skeletonvis->GetOutputPort() );

    vtkSmartPointer<vtkActor> skeletonactor = vtkSmartPointer<vtkActor>::New();
	skeletonactor->SetMapper(polymapper);
    skeletonactor->GetProperty()->SetLineWidth(5);
    
    // SET UP RENDERER
    
    renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->AddActor(actor);
    renderer->AddActor(skeletonactor);
    //renderer->AddActor2D(scalarBar);
    renderer->SetBackground(1,1,1); // Background color green
    
    ui->modelrenderer->renderWindow()->AddRenderer(renderer); 
    
    // Reset sliders 
    ui->slider_angle1->setValue(50); BoneAngleSlider(ui->slider_angle1, 1, false);    
    ui->slider_angle2->setValue(50); BoneAngleSlider(ui->slider_angle2, 2, false);
    ui->slider_angle3->setValue(50); BoneAngleSlider(ui->slider_angle3, 3, false);
    
    pcslidervec = {ui->slider_pc1, ui->slider_pc2, ui->slider_pc3, ui->slider_pc4, ui->slider_pc5, ui->slider_pc6, ui->slider_pc7, ui->slider_pc8, ui->slider_pc9, ui->slider_pc10, ui->slider_pc11, ui->slider_pc12, ui->slider_pc13, ui->slider_pc14, ui->slider_pc15, ui->slider_pc16, ui->slider_pc17, ui->slider_pc18, ui->slider_pc19, ui->slider_pc20, ui->slider_pc21, ui->slider_pc22, ui->slider_pc23, ui->slider_pc24, ui->slider_pc25, ui->slider_pc26, ui->slider_pc27, ui->slider_pc28, ui->slider_pc29, ui->slider_pc30, ui->slider_pc31, ui->slider_pc32, ui->slider_pc33, ui->slider_pc34, ui->slider_pc35, ui->slider_pc36, ui->slider_pc37, ui->slider_pc38, ui->slider_pc39, ui->slider_pc40, ui->slider_pc41, ui->slider_pc42, ui->slider_pc43, ui->slider_pc44, ui->slider_pc45, ui->slider_pc46, ui->slider_pc47, ui->slider_pc48, ui->slider_pc49, ui->slider_pc50, ui->slider_pc51, ui->slider_pc52, ui->slider_pc53, ui->slider_pc54, ui->slider_pc55, ui->slider_pc56, ui->slider_pc57, ui->slider_pc58, ui->slider_pc59, ui->slider_pc60, ui->slider_pc61, ui->slider_pc62, ui->slider_pc63, ui->slider_pc64, ui->slider_pc65, ui->slider_pc66, ui->slider_pc67, ui->slider_pc68};
    
    for(int pc_id=0; pc_id<NUMBER_PCWEIGHTS; pc_id++){
        pcslidervec[pc_id]->setValue(50);
        PcSlider(pc_id, false);
    }
    
    ui->biometric_slider->setValue(50);
        
    RenderModel(false);
}

void MainWindow::on_btnSave_clicked(){
    QString fileName = QFileDialog::getSaveFileName(this,
        tr("Save EquiSim model"), "",
        tr("vtkMultiBlockDataSet (*.vtm)"));
    
    if (fileName.isEmpty())
        return;
    else {
        std::cout << "Writing to " << fileName.toStdString() << std::endl;
        
        vtkSmartPointer<vtkXMLMultiBlockDataWriter> writer = vtkSmartPointer<vtkXMLMultiBlockDataWriter>::New();
        writer->SetInputData(SplitIntoBlocks(true)); //skinning->GetOutput());
        writer->SetFileName(fileName.toStdString().c_str());
        writer->Update();        
    }
    
}

void MainWindow::on_btnScreenshot_clicked(){
    QString fileName = QFileDialog::getSaveFileName(this,
        tr("Save screenshot"), "",
        tr("Image files (*.tif *.png)"));
    
    if (fileName.isEmpty())
        return;
    else {
        std::string filename = fileName.toStdString();
        
        // Capture the display and place in a tiff
        vtkSmartPointer<vtkWindowToImageFilter> w2i = vtkSmartPointer<vtkWindowToImageFilter>::New();
        w2i->SetInput(ui->modelrenderer->renderWindow());
        w2i->Update();
        
        if(filename.substr(filename.find_last_of(".") + 1) == "tif") {
            vtkSmartPointer<vtkTIFFWriter> writer = vtkSmartPointer<vtkTIFFWriter>::New();
            writer->SetInputConnection(w2i->GetOutputPort());
            writer->SetFileName(filename.c_str());
            writer->Write();
        }
        else if(filename.substr(filename.find_last_of(".") + 1) == "png"){
            vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter>::New();
            writer->SetInputConnection(w2i->GetOutputPort());
            writer->SetFileName(filename.c_str());
            writer->Write();
        }
        else{return;}
    }
}

void MainWindow::CheckButton(){
    if( lastcheckedbutton == buttonGroup->checkedId()){nocolorbutton->setChecked(true);}
    ColorModel();
}

void MainWindow::ColorModel(){
    lastcheckedbutton = buttonGroup->checkedId();
    colorselector->SetColorId(lastcheckedbutton);
    if(nocolorbutton->isChecked()){colorselector->SetColorId(-1);}
    RenderModel(true);
}

void MainWindow::on_radioButton_1_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_2_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_3_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_4_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_5_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_6_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_7_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_8_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_9_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_10_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_11_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_12_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_13_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_14_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_15_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_16_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_17_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_18_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_19_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_20_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_21_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_22_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_23_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_24_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_25_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_26_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_27_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_28_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_29_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_30_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_31_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_32_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_33_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_34_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_35_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_36_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_37_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_38_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_39_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_40_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_41_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_42_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_43_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_44_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_45_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_46_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_47_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_48_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_49_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_50_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_51_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_52_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_53_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_54_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_55_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_56_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_57_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_58_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_59_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_60_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_61_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_62_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_63_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_64_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_65_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_66_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_67_clicked(){
    CheckButton();
}

void MainWindow::on_radioButton_68_clicked(){
    CheckButton();
}

void MainWindow::DisplayAngles(){
    
    double boneangle[3]; 
    
    boneanglearray->GetTuple(1, boneangle);
    ui->lcd_lp_a1->display(boneangle[0] * 180/3.1415);
    ui->lcd_lp_a2->display(- boneangle[1] * 180/3.1415);
    ui->lcd_lp_a3->display(boneangle[2] * 180/3.1415);   
    
    boneanglearray->GetTuple(2, boneangle);
    ui->lcd_sp_a1->display(boneangle[0] * 180/3.1415);
    ui->lcd_sp_a2->display(- boneangle[1] * 180/3.1415);
    ui->lcd_sp_a3->display(boneangle[2] * 180/3.1415); 
    
    boneanglearray->GetTuple(3, boneangle);
    ui->lcd_co_a1->display(boneangle[0] * 180/3.1415);
    ui->lcd_co_a2->display(- boneangle[1] * 180/3.1415);
    ui->lcd_co_a3->display(boneangle[2] * 180/3.1415); 
}
    
    
void MainWindow::UpdateStatModelReader(){
    // Update pc weights 
    for(int pc=0; pc<NUMBER_PCWEIGHTS; pc++){
        statmodelreader->SetPcWeights(pc, pcweights[pc]);
    }
    
    statmodelreader->Update();
}

void MainWindow::RenderModel(bool keepcamera){

    // Update bone angles always, also if only PC weights were changed. 
    // Changing PC weights also affects the bone angles slightly,
    // so, they will need to be reset to the chosen angles of the user
    for(int b=1; b<=3; b++){
        double boneangle[3]; boneanglearray->GetTuple(b, boneangle);
        boneangle[1] = -angles[b-1];
        skeleton->SetBoneAnglesById(b,&boneangle);
        skeleton->Modified();
    }
        
    DisplayAngles();
        
    skinning->Update();
    
    // Keep camera direction 
    vtkCamera* camera = NULL;
    if(ui->modelrenderer->renderWindow()->GetRenderers()->GetFirstRenderer() != NULL and keepcamera == true){    
        camera = ui->modelrenderer->renderWindow()->GetRenderers()->GetFirstRenderer()->GetActiveCamera();
        renderer->SetActiveCamera(camera); 
    }
    
    ui->modelrenderer->renderWindow()->Render();
}

void MainWindow::BoneAngleSlider(QSlider* slider, int boneid, bool update)
{
    double min = anglerange->GetComponent(boneid,2);
    double max = anglerange->GetComponent(boneid,3);
    int pos = slider->sliderPosition();
    double angle = min + (max-min)*double(pos)/100.;
    
    angles[boneid-1]=angle;

    if(update==true){RenderModel();}
}

void MainWindow::PcSlider(int pc, bool update)
{
    int pos = pcslidervec[pc]->sliderPosition();
    pcweights[pc] = -1 + 2*double(pos)/100;
    
    if(update==true){
        UpdateStatModelReader();        
        RenderModel();
    }
}

void MainWindow::on_slider_pc1_sliderMoved(){
    ui->radioButton_1->setChecked(true); ColorModel();
    PcSlider(0);    
}

void MainWindow::on_slider_pc2_sliderMoved(){
    ui->radioButton_2->setChecked(true); ColorModel();
    PcSlider(1);    
}

void MainWindow::on_slider_pc3_sliderMoved(){
    ui->radioButton_3->setChecked(true); ColorModel();
    PcSlider(2);    
}

void MainWindow::on_slider_pc4_sliderMoved(){
    ui->radioButton_4->setChecked(true); ColorModel();
    PcSlider(3);    
}

void MainWindow::on_slider_pc5_sliderMoved(){
    ui->radioButton_5->setChecked(true); ColorModel();
    PcSlider(4);    
}

void MainWindow::on_slider_pc6_sliderMoved(){
    ui->radioButton_6->setChecked(true); ColorModel();
    PcSlider(5);    
}

void MainWindow::on_slider_pc7_sliderMoved(){
    ui->radioButton_7->setChecked(true); ColorModel();
    PcSlider(6);    
}

void MainWindow::on_slider_pc8_sliderMoved(){
    ui->radioButton_8->setChecked(true); ColorModel();
    PcSlider(7);    
}

void MainWindow::on_slider_pc9_sliderMoved(){
    ui->radioButton_9->setChecked(true); ColorModel();
    PcSlider(8);     
}

void MainWindow::on_slider_pc10_sliderMoved(){
    ui->radioButton_10->setChecked(true); ColorModel();
    PcSlider(9);    
}

void MainWindow::on_slider_pc11_sliderMoved(){
    ui->radioButton_11->setChecked(true); ColorModel();
    PcSlider(10);    
}

void MainWindow::on_slider_pc12_sliderMoved(){
    ui->radioButton_12->setChecked(true); ColorModel();
    PcSlider(11);    
}

void MainWindow::on_slider_pc13_sliderMoved(){
    ui->radioButton_13->setChecked(true); ColorModel();
    PcSlider(12);    
}

void MainWindow::on_slider_pc14_sliderMoved(){
    ui->radioButton_14->setChecked(true); ColorModel();
    PcSlider(13);    
}

void MainWindow::on_slider_pc15_sliderMoved(){
    ui->radioButton_15->setChecked(true); ColorModel();
    PcSlider(14);    
}

void MainWindow::on_slider_pc16_sliderMoved(){
    ui->radioButton_16->setChecked(true); ColorModel();
    PcSlider(15);    
}

void MainWindow::on_slider_pc17_sliderMoved(){
    ui->radioButton_17->setChecked(true); ColorModel();
    PcSlider(16);    
}

void MainWindow::on_slider_pc18_sliderMoved(){
    ui->radioButton_18->setChecked(true); ColorModel();
    PcSlider(17);    
}

void MainWindow::on_slider_pc19_sliderMoved(){
    ui->radioButton_19->setChecked(true); ColorModel();
    PcSlider(18);     
}

void MainWindow::on_slider_pc20_sliderMoved(){
    ui->radioButton_20->setChecked(true); ColorModel();
    PcSlider(19);    
}

void MainWindow::on_slider_pc21_sliderMoved(){
    ui->radioButton_21->setChecked(true); ColorModel();
    PcSlider(20);    
}

void MainWindow::on_slider_pc22_sliderMoved(){
    ui->radioButton_22->setChecked(true); ColorModel();
    PcSlider(21);    
}

void MainWindow::on_slider_pc23_sliderMoved(){
    ui->radioButton_23->setChecked(true); ColorModel();
    PcSlider(22);    
}

void MainWindow::on_slider_pc24_sliderMoved(){
    ui->radioButton_24->setChecked(true); ColorModel();
    PcSlider(23);    
}

void MainWindow::on_slider_pc25_sliderMoved(){
    ui->radioButton_25->setChecked(true); ColorModel();
    PcSlider(24);    
}

void MainWindow::on_slider_pc26_sliderMoved(){
    ui->radioButton_26->setChecked(true); ColorModel();
    PcSlider(25);    
}

void MainWindow::on_slider_pc27_sliderMoved(){
    ui->radioButton_27->setChecked(true); ColorModel();
    PcSlider(26);    
}

void MainWindow::on_slider_pc28_sliderMoved(){
    ui->radioButton_28->setChecked(true); ColorModel();
    PcSlider(27);    
}

void MainWindow::on_slider_pc29_sliderMoved(){
    ui->radioButton_29->setChecked(true); ColorModel();
    PcSlider(28);    
}

void MainWindow::on_slider_pc30_sliderMoved(){
    ui->radioButton_30->setChecked(true); ColorModel();
    PcSlider(29);     
}

void MainWindow::on_slider_pc31_sliderMoved(){
    ui->radioButton_31->setChecked(true); ColorModel();
    PcSlider(30);    
}

void MainWindow::on_slider_pc32_sliderMoved(){
    ui->radioButton_32->setChecked(true); ColorModel();
    PcSlider(31);    
}

void MainWindow::on_slider_pc33_sliderMoved(){
    ui->radioButton_33->setChecked(true); ColorModel();
    PcSlider(32);    
}

void MainWindow::on_slider_pc34_sliderMoved(){
    ui->radioButton_34->setChecked(true); ColorModel();
    PcSlider(33);    
}

void MainWindow::on_slider_pc35_sliderMoved(){
    ui->radioButton_35->setChecked(true); ColorModel();
    PcSlider(34);    
}

void MainWindow::on_slider_pc36_sliderMoved(){
    ui->radioButton_36->setChecked(true); ColorModel();
    PcSlider(35);    
}

void MainWindow::on_slider_pc37_sliderMoved(){
    ui->radioButton_37->setChecked(true); ColorModel();
    PcSlider(36);    
}

void MainWindow::on_slider_pc38_sliderMoved(){
    ui->radioButton_38->setChecked(true); ColorModel();
    PcSlider(37);    
}

void MainWindow::on_slider_pc39_sliderMoved(){
    ui->radioButton_39->setChecked(true); ColorModel();
    PcSlider(38);    
}

void MainWindow::on_slider_pc40_sliderMoved(){
    ui->radioButton_40->setChecked(true); ColorModel();
    PcSlider(39);     
}

void MainWindow::on_slider_pc41_sliderMoved(){
    ui->radioButton_41->setChecked(true); ColorModel();
    PcSlider(40);    
}

void MainWindow::on_slider_pc42_sliderMoved(){
    ui->radioButton_42->setChecked(true); ColorModel();
    PcSlider(41);    
}

void MainWindow::on_slider_pc43_sliderMoved(){
    ui->radioButton_43->setChecked(true); ColorModel();
    PcSlider(42);    
}

void MainWindow::on_slider_pc44_sliderMoved(){
    ui->radioButton_44->setChecked(true); ColorModel();
    PcSlider(43);    
}

void MainWindow::on_slider_pc45_sliderMoved(){
    ui->radioButton_45->setChecked(true); ColorModel();
    PcSlider(44);    
}

void MainWindow::on_slider_pc46_sliderMoved(){
    ui->radioButton_46->setChecked(true); ColorModel();
    PcSlider(45);    
}

void MainWindow::on_slider_pc47_sliderMoved(){
    ui->radioButton_47->setChecked(true); ColorModel();
    PcSlider(46);    
}

void MainWindow::on_slider_pc48_sliderMoved(){
    ui->radioButton_48->setChecked(true); ColorModel();
    PcSlider(47);    
}

void MainWindow::on_slider_pc49_sliderMoved(){
    ui->radioButton_49->setChecked(true); ColorModel();
    PcSlider(48);    
}

void MainWindow::on_slider_pc50_sliderMoved(){
    ui->radioButton_50->setChecked(true); ColorModel();
    PcSlider(49);    
}

void MainWindow::on_slider_pc51_sliderMoved(){
    ui->radioButton_51->setChecked(true); ColorModel();
    PcSlider(50);     
}

void MainWindow::on_slider_pc52_sliderMoved(){
    ui->radioButton_52->setChecked(true); ColorModel();
    PcSlider(51);    
}

void MainWindow::on_slider_pc53_sliderMoved(){
    ui->radioButton_53->setChecked(true); ColorModel();
    PcSlider(52);    
}

void MainWindow::on_slider_pc54_sliderMoved(){
    ui->radioButton_54->setChecked(true); ColorModel();
    PcSlider(53);    
}

void MainWindow::on_slider_pc55_sliderMoved(){
    ui->radioButton_55->setChecked(true); ColorModel();
    PcSlider(54);    
}

void MainWindow::on_slider_pc56_sliderMoved(){
    ui->radioButton_56->setChecked(true); ColorModel();
    PcSlider(55);    
}

void MainWindow::on_slider_pc57_sliderMoved(){
    ui->radioButton_57->setChecked(true); ColorModel();
    PcSlider(56);    
}

void MainWindow::on_slider_pc58_sliderMoved(){
    ui->radioButton_58->setChecked(true); ColorModel();
    PcSlider(57);    
}

void MainWindow::on_slider_pc59_sliderMoved(){
    ui->radioButton_59->setChecked(true); ColorModel();
    PcSlider(58);    
}

void MainWindow::on_slider_pc60_sliderMoved(){
    ui->radioButton_60->setChecked(true); ColorModel();
    PcSlider(59);    
}

void MainWindow::on_slider_pc61_sliderMoved(){
    ui->radioButton_61->setChecked(true); ColorModel();
    PcSlider(60);     
}

void MainWindow::on_slider_pc62_sliderMoved(){
    ui->radioButton_62->setChecked(true); ColorModel();
    PcSlider(61);    
}

void MainWindow::on_slider_pc63_sliderMoved(){
    ui->radioButton_63->setChecked(true); ColorModel();
    PcSlider(62);    
}

void MainWindow::on_slider_pc64_sliderMoved(){
    ui->radioButton_64->setChecked(true); ColorModel();
    PcSlider(63);    
}

void MainWindow::on_slider_pc65_sliderMoved(){
    ui->radioButton_65->setChecked(true); ColorModel();
    PcSlider(64);    
}

void MainWindow::on_slider_pc66_sliderMoved(){
    ui->radioButton_66->setChecked(true); ColorModel();
    PcSlider(65);    
}

void MainWindow::on_slider_pc67_sliderMoved(){
    ui->radioButton_67->setChecked(true); ColorModel();
    PcSlider(66);    
}

void MainWindow::on_slider_pc68_sliderMoved(){
    ui->radioButton_68->setChecked(true); ColorModel();
    PcSlider(67);    
}






void MainWindow::on_slider_angle1_sliderMoved(){
    BoneAngleSlider(ui->slider_angle1, 1);    
}

void MainWindow::on_slider_angle2_sliderMoved(){
    BoneAngleSlider(ui->slider_angle2, 2);
}

void MainWindow::on_slider_angle3_sliderMoved(){
    BoneAngleSlider(ui->slider_angle3, 3);
}

void MainWindow::indexChanged(int index)
{
    // reset the slider to the mean and update the model for this average metric
    ui->biometric_slider->setValue(50);
    ChangeBiometric(); 
}

void MainWindow::on_biometric_slider_sliderMoved(){
    ChangeBiometric(); 
}

void MainWindow::ChangeBiometric(){
    std::string metricname = ui->biometric_combobox->currentText().toStdString(); 
    
    if ( biomap.find(metricname) == biomap.end() ) { std::cout << "metric name is unknown in the correlation file" << std::endl; return; }
    
    double mean = biomap[metricname]->mean;
    double stdev = 3*biomap[metricname]->stdev;
    
    int pos = ui->biometric_slider->sliderPosition();
    double biometric_value = mean-stdev + 2*stdev*double(pos)/100.;
    
    // Convert rad->degrees and cm->mm for printing the biometric values
    double printbioval = biometric_value;
    if(metricname=="Capsule toeangle" || metricname=="Capsule heelangle" || metricname=="Capsule underrun" || metricname=="Coffin coffinangle" || metricname=="Coffin capsuledev" || metricname=="Coffin palmarangle" ||metricname=="navicular angle"){printbioval *= 180./vtkMath::Pi(); }
    else if(metricname!="Coffin toeheelsupport"){printbioval *= 10;} 
    
    ui->biometric_display->display(printbioval);
    
    // number of pc weights we can change here
    int Npc = std::min(biomap[metricname]->corr_alpha.size(), biomap[metricname]->corr_beta.size());
    
    for(int pc=0; pc<std::min(NUMBER_PCWEIGHTS, Npc); pc++){
        double pc_weight = biomap[metricname]->corr_alpha[pc] + biomap[metricname]->corr_beta[pc] * biometric_value;
        pcweights[pc] = pc_weight;
        pcslidervec[pc]->setValue((pc_weight+1)*50);
    }
    
    UpdateStatModelReader();        
    RenderModel();
}

vtkSmartPointer<vtkStringArray> MainWindow::GetFieldNameArray(int segmidx){
    std::vector<std::string> labels = {"cannon", "mc2", "mc4", "longpastern", "sesamoidleft", "sesamoidright", "shortpastern", "coffin", "navicular", "hoofcapsule"};
    
    vtkSmartPointer<vtkStringArray> namearray = vtkSmartPointer<vtkStringArray>::New();
    namearray->SetNumberOfValues(1);
    namearray->SetValue(0,labels[segmidx]);
    namearray->SetName("ModelName");
    
    return namearray;
}

vtkSmartPointer<vtkPolyData> MainWindow::AddLandmarkArray(vtkSmartPointer<vtkPolyData> poly){
    
    for(int i=0; i<poly->GetPointData()->GetNumberOfArrays(); i++){
        poly->GetPointData()->RemoveArray(i);
    }
    
    std::string modelname = vtkStringArray::SafeDownCast(poly->GetFieldData()->GetAbstractArray("ModelName"))->GetValue(0);
    
    std::vector<int> idlist;
    if(modelname=="hoofcapsule"){ idlist = {7414, 9592, 8, 9848, 1967}; } // model version 1: idlist = {6577,2085,6201,2284,11748}; }
    else if(modelname=="navicular"){ idlist = {175,176,146,147,71,69,0,1,3,5,9,12,13,25,17,22,93,126,153}; } // model version 1
    else{return poly;}
    
    vtkSmartPointer<vtkDoubleArray> array = vtkSmartPointer<vtkDoubleArray>::New();
    array->SetNumberOfComponents(3);
    array->SetNumberOfTuples(idlist.size());
    array->SetName("Landmarks");
    
    for(int i=0; i<idlist.size(); i++){
        double p[3]; poly->GetPoint(idlist[i], p);
        array->SetTuple(i, p);
    }
    
    poly->GetFieldData()->AddArray(array);
    
    return poly;
}



vtkSmartPointer<vtkMultiBlockDataSet> MainWindow::SplitIntoBlocks(bool linstructure){
    
    vtkSmartPointer<vtkMultiBlockDataSet> multiblock = vtkMultiBlockDataSet::SafeDownCast(skinning->GetOutput());

    vtkSmartPointer<vtkMultiBlockDataSet> output = vtkSmartPointer<vtkMultiBlockDataSet>::New();

    int labelcounter=0;
    
    
    for(int block_id=0; block_id<multiblock->GetNumberOfBlocks(); block_id++){
        
        vtkSmartPointer<vtkPolyData> poly = vtkPolyData::SafeDownCast(multiblock->GetBlock(block_id));    
        
    
        vtkSmartPointer<vtkPolyDataConnectivityFilter> connectivityFilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
        connectivityFilter->SetInputData(poly);
        connectivityFilter->SetExtractionModeToAllRegions();
        connectivityFilter->Update();
        
        int Nsegm = connectivityFilter->GetNumberOfExtractedRegions();
        
        
        if(Nsegm==1){  // only 1 segment
            poly->GetFieldData()->AddArray(GetFieldNameArray(labelcounter++));
            output->SetBlock(output->GetNumberOfBlocks(), poly);
        }
        else{  // In case of multiple segments, split them and order them in multiblock
            vtkSmartPointer<vtkMultiBlockDataSet> outputblock = vtkSmartPointer<vtkMultiBlockDataSet>::New();
            for(int segm=0; segm<Nsegm; segm++){
                vtkSmartPointer<vtkPolyDataConnectivityFilter> segmconnectivityFilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
                segmconnectivityFilter->SetInputData(poly);
                segmconnectivityFilter->SetExtractionModeToSpecifiedRegions();
                segmconnectivityFilter->AddSpecifiedRegion(segm); //select the region to extract here
                segmconnectivityFilter->Update();
                
                vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
                cleaner->SetInputData(segmconnectivityFilter->GetOutput());
                cleaner->Update();
                
                cleaner->GetOutput()->GetFieldData()->AddArray(GetFieldNameArray(labelcounter++));
                
                vtkSmartPointer<vtkPolyData> segmpoly = AddLandmarkArray(cleaner->GetOutput());
                
                if(linstructure){ output->SetBlock(output->GetNumberOfBlocks(), segmpoly); }
                else{ outputblock->SetBlock(outputblock->GetNumberOfBlocks(), segmpoly); }
            }
            if(!linstructure){ output->SetBlock(output->GetNumberOfBlocks(), outputblock); }
        }
    }
    
    return output;
}




