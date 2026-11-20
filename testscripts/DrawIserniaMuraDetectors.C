void DrawIserniaMuraDetectors(){

  float m = 100.;
  float cm =  1.;
  const float deg = 1.;

  // detector coordinates
  const double  xcA =   143.5*m;     // A location detector center abscissa
  const double xcB =   268.5*m;     // B location detector center abscissa

  const double ycA = 47.0*m;      // A location detector center ordinate
  const double ycB = 20*m;      // B location detector center ordinate

  const double zcA = 550.321*m;    // A location detector center applicate
  const double zcB = 573.578*m;    // B location detector center applicate

  
  const double EmulsionX = 19.0 *cm;
  const double EmulsionY = 19.0 *cm;
  const double PassiveSlabThickness = 0.1*cm;
  const double EmPlateWidth = 0.03 *cm;

  const double SensorHL = 60.0/2.*cm; // detector half length - X

  const double SensorHTh = 0.005/2.*cm; // detector half thickness - Y
  const double SensorHW = 20.0/2.*cm; // detector half width - Z

  TGeoManager * mygeometry = new TGeoManager("test","Test geometry");

  TGeoMaterial *mat = new TGeoMaterial("Vacuum",0,0,0);
  TGeoMedium   *med = new TGeoMedium("Vacuum",1,mat);
 
  TGeoBBox *worldSolid = new TGeoBBox("World",5000.*m,5000.*m,3000.*m);
  TGeoVolume *top = new TGeoVolume("Top",worldSolid,med);
  mygeometry->SetTopVolume(top);

  TGeoBBox* solidSensor1 = new TGeoBBox("Sensor1",SensorHL, SensorHW, SensorHTh);
  TGeoVolume* logicSensor1 = new TGeoVolume("logicSensor1",solidSensor1,med);
  TGeoRotation* rotm1 = new TGeoRotation("rotm1");
  rotm1->RotateX(80.0*deg);
  rotm1->RotateY(-28.52917855*deg);
  TGeoTranslation *trans = new TGeoTranslation("trans",xcA,ycA,zcA+SensorHW - 1.0*m);
  TGeoCombiTrans *combi = new TGeoCombiTrans(*trans,*rotm1);
  combi->RegisterYourself();
  top->AddNode(logicSensor1,1,combi); //Box ends at origin

  cout << "  *** physiSensor1  placed !"<< endl;
  

  // 2nd detector
  TGeoBBox* solidSensor2 = new TGeoBBox("Sensor2", SensorHL, SensorHW, SensorHTh);
  TGeoVolume* logicSensor2 = new TGeoVolume("Sensor2",solidSensor2, med);

  TGeoRotation* rotm2 = new TGeoRotation("rotm2");
  rotm2->RotateX(90.0*deg);
  TGeoTranslation *transB = new TGeoTranslation("transB",xcB,ycB,zcB+SensorHW - 1.0*m);
  TGeoCombiTrans *combiB = new TGeoCombiTrans(*transB,*rotm2);
  combiB->RegisterYourself();
  top->AddNode(logicSensor2,1,combiB);

  cout << "  *** physiSensor2  placed !"<< endl;
/*
//Envelope definition

   G4Box* Env1Solid = new G4Box("Env1",4000.*m,4000.*m,2000.*m);
   G4Box* solidSensor22 = new G4Box("Sensor22", 1.1*SensorHL, 1.1*SensorHTh, 1.1*SensorHW);
   G4SubtractionSolid* Env2Solid = new G4SubtractionSolid("Env2",Env1Solid,solidSensor22,0,G4ThreeVector(xcB,ycB,zcB + SensorHW - 1.0*m));
   G4RotationMatrix* rotm4 = new G4RotationMatrix();
   rotm4->rotateX(-10.0*deg);
   rotm4->rotateY(-28.52917855*deg);
   G4SubtractionSolid* EnvSolid = new G4SubtractionSolid("Envelope",Env2Solid,solidSensor22,rotm4,G4ThreeVector(xcA,ycA,zcA + SensorHW - 1.0*m));
   G4LogicalVolume* EnvLogical= new G4LogicalVolume(EnvSolid,nist->FindOrBuildMaterial("G4_Galactic"),"Envelope");
   G4VPhysicalVolume* EnvPhysical= new G4PVPlacement(0,G4ThreeVector(),EnvLogical,"Envelope",worldLogical,0,0,check);
  */
  mygeometry->CloseGeometry();

  mygeometry->GetTopVolume()->Draw("ogl");

  gGeoManager->CheckOverlaps(0.001);
} 