void testSNDMS(){

  const Double_t cm = 1.;
  const Double_t m = 100 * cm;

  
  
  TGeoManager * mygeometry = new TGeoManager("test","Test geometry");
  //declare some "dummy" media for this test
  TGeoMaterial *mat = new TGeoMaterial("Vacuum",0,0,0);
  TGeoMedium   *med = new TGeoMedium("Vacuum",1,mat);

  TGeoMedium *Silicon = new TGeoMedium("Silicon",2,mat);
  TGeoMedium *Tungsten = new TGeoMedium("W",3,mat);
  TGeoMedium *Iron = new TGeoMedium("Fe",4,mat);
 
  TGeoBBox *mybox = new TGeoBBox("mybox",1000.*cm,1000.*cm,1000.*cm);
  TGeoVolume *top = new TGeoVolume("Top",mybox,med);
  gGeoManager->SetTopVolume(top);
  const Double_t CenterZ = 0.;
  //parameters

  const Int_t nbricks = 5;
  cout<<"Test using "<<nbricks<< " brick "<<endl;
   
  const Int_t n_plates = 36;

  //emulsion film parameters
  const Double_t EmTh = 0.0070 * cm;
  const Double_t EmX = 40. * cm;
  const Double_t EmY = 40. * cm;
  const Double_t PbTh = 0.0175 * cm;
  const Double_t PassiveTh = 0.1 * cm;
  const Double_t EPlW = 2* EmTh + PbTh;
  const Double_t AllPW = PassiveTh + EPlW;
  //brick packaging borders
  const Double_t BrPackX = 2*0.0 *cm; // 1 mm gap (for now removed)
  const Double_t BrPackY = 2*0.0 *cm;
  const Double_t BrPackZ = 0. * cm; //no z border at the moment
  
  const Double_t BrickX = EmX + BrPackX;
  const Double_t BrickY = EmY + BrPackY;

  const Double_t BrickZ = n_plates * AllPW + EPlW + BrPackZ;
  
  const Double_t TTrackerZ = 2. * cm;
 
  //const Double_t ZDimension = 22 *cm;
  const Double_t XDimension = EmX;
  const Double_t YDimension = EmY;
  const Double_t EmTargetZDimension = nbricks* BrickZ + nbricks *TTrackerZ;
  
  //declaring volumes
  TGeoBBox *EmTargetBox = new TGeoBBox("EmTargetBox",XDimension/2, YDimension/2, EmTargetZDimension/2);
  TGeoVolume *volEmTarget = new TGeoVolume("volEmTarget",EmTargetBox, med);

  TGeoBBox *Brick = new TGeoBBox("brick", BrickX/2, BrickY/2, BrickZ/2);
  TGeoVolume *volBrick = new TGeoVolume("Brick",Brick,med);
  volBrick->SetLineColor(kCyan);
  //volBrick->SetTransparency(1);   

  //adding ECC bricks
  Double_t d_cl_z = - EmTargetZDimension/2;
  //adding also TT heres
  TGeoBBox *TT = new TGeoBBox("TT", XDimension/2, YDimension/2, (TTrackerZ)/2);
  TGeoVolume *volTT = new TGeoVolume("TargetTracker",TT,med); //TOP

  volTT->SetLineColor(kBlue);

  for(int l = 0; l < nbricks; l++)
	{
	  volEmTarget->AddNode(volBrick,l,new TGeoTranslation(0, 0, d_cl_z +BrickZ/2));
    volEmTarget->AddNode(volTT,l,new TGeoTranslation(0, 0, d_cl_z +BrickZ+TTrackerZ/2));
        
	  //6 cm is the distance between 2 columns of consecutive Target for TT placement
	  d_cl_z += BrickZ + TTrackerZ;
	}	

  
  //emulsion films
  TGeoBBox *EmulsionFilm = new TGeoBBox("EmulsionFilm", EmX/2, EmY/2, EmTh/2);
  TGeoVolume *volEmulsionFilm = new TGeoVolume("Emulsion",EmulsionFilm,med); //TOP
  TGeoVolume *volEmulsionFilm2 = new TGeoVolume("Emulsion2",EmulsionFilm,med); //BOTTOM
  volEmulsionFilm->SetLineColor(kBlue);
  volEmulsionFilm2->SetLineColor(kBlue);
  //plastic base

  TGeoBBox *PlBase = new TGeoBBox("PlBase", EmX/2, EmY/2, PbTh/2);
  TGeoVolume *volPlBase = new TGeoVolume("PlasticBase",PlBase,med);
  volPlBase->SetLineColor(kYellow-4);

  //passive tungsten layers
  TGeoBBox *Passive = new TGeoBBox("Passive", EmX/2, EmY/2, PassiveTh/2);
  TGeoVolume *volPassive = new TGeoVolume("volPassive",Passive,med);
  volPassive->SetLineColor(kGray);

    
  for(Int_t n=0; n<n_plates; n++)
    {
      volBrick->AddNode(volPassive, n, new TGeoTranslation(0,0,-BrickZ/2+BrPackZ/2+ EPlW + PassiveTh/2 + n*AllPW)); //LEAD
    }

  for(Int_t n=0; n<n_plates+1; n++)
    {
      volBrick->AddNode(volEmulsionFilm2, n, new TGeoTranslation(0,0,-BrickZ/2+BrPackZ/2+ EmTh/2 + n*AllPW)); //BOTTOM
      volBrick->AddNode(volEmulsionFilm, n, new TGeoTranslation(0,0,-BrickZ/2+BrPackZ/2+3*EmTh/2+PbTh+n*AllPW)); //TOP
      volBrick->AddNode(volPlBase, n, new TGeoTranslation(0,0,-BrickZ/2+BrPackZ/2+EmTh+PbTh/2+n*AllPW)); //PLASTIC BASE
    }  

  auto *volTarget = new TGeoVolumeAssembly("volTarget");

  Double_t dZ_subtargets = 0.; //avoid repeating the previous translations;
  
  
  top->AddNode(volTarget,1,new TGeoTranslation(0,0,0.));

  dZ_subtargets += CenterZ+EmTargetZDimension/2.;
  volTarget->AddNode(volEmTarget, 1, new TGeoTranslation(0.,0.,dZ_subtargets));
  dZ_subtargets += EmTargetZDimension/2.;
  
  //*****Silicon Target******//

  const Double_t EmTarget_SiTarget_Gap = 10.; //gap with Emulsion Target upstream

  const Double_t SiX = XDimension;
  const Double_t SiY = YDimension;
  const Double_t SiZ = 0.8;

  const Double_t TungstenX = XDimension;
  const Double_t TungstenY = YDimension;
  const Double_t TungstenZ = 0.7;

  const Int_t nlayers_SiTarget = 58;

  const Double_t SiTargetX = XDimension;
  const Double_t SiTargetY = YDimension;
  const Double_t SiTargetZ = nlayers_SiTarget * (SiZ + TungstenZ);

  auto *SiTargetBox = new TGeoBBox("SiTargetBox",SiTargetX/2.,SiTargetY/2.,SiTargetZ/2.);
  auto *volSiTarget = new TGeoVolume("volSiTarget",SiTargetBox,med);

  auto * SNDTargetSiliconLayer = new TGeoBBox("SNDTargetSiliconLayer", SiX/2., SiY/2., SiZ/2.);
  auto * volSNDTargetSiliconLayer = new TGeoVolume("volSNDTargetSiliconLayer",SNDTargetSiliconLayer,Silicon); //TOP
  volSNDTargetSiliconLayer->SetLineColor(kGreen);
  //AddSensitiveVolume(volSNDTargetSiliconLayer) //uncomment when copying in actual class!

  auto * SNDTargetTungstenBlock = new TGeoBBox("SNDTargetTungstenBlock", TungstenX/2., TungstenY/2., TungstenZ/2.);
  auto * volSNDTargetTungstenBlock = new TGeoVolume("volSNDTargetTungstenBlock",SNDTargetTungstenBlock,Tungsten); //TOP
  volSNDTargetTungstenBlock->SetLineColor(kGray);

  for(Int_t n=0; n<nlayers_SiTarget; n++)
    {
      volSiTarget->AddNode(volSNDTargetTungstenBlock, n, new TGeoTranslation(0,0, -SiTargetZ/2. + n *(SiZ+TungstenZ) + TungstenZ/2. )); //W
      volSiTarget->AddNode(volSNDTargetSiliconLayer, n, new TGeoTranslation(0,0,-SiTargetZ/2. + n *(SiZ+TungstenZ) + TungstenZ + SiZ/2 )); //Silicon
    }

  dZ_subtargets += EmTarget_SiTarget_Gap +SiTargetZ/2.;
  volTarget->AddNode(volSiTarget,1, new TGeoTranslation(0,0,dZ_subtargets));
  dZ_subtargets += SiTargetZ/2.;

  //*****Magnetic HCAL******//
  
  const Double_t SiTarget_MagHCAL_Gap = 10.; //gap with Emulsion Target upstream

  const Double_t FeX = XDimension;
  const Double_t FeY = YDimension;
  const Double_t FeZ = 5.;

  const Int_t nlayers_MagHCAL = 34;

  const Double_t MagHCALX = XDimension;
  const Double_t MagHCALY = YDimension;
  const Double_t MagHCALZ = nlayers_MagHCAL * (SiZ + FeZ);

  auto *MagHCALBox = new TGeoBBox("MagHCALBox", MagHCALX/2., MagHCALY/2., MagHCALZ/2.);
  auto *volMagHCAL = new TGeoVolume("volMagHCAL",MagHCALBox,med);

  auto * SNDTargetIronBlock = new TGeoBBox("SNDTargetIronBlock", FeX/2., FeY/2., FeZ/2.);
  auto * volSNDTargetIronBlock = new TGeoVolume("volSNDTargetIronBlock",SNDTargetIronBlock,Iron); //TOP
  volSNDTargetIronBlock->SetLineColor(kRed);

  for(Int_t n=0; n<nlayers_MagHCAL; n++)
    {
      volMagHCAL->AddNode(volSNDTargetIronBlock, n, new TGeoTranslation(0,0, -MagHCALZ/2. + n *(SiZ + FeZ) + FeZ/2. )); //W
      volMagHCAL->AddNode(volSNDTargetSiliconLayer, n+nlayers_SiTarget, new TGeoTranslation(0,0,-MagHCALZ/2. + n *(SiZ + FeZ) + FeZ + SiZ/2. )); //Silicon
    }

  dZ_subtargets += SiTarget_MagHCAL_Gap + MagHCALZ/2.;
  volTarget->AddNode(volMagHCAL,1, new TGeoTranslation(0,0,dZ_subtargets));
  dZ_subtargets += MagHCALZ/2.;

  //*****Magnetic Target + HCAL******//

  const Double_t MagHCAL_MagTargetHCAL_Gap = 10.; //gap with Emulsion Target upstream


  const Int_t nlayers_MagTargetHCAL = 50;

  const Double_t SensX = XDimension;
  const Double_t SensY = YDimension;
  const Double_t SensZ = 2.;

  const Double_t MagTargetHCALX = XDimension;
  const Double_t MagTargetHCALY = YDimension;
  const Double_t MagTargetHCALZ = nlayers_MagTargetHCAL * (SensZ + FeZ);

  auto *MagTargetHCALBox = new TGeoBBox("MagTargetHCALBox", MagTargetHCALX/2., MagTargetHCALY/2., MagTargetHCALZ/2.);
  auto *volMagTargetHCAL = new TGeoVolume("volMagTargetHCAL",MagTargetHCALBox,med);

  auto * SNDSensitiveLayer = new TGeoBBox("SNDSensitiveLayer", SensX/2., SensY/2., SensZ/2.);
  auto * volSNDSensitiveLayer = new TGeoVolume("volSNDSensitiveLayer",SNDSensitiveLayer, med); //TOP
  volSNDSensitiveLayer->SetLineColor(kCyan);
  //AddSensitiveVolume(volSNDSensitiveLayer);//to be uncommented in the class geometry code

  for(Int_t n=0; n<nlayers_MagTargetHCAL; n++)
    {
      volMagTargetHCAL->AddNode(volSNDTargetIronBlock, n+nlayers_MagHCAL, new TGeoTranslation(0,0, -MagTargetHCALZ/2. + n *(SensZ + FeZ) + FeZ/2. )); //W
      volMagTargetHCAL->AddNode(volSNDSensitiveLayer, n, new TGeoTranslation(0,0,-MagTargetHCALZ/2. + n *(SensZ + FeZ) + FeZ + SensZ/2. )); //Silicon
    }

  dZ_subtargets +=  MagHCAL_MagTargetHCAL_Gap + MagTargetHCALZ/2.;
  volTarget->AddNode(volMagTargetHCAL,1, new TGeoTranslation(0,0,dZ_subtargets));
  dZ_subtargets += MagTargetHCALZ/2.;

  //composition of sensitive volumes
  const Double_t SciFiX = XDimension;
  const Double_t SciFiY = YDimension;
  const Double_t SciFiZ = 0.5;

  const Double_t ScintX = XDimension;
  const Double_t ScintY = YDimension;
  const Double_t ScintZ = 1.5;

  auto * SNDSciFi = new TGeoBBox("SNDSciFi", SciFiX/2., SciFiY/2., SciFiZ/2.);
  auto * volSNDSciFi = new TGeoVolume("volSNDSciFi",SNDSciFi, med); //TOP
  volSNDSciFi->SetLineColor(kGreen);

  auto * SNDScint = new TGeoBBox("SNDScint", ScintX/2., ScintY/2., ScintZ/2.);
  auto * volSNDScint = new TGeoVolume("volSNDScint",SNDScint, med); //TOP
  volSNDScint->SetLineColor(kBlue);

  volSNDSensitiveLayer->AddNode(volSNDSciFi,0, new TGeoTranslation(0.,0.,-SensZ/2. + SciFiZ/2.));
  volSNDSensitiveLayer->AddNode(volSNDScint,0, new TGeoTranslation(0.,0.,-SensZ/2. + SciFiZ + ScintZ/2.));

  //end of building geometry, close it and draw it

  mygeometry->CloseGeometry();

  volTarget->InspectShape();

  mygeometry->GetTopVolume()->Draw("ogl");
  //mygeometry->GetVolume("volTarget")->Draw("ogl");
  gGeoManager->SetVisLevel(4);

  gGeoManager->CheckOverlaps(0.001);
  gGeoManager->PrintOverlaps();
} 
