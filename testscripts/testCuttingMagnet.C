using namespace ROOT;
//It is time to CUT OUR LOSSSES (and our muon shield) - A. I. 30 October 2024
TGeoArb8* BuildArb8(const char* name, Double_t *xvertices, Double_t *yvertices, Double_t dz){
    auto *arb = new TGeoArb8(name,dz);
    arb->SetVertex(0,xvertices[0],yvertices[0]);
    arb->SetVertex(1,xvertices[1],yvertices[1]);
    arb->SetVertex(2,xvertices[2],yvertices[2]);
    arb->SetVertex(3,xvertices[3],yvertices[3]);
    arb->SetVertex(4,xvertices[4],yvertices[4]);
    arb->SetVertex(5,xvertices[5],yvertices[5]);
    arb->SetVertex(6,xvertices[6],yvertices[6]);
    arb->SetVertex(7,xvertices[7],yvertices[7]);
    return arb;
}

TString PrepareBooleanOperation(TGeoVolumeAssembly* TargetVolume, TGeoVolume * SensitiveVolume, const int nlayers, const char *sensname, const char* shieldname, Double_t dzshield, Double_t sensdz, Double_t passdz){
    //cut shield inserts and place sensitive layers inside
    TGeoTranslation *cutT;
    TString BooleanUnionShapes("(");
    //union of the sensitive layers positions to be inserted
    for (int ilayer = 0; ilayer < nlayers; ilayer++){
        TString TranslationName("cutT");
        TranslationName += TString(shieldname);

        cutT = new TGeoTranslation(TranslationName.Data(),0,0,+dzshield/2. - sensdz/2. - ilayer*(sensdz+passdz));
        cutT->RegisterYourself();
        BooleanUnionShapes += TString(sensname)+TString(":")+TranslationName;
        if (ilayer < nlayers-1) BooleanUnionShapes += " + ";

        TargetVolume->AddNode(SensitiveVolume,ilayer,cutT);
    }
    BooleanUnionShapes += ")";
    //cutting the shapes from the muon shield
    return (TString(shieldname)+TString(" - ")+BooleanUnionShapes);
}

void testCuttingMagnet(){
    const Double_t cm = 1.;
    const Double_t m = 100 * cm;

    TGeoManager * mygeometry = new TGeoManager("cutmagnet","Geometry with cut muon shield");
    //declare some "dummy" media for this test
    TGeoMaterial *mat = new TGeoMaterial("Vacuum",0,0,0);
    TGeoMedium   *med = new TGeoMedium("Vacuum",1,mat);
    TGeoMedium *Silicon = new TGeoMedium("Silicon",2,mat);
    TGeoMedium *Tungsten = new TGeoMedium("Tungsten",3,mat);

    TGeoBBox *mybox = new TGeoBBox("mybox",1000.*cm,1000.*cm,1000.*cm);
    TGeoVolume *top = new TGeoVolume("Top",mybox,med);
    gGeoManager->SetTopVolume(top);

    const Double_t XDimension = 40 * cm;
    const Double_t YDimension = 40 * cm;
    //layers to be cut
    const Double_t SensX = XDimension;
    const Double_t SensY = YDimension;
    const Double_t SensZ = 2.;

    const Double_t FeZ = 5.;

    Double_t dz_6R = 242. *cm *2; //better to define dZ as the full length, as usual

    auto *volMagTargetHCAL = new TGeoVolumeAssembly("volMagTargetHCAL");
    top->AddNode(volMagTargetHCAL,0,new TGeoTranslation(0,0,-dz_6R/2.));

    auto * SNDSensitiveLayer = new TGeoBBox("SNDSensitiveLayer", SensX/2., SensY/2., SensZ/2.);
    auto * volSNDSensitiveLayer = new TGeoVolume("volSNDSensitiveLayer",SNDSensitiveLayer, med); //TOP
    volSNDSensitiveLayer->SetLineColor(kCyan);

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

    const int nvertices = 8;
    //Magn6 MiddleMagR

    Double_t xvertices_6R[nvertices] = {-0.,-0.,-33.,-33.,-0.,-0.,-77.,-77.};  
    Double_t yvertices_6R[nvertices] = {117.9,-117.9,-84.9,84.9,317.9,-317.9,-240.9,240.9};

    TGeoArb8 *arb = BuildArb8("arb",xvertices_6R,yvertices_6R,dz_6R/2.);
    //start cutting away!!!!
    auto Magn6_MiddleMagR = new TGeoVolume("Magn6_MiddleMagR",arb,med); //we first build it, then we use it
    const Int_t nlayers_MagTargetHCAL = 50;

    TString BooleanCut = PrepareBooleanOperation(volMagTargetHCAL, volSNDSensitiveLayer, nlayers_MagTargetHCAL,  "SNDSensitiveLayer" , "arb", dz_6R, SensZ, FeZ);

    cout<<"Check of boolean operation "<<BooleanCut.Data()<<endl;

    TGeoCompositeShape *csR = new TGeoCompositeShape("csR",BooleanCut.Data());

    Magn6_MiddleMagR->SetShape(csR);
    Magn6_MiddleMagR->SetTransparency(1);

    Double_t dz_6L = dz_6R; 
    Double_t xvertices_6L[nvertices] = {-0.,-0.,33.,33.,-0.,-0.,77.,77.};  
    Double_t yvertices_6L[nvertices] = {117.9,-117.9,-84.9,84.9,317.9,-317.9,-240.9,240.9};

    TGeoArb8 *arb_L = BuildArb8("arb_L",xvertices_6L,yvertices_6L,dz_6L/2.);

    BooleanCut = PrepareBooleanOperation(volMagTargetHCAL, volSNDSensitiveLayer, nlayers_MagTargetHCAL, "SNDSensitiveLayer" , "arb_L", dz_6R, SensZ, FeZ);

    TGeoCompositeShape *csL = new TGeoCompositeShape("csL",BooleanCut.Data());
    auto Magn6_MiddleMagL = new TGeoVolume("Magn6_MiddleMagL",csL,med);
    Magn6_MiddleMagL->SetTransparency(1);


    Magn6_MiddleMagR->SetLineColor(kRed);
    Magn6_MiddleMagL->SetLineColor(kRed);
    top->AddNode(Magn6_MiddleMagR,0,new TGeoTranslation(0.,0.,-dz_6L/2.));
    top->AddNode(Magn6_MiddleMagL,0,new TGeoTranslation(0.,0.,-dz_6L/2.));

    const Double_t SiX = XDimension;
    const Double_t SiY = YDimension;
    const Double_t SiZ = 0.8;

    auto * SNDTargetSiliconLayer = new TGeoBBox("SNDTargetSiliconLayer", SiX/2., SiY/2., SiZ/2.);
    auto * volSNDTargetSiliconLayer = new TGeoVolume("volSNDTargetSiliconLayer",SNDTargetSiliconLayer,Silicon);
    volSNDTargetSiliconLayer->SetLineColor(kGreen);

    auto *volMagHCAL = new TGeoVolumeAssembly("volMagHCAL");

    const Int_t nlayers_MagHCAL = 34;
    Double_t dz_5L = 2*305. *cm;
    Double_t dz_5R = dz_5L;
    //building replicas of inner volumes of magnet 5 from FairShip
    RVec<Double_t> xvertices_5L = {0.,0.,22.,22.,0.,0.,32.,32.};  
    RVec<Double_t> yvertices_5L = {-230.9,230.9,-208.9,208.9,-66.9,66.9,34.9,-34.9};
    TGeoArb8 *arb_5L = BuildArb8("arb_5L",xvertices_5L.data(),yvertices_5L.data(),dz_5L/2.);
    auto Magn5_MiddleMagL = new TGeoVolume("Magn5_MiddleMagL",arb_5L,med); //we first build it, then we use it

    RVec<Double_t> xvertices_5R = xvertices_5L * -1;
    RVec<Double_t> yvertices_5R = yvertices_5L * -1;

    cout<<xvertices_5R<<endl;
    cout<<yvertices_5R<<endl;
    TGeoArb8 *arb_5R = BuildArb8("arb_5R",xvertices_5R.data(),yvertices_5R.data(),dz_5R/2.);
    auto Magn5_MiddleMagR = new TGeoVolume("Magn5_MiddleMagR",arb_5R,med); //we first build it, then we use it
    //cutting them and inserting sensitive volumes
    BooleanCut = PrepareBooleanOperation(volMagHCAL, volSNDTargetSiliconLayer, nlayers_MagHCAL,  "SNDSensitiveLayer" , "arb_5R", dz_5R, SiZ, FeZ);
    cout<<"Check of boolean operation magnet5"<<BooleanCut.Data()<<endl;
    TGeoCompositeShape *cs5R = new TGeoCompositeShape("cs5R",BooleanCut.Data());

    BooleanCut = PrepareBooleanOperation(volMagHCAL, volSNDTargetSiliconLayer, nlayers_MagHCAL,  "SNDSensitiveLayer" , "arb_5L", dz_5L, SiZ, FeZ);
    TGeoCompositeShape *cs5L = new TGeoCompositeShape("cs5L",BooleanCut.Data());

    //inserting the volumes
    Magn5_MiddleMagR->SetLineColor(kRed);
    Magn5_MiddleMagL->SetLineColor(kRed);
    top->AddNode(volMagHCAL,0,new TGeoTranslation(0,0,-dz_6L-dz_5R/2.-10.));

    //Emulsion Target

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

    

    //*****Silicon Target******//

    const Double_t EmTarget_SiTarget_Gap = 10.; //gap with Emulsion Target upstream

    const Double_t TungstenX = XDimension;
    const Double_t TungstenY = YDimension;
    const Double_t TungstenZ = 0.7;

    const Int_t nlayers_SiTarget = 58;

    const Double_t SiTargetX = XDimension;
    const Double_t SiTargetY = YDimension;
    const Double_t SiTargetZ = nlayers_SiTarget * (SiZ + TungstenZ);

    auto *SiTargetBox = new TGeoBBox("SiTargetBox",SiTargetX/2.,SiTargetY/2.,SiTargetZ/2.);
    auto *volSiTarget = new TGeoVolume("volSiTarget",SiTargetBox,med);

    //AddSensitiveVolume(volSNDTargetSiliconLayer) //uncomment when copying in actual class!

    auto * SNDTargetTungstenBlock = new TGeoBBox("SNDTargetTungstenBlock", TungstenX/2., TungstenY/2., TungstenZ/2.);
    auto * volSNDTargetTungstenBlock = new TGeoVolume("volSNDTargetTungstenBlock",SNDTargetTungstenBlock,Tungsten); //TOP
    volSNDTargetTungstenBlock->SetLineColor(kGray);

    for(Int_t n=0; n<nlayers_SiTarget; n++)
    {
      volSiTarget->AddNode(volSNDTargetTungstenBlock, n, new TGeoTranslation(0,0, -SiTargetZ/2. + n *(SiZ+TungstenZ) + TungstenZ/2. )); //W
      volSiTarget->AddNode(volSNDTargetSiliconLayer, n*1000, new TGeoTranslation(0,0,-SiTargetZ/2. + n *(SiZ+TungstenZ) + TungstenZ + SiZ/2 )); //Silicon
    }

    
    TGeoShapeAssembly * MagHCAL = static_cast<TGeoShapeAssembly*> (volMagHCAL->GetShape());
    MagHCAL->ComputeBBox(); //for an assembly needs to be computed
    Double_t dZ_MagHCAL = MagHCAL->GetDZ();
    cout<<"TEST "<<dZ_MagHCAL<<endl;
    //cutting the holes in the magnet for the big targets

    TGeoTranslation * T_SiTarget = new TGeoTranslation("T_SiTarget",0,0,+SiTargetZ/2.); 
    T_SiTarget->RegisterYourself();
    TGeoTranslation * T_EmTarget = new TGeoTranslation("T_EmTarget",0,0,-SiTargetZ/2.+EmTargetZDimension/2.); 
    T_EmTarget->RegisterYourself();

    TGeoCompositeShape *cs5L_si = new TGeoCompositeShape("cs5L_si","cs5L-SiTargetBox:T_SiTarget");
    TGeoCompositeShape *cs5L_siem = new TGeoCompositeShape("cs5L_siem","cs5L_si-EmTargetBox:T_EmTarget");
    Magn5_MiddleMagL->SetShape(cs5L_siem);
    Magn5_MiddleMagL->SetTransparency(1);

    TGeoCompositeShape *cs5R_si = new TGeoCompositeShape("cs5R_si","cs5R-SiTargetBox:T_SiTarget");
    TGeoCompositeShape *cs5R_siem = new TGeoCompositeShape("cs5R_siem","cs5R_si-EmTargetBox:T_EmTarget");
    Magn5_MiddleMagR->SetShape(cs5R_siem);
    Magn5_MiddleMagR->SetTransparency(1);

    top->AddNode(Magn5_MiddleMagR,0,new TGeoTranslation(0.,0.,-dz_6L-dz_5R/2.-10.));
    top->AddNode(Magn5_MiddleMagL,0,new TGeoTranslation(0.,0.,-dz_6L-dz_5R/2.-10.));

    top->AddNode(volSiTarget,0,new TGeoTranslation(0,0,-dz_6L-dz_5R/2.-10.+SiTargetZ/2.));
    top->AddNode(volEmTarget,0,new TGeoTranslation(0,0,-dz_6L-dz_5R/2.-10.-SiTargetZ/2.+EmTargetZDimension/2.));

    mygeometry->CloseGeometry();

    mygeometry->GetTopVolume()->Draw("ogl");

    //gGeoManager->CheckOverlaps();
    //gGeoManager->PrintOverlaps();

}