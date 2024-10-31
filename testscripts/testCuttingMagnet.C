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

        cutT = new TGeoTranslation(TranslationName.Data(),0,0,+dzshield - sensdz/2. - ilayer*(sensdz+passdz));
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

    auto *volMagTargetHCAL = new TGeoVolumeAssembly("volMagTargetHCAL");
    top->AddNode(volMagTargetHCAL,0,new TGeoTranslation(0,0,0));

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
    Double_t dz_6R = 242. *cm;

    Double_t xvertices_6R[nvertices] = {-0.,-0.,-33.,-33.,-0.,-0.,-77.,-77.};  
    Double_t yvertices_6R[nvertices] = {117.9,-117.9,-84.9,84.9,317.9,-317.9,-240.9,240.9};

    TGeoArb8 *arb = BuildArb8("arb",xvertices_6R,yvertices_6R,dz_6R);
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

    TGeoArb8 *arb_L = BuildArb8("arb_L",xvertices_6L,yvertices_6L,dz_6L);

    BooleanCut = PrepareBooleanOperation(volMagTargetHCAL, volSNDSensitiveLayer, nlayers_MagTargetHCAL, "SNDSensitiveLayer" , "arb_L", dz_6R, SensZ, FeZ);

    TGeoCompositeShape *csL = new TGeoCompositeShape("csL",BooleanCut.Data());
    auto Magn6_MiddleMagL = new TGeoVolume("Magn6_MiddleMagL",csL,med);
    Magn6_MiddleMagL->SetTransparency(1);

    top->AddNode(Magn6_MiddleMagR,0,new TGeoTranslation(0.,0.,0.));
    top->AddNode(Magn6_MiddleMagL,0,new TGeoTranslation(0.,0.,0.));

    const Double_t SiX = XDimension;
    const Double_t SiY = YDimension;
    const Double_t SiZ = 0.8;

    auto * SNDTargetSiliconLayer = new TGeoBBox("SNDTargetSiliconLayer", SiX/2., SiY/2., SiZ/2.);
    auto * volSNDTargetSiliconLayer = new TGeoVolume("volSNDTargetSiliconLayer",SNDTargetSiliconLayer,Silicon);
    volSNDTargetSiliconLayer->SetLineColor(kGreen);

    auto *volMagHCAL = new TGeoVolumeAssembly("volMagHCAL");

    const Int_t nlayers_MagHCAL = 34;
    Double_t dz_5L = 305. *cm;
    Double_t dz_5R = dz_5L;
    //building replicas of inner volumes of magnet 5 from FairShip
    RVec<Double_t> xvertices_5L = {0.,0.,22.,22.,0.,0.,32.,32.};  
    RVec<Double_t> yvertices_5L = {-230.9,230.9,-208.9,208.9,-66.9,66.9,34.9,-34.9};
    TGeoArb8 *arb_5L = BuildArb8("arb_5L",xvertices_5L.data(),yvertices_5L.data(),dz_5L);
    auto Magn5_MiddleMagL = new TGeoVolume("Magn5_MiddleMagL",arb_5L,med); //we first build it, then we use it

    RVec<Double_t> xvertices_5R = xvertices_5L * -1;
    RVec<Double_t> yvertices_5R = yvertices_5L * -1;

    cout<<xvertices_5R<<endl;
    cout<<yvertices_5R<<endl;
    TGeoArb8 *arb_5R = BuildArb8("arb_5R",xvertices_5R.data(),yvertices_5R.data(),dz_5R);
    auto Magn5_MiddleMagR = new TGeoVolume("Magn5_MiddleMagR",arb_5R,med); //we first build it, then we use it
    //cutting them and inserting sensitive volumes
    BooleanCut = PrepareBooleanOperation(volMagHCAL, volSNDTargetSiliconLayer, nlayers_MagHCAL,  "SNDSensitiveLayer" , "arb_5R", dz_5R, SiZ, FeZ);
    cout<<"Check of boolean operation magnet5"<<BooleanCut.Data()<<endl;
    TGeoCompositeShape *cs5R = new TGeoCompositeShape("cs5R",BooleanCut.Data());
    Magn5_MiddleMagR->SetShape(cs5R);
    Magn5_MiddleMagR->SetTransparency(1);

    BooleanCut = PrepareBooleanOperation(volMagHCAL, volSNDTargetSiliconLayer, nlayers_MagHCAL,  "SNDSensitiveLayer" , "arb_5L", dz_5L, SiZ, FeZ);
    TGeoCompositeShape *cs5L = new TGeoCompositeShape("cs5L",BooleanCut.Data());
    Magn5_MiddleMagL->SetShape(cs5L);
    Magn5_MiddleMagL->SetTransparency(1);

    //inserting the volumes
    top->AddNode(Magn5_MiddleMagR,0,new TGeoTranslation(0.,0.,-dz_6L-dz_5R-10.));
    top->AddNode(Magn5_MiddleMagL,0,new TGeoTranslation(0.,0.,-dz_6L-dz_5R-10.));
    top->AddNode(volMagHCAL,0,new TGeoTranslation(0,0,-dz_6L-dz_5R-10.));

    mygeometry->CloseGeometry();

    mygeometry->GetTopVolume()->Draw("ogl");

}