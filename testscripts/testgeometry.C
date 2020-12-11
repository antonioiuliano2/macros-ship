void testgeometry(double stereoangle){

  float cm = 1.;

  
  
  TGeoManager * mygeometry = new TGeoManager("test","Test geometry");

  TGeoMaterial *mat = new TGeoMaterial("Vacuum",0,0,0);
  TGeoMedium   *med = new TGeoMedium("Vacuum",1,mat);
 
  TGeoBBox *mybox = new TGeoBBox("mybox",1000.*cm,1000.*cm,1000.*cm);
  TGeoVolume *top = new TGeoVolume("Top",mybox,med);
  gGeoManager->SetTopVolume(top);

  const Double_t mufilterXDIM = 62.0 *cm;
  const Double_t mufilterYDIM = 60.5 *cm;
  const Double_t mufilterZDIM = 22 *cm;

  const Double_t ironXDIM = mufilterXDIM;
  const Double_t ironYDIM = mufilterYDIM;
  const Double_t ironZDIM = 20 *cm;

  const int nbars = 77;
  const Double_t overlapbars = 0. *cm;

  const Double_t stereoXDIM = mufilterXDIM/TMath::Cos(TMath::DegToRad()*stereoangle);
  const Double_t stereoYDIM = mufilterYDIM+mufilterXDIM*TMath::Tan(TMath::DegToRad()*stereoangle);

  const Double_t stereoZDIM = 2 * cm;

  const Double_t stereodiagonalDIM = TMath::Sqrt(stereoXDIM * stereoXDIM + stereoYDIM * stereoYDIM);
  
  const Double_t barXDIM = stereoXDIM; //larger than anywhere around X;
  const Double_t barYDIM = (stereoYDIM + overlapbars/TMath::Cos(TMath::DegToRad()*stereoangle) * (nbars - 1))/nbars * TMath::Cos(TMath::DegToRad()*stereoangle);
  const Double_t barZDIM = 1*cm;

  cout<<barYDIM<<endl;
  //motherbox
  TGeoBBox *mufilterbox = new TGeoBBox("mufilterbox",mufilterXDIM/2., mufilterYDIM/2., mufilterZDIM/2.);
  TGeoVolume *volmufilter = new TGeoVolume("volmufilter",mufilterbox,med);
  volmufilter->SetLineColor(kGray);

  top->AddNode(volmufilter, 1, new TGeoTranslation(0,0,0));

  //first simple box, then test stereo box
  TGeoBBox *ironbox = new TGeoBBox("ironbox",ironXDIM/2., ironYDIM/2., ironZDIM/2.);
  TGeoVolume *voliron = new TGeoVolume("voliron",ironbox,med);
  voliron->SetLineColor(kGray);

  TGeoVolumeAssembly *StereoPlane = new TGeoVolumeAssembly("StereoPlane");

  volmufilter->AddNode(voliron, 1, new TGeoTranslation(0,0,-mufilterZDIM/2. + ironZDIM/2.));

  TGeoTranslation * stereotrans = new TGeoTranslation(0,0,-mufilterZDIM/2. + ironZDIM + stereoZDIM/2.);
  TGeoRotation *stereorot = new TGeoRotation("stereorot",stereoangle,0,0);
  stereorot->SetName("stereorot");
  stereorot->RegisterYourself();
  TGeoCombiTrans * stereo = new TGeoCombiTrans(*stereotrans, *stereorot);
  volmufilter->AddNode(StereoPlane,1,stereotrans);
  
  //adding bars

  TGeoBBox *barbox =  new TGeoBBox("barbox",barXDIM/2.,barYDIM/2.,barZDIM/2.);
  
  TGeoVolume *volbar[nbars];
  TGeoCompositeShape *cutbar[nbars];
 
  for (int ibar = 0; ibar < nbars; ibar++){    
    
    Double_t dy_bar = -stereoYDIM/2. + barYDIM/2. + (barYDIM-overlapbars)/TMath::Cos(TMath::DegToRad()*stereoangle)*ibar; 
    Double_t dz_bar_hor = -stereoZDIM/2. + barZDIM/2. * (2 *(ibar%2) + 1.); //on the left or right side of the volume

    TGeoTranslation *yztrans = new TGeoTranslation(barYDIM*TMath::Cos(TMath::DegToRad()*stereoangle)/2.,dy_bar,dz_bar_hor);
    yztrans->SetName(Form("yztrans[%i]",ibar));
    yztrans->RegisterYourself();
    TGeoCombiTrans *yzstereo = new TGeoCombiTrans(*yztrans, *stereorot);
    yzstereo->SetName(Form("yzstereo[%i]",ibar));
    yzstereo->RegisterYourself();

    cutbar[ibar] = new TGeoCompositeShape("cutbar",Form("(barbox:yzstereo[%i])*mufilterbox",ibar));
    volbar[ibar] = new TGeoVolume(Form("volbar[%i]",ibar),cutbar[ibar],med);
    volbar[ibar]->SetLineColor(kRed);
    
    StereoPlane->AddNode(volbar[ibar],ibar+1E+3,new TGeoTranslation(0,0,0));
  }
  
  //end of building geometry, close it and draw it

  mygeometry->CloseGeometry();

  mygeometry->GetVolume("volmufilter")->Draw("ogl");

  gGeoManager->CheckOverlaps(0.001);
} 

void testallangles(){
  const double stereomin = 10.;
  const double stereomax = 80.;
  const double stereostep = 5.;
  for (int istep = stereomin; istep<=stereomax;istep++){
    testgeometry(istep);
    TGLViewer * v = (TGLViewer *)gPad->GetViewer3D();
    v->SetCurrentCamera(TGLViewer::kCameraOrthoXOY );
    v->SetStyle(TGLRnrCtx::kWireFrame);
    if (istep == stereomin) v->SavePicture("plots/stereoviews.gif");
    else v->SavePicture("plots/stereoviews.gif+");
  }

}
