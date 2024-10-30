void testrotation(double angle=0.){

  float cm = 1.;

  
  const double EmulsionX = 12.5 *cm;
  const double EmulsionY = 10.0 *cm;
  const double PassiveSlabThickness = 0.1*cm;
  const double EmPlateWidth = 0.03 *cm;

  TGeoManager * mygeometry = new TGeoManager("test","Test geometry");

  TGeoMaterial *mat = new TGeoMaterial("Vacuum",0,0,0);
  TGeoMedium   *med = new TGeoMedium("Vacuum",1,mat);
 
  TGeoBBox *mybox = new TGeoBBox("mybox",1000.*cm,1000.*cm,1000.*cm);
  TGeoVolume *top = new TGeoVolume("Top",mybox,med);
  mygeometry->SetTopVolume(top);

  TGeoBBox *EmulsionFilm = new TGeoBBox("EmulsionFilm", EmulsionX/2, EmulsionY/2, EmPlateWidth/2);
  TGeoVolume *volEmulsionFilm = new TGeoVolume("Emulsion",EmulsionFilm,med); //TOP
  volEmulsionFilm->SetLineColor(kBlue);

  TGeoVolumeAssembly *volTarget = new TGeoVolumeAssembly("volTarget");
  TGeoBBox *Passiveslab = new TGeoBBox("Passiveslab", EmulsionX/2, EmulsionY/2, PassiveSlabThickness/2);
  TGeoVolume *volPassiveslab = new TGeoVolume("volPassiveslab",Passiveslab,med);
  volPassiveslab->SetLineColor(kGray);
  double AllPlateWidth = EmPlateWidth + PassiveSlabThickness;
  //end of building geometry, close it and draw it
  Int_t nfilm = 1, npassiveslab = 1;
  const int fNPlates = 3;
  Double_t zpoint = 0.;
  TGeoRotation *rot = new TGeoRotation("rot");
  rot->RotateX(angle);
  TGeoTranslation *trans = new TGeoTranslation("trans",0,0,0);
  TGeoCombiTrans *combi = new TGeoCombiTrans(*trans,*rot);
  combi->RegisterYourself();
  top->AddNode(volTarget,1,combi); //Box ends at origin
  //adding nodes
	//adding emulsions
  for(Int_t n=0; n<fNPlates+1; n++)
	    {
	      volTarget->AddNode(volEmulsionFilm, nfilm, new TGeoTranslation(0,0,zpoint + EmPlateWidth/2 + n*AllPlateWidth));
	      nfilm++;
	    }           
	for(Int_t n=0; n<fNPlates; n++) //adding 1 mm lead plates
	    {
        volTarget->AddNode(volPassiveslab, npassiveslab, new TGeoTranslation(0,0,zpoint + EmPlateWidth + PassiveSlabThickness/2 + n*AllPlateWidth));
        npassiveslab++;
	    }	
  mygeometry->CloseGeometry();

  mygeometry->GetTopVolume()->Draw("ogl");

  gGeoManager->CheckOverlaps(0.001);
} 