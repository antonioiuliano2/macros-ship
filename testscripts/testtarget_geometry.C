void testtarget_geometry(){

  const Double_t cm = 1.;
  const Double_t m = 100 * cm;

  
  
  TGeoManager * mygeometry = new TGeoManager("test","Test geometry");

  TGeoMaterial *mat = new TGeoMaterial("Vacuum",0,0,0);
  TGeoMedium   *med = new TGeoMedium("Vacuum",1,mat);
 
  TGeoBBox *mybox = new TGeoBBox("mybox",1000.*cm,1000.*cm,1000.*cm);
  TGeoVolume *top = new TGeoVolume("Top",mybox,med);
  gGeoManager->SetTopVolume(top);
  const Double_t CenterZ = 0.;
  //parameters
  const Int_t nrows = 4; //double with respect to SNDLHC
  const Int_t ncols = 2; //as with SNDLHC
  const Int_t nwalls = 13; //start from five, then see how changes

  const Int_t nbricks = nrows * ncols * nwalls;
  cout<<"Test using "<<nbricks<< " brick "<<endl;
   
  const Int_t n_plates = 59;
  
  const Double_t Ydist = 0.0 * cm; //already within BrPackY
  //emulsion film parameters
  const Double_t EmTh = 0.0070 * cm;
  const Double_t EmX = 19.2 * cm;
  const Double_t EmY = 19.2 * cm;
  const Double_t PbTh = 0.0175 * cm;
  const Double_t PassiveTh = 0.1 * cm;
  const Double_t EPlW = 2* EmTh + PbTh;
  const Double_t AllPW = PassiveTh + EPlW;
  //brick packaging borders
  const Double_t BrPackX = 2*0.05 *cm; // 1 mm gap
  const Double_t BrPackY = 2*0.05 *cm;
  const Double_t BrPackZ = 0. * cm; //no z border at the moment
  
  const Double_t BrickX = EmX + BrPackX;
  const Double_t BrickY = EmY + BrPackY;

  const Double_t BrickZ = n_plates * AllPW + EPlW + BrPackZ;
 

  //wall sizes
  const Double_t WallXDim = ncols*BrickX;
  const Double_t WallYDim = nrows*BrickY+(nrows-1)*Ydist;
  const Double_t WallZDim = 10. * cm; // including borders
  
  const Double_t TTrackerZ = 3. * cm;
 
  const Double_t XDimension = WallXDim; //for now, same as wall
  const Double_t YDimension = WallYDim; 
  //const Double_t ZDimension = 22 *cm;
  const Double_t ZDimension = nwalls* WallZDim + nwalls *TTrackerZ;
  
  //mufilter 
  const Double_t FeBlockX = 1.950*m;
  const Double_t FeBlockY = 3.850*m;
  const Double_t FeThickZ = 30.* cm; //everything thin, to start with  
  const Double_t FeThinZ = 10.* cm; //everything thin, to start with
  
  const Int_t nblocksthick = 2;
  const Int_t nblocksthin = 4;


  const Double_t XRpc = FeBlockX;
  const Double_t YRpc = FeBlockY;
  const Double_t ZRpc = 8.*cm;

  const Int_t nrpc = nblocksthick + nblocksthin;

  const Double_t MuFilterX = XRpc;
  const Double_t MuFilterY = YRpc;
  const Double_t MuFilterZ = nrpc * ZRpc + nblocksthick * FeThickZ + nblocksthin * FeThinZ;

  const Double_t EmuTargetGapMuFilter = 25 * cm;
  //declaring volumes
  TGeoBBox *TargetBox = new TGeoBBox("TargetBox",XDimension/2, YDimension/2, ZDimension/2);
  TGeoVolume *volTarget = new TGeoVolume("volTarget",TargetBox, med);
  
  TGeoBBox *Wall = new TGeoBBox("wall",XDimension/2, YDimension/2, WallZDim/2);
  TGeoVolume *volWall = new TGeoVolume("Wall",Wall,med);

  volWall->SetLineColor(kGray);
  
  TGeoBBox *Row = new TGeoBBox("row",XDimension/2, BrickY/2, BrickZ/2);
  TGeoVolume *volRow = new TGeoVolume("Row",Row,med);

  volRow->SetLineColor(kGray);

  TGeoBBox *Brick = new TGeoBBox("brick", BrickX/2, BrickY/2, BrickZ/2);
  TGeoVolume *volBrick = new TGeoVolume("Brick",Brick,med);
  volBrick->SetLineColor(kCyan);
  //volBrick->SetTransparency(1);   
    
  //adding volumes , bricks, rows of bricks, walls 
  Double_t d_cl_x = -WallXDim/2;
  for(int j= 0; j < ncols; j++)
	{
	  volRow->AddNode(volBrick,j,new TGeoTranslation(d_cl_x+BrickX/2, 0, 0));
	  d_cl_x += BrickX;
	}

  Double_t d_cl_y = -WallYDim/2;
  for(int k= 0; k< nrows; k++)
	{
	  volWall->AddNode(volRow,k,new TGeoTranslation(0, d_cl_y + BrickY/2, 0));
        
	  // 2mm is the distance for the structure that holds the brick
	  d_cl_y += BrickY + Ydist;
	}

  //adding walls
  Double_t d_cl_z = - ZDimension/2;
  //adding also TT heres
  TGeoBBox *TT = new TGeoBBox("TT", XDimension/2, YDimension/2, (TTrackerZ)/2);
  TGeoVolume *volTT = new TGeoVolume("TargetTracker",TT,med); //TOP

  volTT->SetLineColor(kYellow);

  for(int l = 0; l < nwalls; l++)
	{
	  volTarget->AddNode(volWall,l,new TGeoTranslation(0, 0, d_cl_z +WallZDim/2));
    volTarget->AddNode(volTT,l,new TGeoTranslation(0, 0, d_cl_z +WallZDim+TTrackerZ/2));
        
	  //6 cm is the distance between 2 columns of consecutive Target for TT placement
	  d_cl_z += WallZDim + TTrackerZ;
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
  volPassive->SetTransparency(1);
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

  top->AddNode(volTarget,1,new TGeoTranslation(0,0,CenterZ+ZDimension/2.));
  //end of building geometry, close it and draw it

  //moving towards mufilter

  TGeoBBox *FeBlockThick = new TGeoBBox("FeBlockThick", FeBlockX/2, FeBlockY/2, FeThickZ/2);
  TGeoVolume *volFeBlockThick = new TGeoVolume("FeBlockThick",FeBlockThick,med); //THICK
  volFeBlockThick->SetLineColor(kGreen);

  TGeoBBox *FeBlockThin = new TGeoBBox("FeBlockThin", FeBlockX/2, FeBlockY/2, FeThinZ/2);
  TGeoVolume *volFeBlockThin = new TGeoVolume("FeBlockThin",FeBlockThin,med); //THIN
  volFeBlockThin->SetLineColor(kGreen);

  TGeoBBox *RPC = new TGeoBBox("RPC", XRpc/2.,YRpc/2.,ZRpc/2.);
  TGeoVolume *volRPC = new TGeoVolume("volRPC",RPC,med); //RPC
  volRPC->SetLineColor(kRed);

  TGeoBBox *MuFilter = new TGeoBBox("MuFilter",MuFilterX/2., MuFilterY/2., MuFilterZ/2.);
  TGeoVolume *volMuFilter = new TGeoVolume("volMuFilter",MuFilter,med); //MuFilter

  volMuFilter->SetLineColor(kGray);
  
  //top->AddNode()

  //adding volumes

  Double_t d_mu_z = -MuFilterZ/2.;
  for (int n=0;n< nblocksthick;n++){
    volMuFilter->AddNode(volFeBlockThick,n,new TGeoTranslation(0,0,d_mu_z + FeThickZ/2.));
    volMuFilter->AddNode(volRPC,n,new TGeoTranslation(0,0,d_mu_z + FeThickZ + ZRpc/2.));

    d_mu_z = d_mu_z + FeThickZ + ZRpc; 
  }
  
  for (int n=0;n< nblocksthin;n++){
    volMuFilter->AddNode(volFeBlockThin,n,new TGeoTranslation(0,0,d_mu_z + FeThinZ/2.));
    volMuFilter->AddNode(volRPC,n + nblocksthick,new TGeoTranslation(0,0,d_mu_z + FeThinZ + ZRpc/2.));

    d_mu_z = d_mu_z + FeThinZ + ZRpc; 
  }

  top->AddNode(volMuFilter,1,new TGeoTranslation(0,0,CenterZ + ZDimension + EmuTargetGapMuFilter + MuFilterZ/2.));



  mygeometry->CloseGeometry();

  volTarget->InspectShape();
  volMuFilter->InspectShape();

  mygeometry->GetTopVolume()->Draw("ogl");
  //mygeometry->GetVolume("volTarget")->Draw("ogl");
  gGeoManager->SetVisLevel(2);

  gGeoManager->CheckOverlaps(0.001);
  gGeoManager->PrintOverlaps();
} 
