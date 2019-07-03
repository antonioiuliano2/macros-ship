/*simple script to loop over events, MC tracks and hit points
Works without any include, provided that we are in FairShip environment,
just launch root -l rpc_loop() in the folder with the simulation output*/

void rpc_loop(TString *inputfile);

void rpc_loop(TString *inputfile){

 TFile *file = TFile::Open(inputfile->Data()); 
 if (!file) return;
 TTreeReader reader("cbmsim",file);

 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
 TTreeReaderArray<ShipRpcPoint> rpcpoints(reader,"ShipRpcPoint");

 TDatabasePDG *mypdg = TDatabasePDG::Instance(); //to know information about particles

 //**************variables*************************
 double vx,vy,vz; //neutrino interaction vertex position

 double x,y,z; //hit positions
 double px, py, pz, p, pt;//momenta
 int charge;
 int pdgcode, detectorID, trackID;
 //*****************DECLARING HISTOGRAMS TO BE FILLED********************
 TH1D *hweight = new TH1D("hweight","Weight of events",250,0,2500);
 TProfile2D *hweightmap = new TProfile2D("hweightmap", "Weight map in the yz plane",400,-3400,-3000,440,-220,220);

 TH2D * hyz =new TH2D("hyz","yz distribution of rpcpoints",400,-2800,-2400,440,-220,220);
 TH2D * hxz =new  TH2D("hxz","xz distribution of rpcpoints",400,-2800,-2400,440,-220,220);

 TH2D * hppt =new TH2D("hppt","momentum vs transverse momentum",200,0,200,100,0,10);

 int nevents = reader.GetEntries(false);
 cout<<"Number of events: "<<nevents<<endl;
 int ientry = 0;
 double totalweight = 0.;
 double weight;
 //*****************START MAIN EVENT LOOP********************************************
 while (reader.Next()){
     ientry = reader.GetCurrentEntry();// keeps track of the number of event (from 0 to Nevents - 1)    
     //getting weight information
     if (ientry%100==0) cout<<"arrived at event "<<ientry<<endl;

     ShipMCTrack firsttrack = tracks[0]; //weight is the same for all the tracks of the same event
     vx = firsttrack.GetStartX();
     vy = firsttrack.GetStartY();
     vz = firsttrack.GetStartZ();

     weight = firsttrack.GetWeight(); //weight defined as sum (x_i pho_i)
     totalweight += weight; //it can be useful to know total weight for normalization (i.e. how many events have surpassed my cuts)

     hweight->Fill(weight);
     hweightmap->Fill(vz,vy,weight);
     //access the hits       
     for (const ShipRpcPoint& rpcpoint: rpcpoints){       

         x = rpcpoint.GetX();
         y = rpcpoint.GetY();
         z = rpcpoint.GetZ();

         px = rpcpoint.GetPx();
         py = rpcpoint.GetPy();
         pz = rpcpoint.GetPz();
         p = pow(pow(px,2)+pow(py,2)+pow(pz,2),0.5);
         pt = pow(pow(px,2)+pow(py,2),0.5); //beam is along z

         detectorID = rpcpoint.GetDetectorID();
         trackID = rpcpoint.GetTrackID();         
         pdgcode = rpcpoint.PdgCode();
         
         charge = 0;
         if(mypdg->GetParticle(pdgcode) != NULL) charge = mypdg->GetParticle(pdgcode)->Charge(); //unknown particles lead to NULL pointer
         if (TMath::Abs(charge) > 0){ //we limit ourselves to charged particles
          //filling weighted histograms
 
          hyz->Fill(z,y,weight);
          hxz->Fill(z,x,weight);

          hppt->Fill(p,pt,weight);
          }
         } //end loop on rpc points
         
     }//end of loop on events
 //setting statistics box options and position (by default it covers the color bar)
 gStyle->SetStatX(0.3);
 gStyle->SetStatY(0.9);  
 gStyle->SetOptStat(11111);

 TCanvas *cweight = new TCanvas();
 cweight->Divide(1,2);
 cweight->cd(1);
 hweight->Draw();
 cweight->cd(2);
 hweightmap->Draw("COLZ");

 cout<<"We have simulated "<<nevents<<" events, after weighting "<< totalweight<<endl;

 TCanvas *cz = new TCanvas();
 cz->Divide(1,2);
 cz->cd(2);
 hyz->GetXaxis()->SetTitle("z[cm]");
 hyz->GetYaxis()->SetTitle("y[cm]");
 hyz->Draw("COLZ");
 cz->cd(1);
 hxz->GetXaxis()->SetTitle("z[cm]");
 hxz->GetYaxis()->SetTitle("x[cm]");
 hxz->Draw("COLZ"); 

 TCanvas *cppt = new TCanvas();
 hppt->GetXaxis()->SetTitle("p[GeV]");
 hppt->GetYaxis()->SetTitle("pt[GeV]");
 hppt->Draw("COLZ");

 }

void rpc_loop(){
 TString *inputfilename = new TString("/eos/user/a/aiuliano/public/sims_FairShip/sim_nutau/neutrinos2019/8/nu_mu/ship.conical.Genie-TGeant4.root");
 rpc_loop(inputfilename);
}
