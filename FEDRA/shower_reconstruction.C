/*void print_help_message()
{
    cout<< "\nUsage: \n\t  emshowerrec -set=ID [ -o=DATA_DIRECTORY -v=DEBUG] \n\n";
    cout<< "\t\t  -set=ID   - id of the dataset formed as BRICK.PLATE.MAJOR.MINOR \n";
    cout<< "\t\t  DEBUG     - verbosity level: 0-print nothing, 1-errors only, 2-normal, 3-print all messages\n";
    cout<< "\t\t  -vtr    - start showering from linked tracks attached to vertices with IP<250 (not yet included)\n";
    cout<< "\t\t  -ltr    - start showering from linked tracks (all; standard)\n";
    cout<< "\t\t  -btr    - start showering from basetracks (all; useful if no linked tracks done yet)\n";
    cout<< "\t\t  -otr    - only tracking, no showering\n";
    cout<< "\nExample: \n";
    cout<< "\t  emshow -id=4554.10.1.0 -o/scratch/BRICKS -lt\n";
    cout<< "\n  The data location directory if not explicitly defined will be taken from .rootrc as: \n";
    cout<< "\t  emrec.outdir:      /scratch/BRICKS \n";
    cout<< "\t  emrec.EdbDebugLevel:      1\n";
    cout<<endl;
}*/

void set_default(TEnv &cenv)
{
    // default parameters for shower reconstruction
    // determined by experimental and simulation studies
    cenv.SetValue("showerrec.cpcut","s.eW>13&&eCHI2P<2.5&&s1.eFlag>=0&&s2.eFlag>=0&&eN1==1&&eN2==1");
    cenv.SetValue("showerrec.trkcut","nseg>1");
    cenv.SetValue("showerrec.firstsegment",0);
    cenv.SetValue("showerrec.ConeRadius", 800);
    cenv.SetValue("showerrec.ConeAngle", 0.1);
    cenv.SetValue("showerrec.ConnectionDR", 150);
    cenv.SetValue("showerrec.ConnectionDT", 0.15);
    cenv.SetValue("showerrec.NPropagation", 3);
    cenv.SetValue("showerrec.outdir", "..");
    cenv.SetValue("showerrec.env", "showerrec.rootrc");
    cenv.SetValue("showerrec.EdbDebugLevel", 1);
}

void shower_reconstruction()
{
/*    if (argc < 2)   {
        print_help_message();
        return 0;
}*/
    int brick = 5;
    int nplate = 0;
    int major = 0;
    int minor = 0;
    
    TEnv cenv("showerrecenv");
    set_default(cenv);
    gEDBDEBUGLEVEL        = cenv.GetValue("showerrec.EdbDebugLevel", 1);
    const char *env       = cenv.GetValue("showerrec.env", "showerrec.rootrc");
    const char *outdir    = cenv.GetValue("showerrec.outdir", "..");

    cenv.SetValue("showerrec.env", env);        
    
    cenv.ReadFile( cenv.GetValue("showerrec.env", "showerrec.rootrc"),kEnvLocal);
    cenv.SetValue("showerrec.outdir", outdir);

    cenv.WriteFile("showerrec.save.rootrc");

    
    EdbScanProc sproc;
    sproc.eProcDirClient=outdir;
    printf("\n----------------------------------------------------------------------------\n");
    printf("tracking set %d.%d.%d.%d\n", brick,nplate, major,minor);
    printf("----------------------------------------------------------------------------\n\n");

    EdbID id(brick,nplate,major,minor);
    EdbScanSet *ss = sproc.ReadScanSet(id);
    ss->Brick().SetID(brick);
  
    EdbPVRec * eEdbPVRec = new EdbPVRec();
    TCut c = cenv.GetValue("showerrec.cpcut","1");
    TCut trackcut = cenv.GetValue("showerrec.trkcut","1");
    bool do_showerfrom_lt = true;
    if (do_showerfrom_lt){
                //loop on plates
          int npl = ss->eIDS.GetEntries();

          for(int i=0; i<npl; i++) {
           EdbID *idplate = ss->GetID(i);
      
           EdbPlateP *plate = ss->GetPlate(idplate->ePlate);
           //read pattern information
           EdbPattern *p = new EdbPattern();
           sproc.ReadPatCPnopar(*p,*idplate, c);
           p->SetZ(plate->Z());
           p->SetSegmentsZ();
           p->SetID(i);
           p->SetPID(i);
           p->SetSegmentsPID();
          //plate->Print();
           p->Transform(    plate->GetAffineXY()   );
           p->TransformShr( plate->Shr() );
           p->TransformA(   plate->GetAffineTXTY() );
           p->SetSegmentsPlate(idplate->ePlate);
           eEdbPVRec->AddPattern(p);
          } //end of loop on patt
          //get tracks tree and apply selections
          TFile *inputfile = TFile::Open("linked_tracks.root");
          TTree *trackstree = (TTree*) inputfile->Get("tracks");

          TFile *outputfile = new TFile("selected_linkedtracks.root","RECREATE");
          TTree *selectedtrackstree = trackstree->CopyTree(trackcut);
          selectedtrackstree->Write();
          std::cout<<"Selected "<<selectedtrackstree->GetEntries()<<" tracks with cut "<<trackcut<<std::endl;

          float ConeRadius = cenv.GetValue("showerrec.ConeRadius", 800);
          float ConeAngle = cenv.GetValue("showerrec.ConeAngle", 0.1);
          float ConnectionDR = cenv.GetValue("showerrec.ConnectionDR", 150);
          float ConnectionDT = cenv.GetValue("showerrec.ConnectionDT", 0.15);
          float NPropagation = cenv.GetValue("showerrec.NPropagation", 3);

          EdbShowerRec *eShowerRec = new EdbShowerRec();
          //Setting parameters
          eShowerRec->SetAlgoParameter(ConeRadius,0);
          eShowerRec->SetAlgoParameter(ConeAngle,1);
          eShowerRec->SetAlgoParameter(ConnectionDR,2);
          eShowerRec->SetAlgoParameter(ConnectionDT,3);
          eShowerRec->SetAlgoParameter(NPropagation,4);


               // Print parameters
          eShowerRec->PrintParameters();
    
         // Create Initiator BT array:
          TObjArray * eInBTArray=new TObjArray();

    // Reset eShowerRec Arrays: InBTArray and RecoShowerArray....
          eShowerRec->ResetInBTArray();
          eShowerRec->ResetRecoShowerArray();

          TClonesArray *seg  = new TClonesArray("EdbSegP", 60);
          int nseg;       
          
          selectedtrackstree->SetBranchAddress("s",  &seg);
          selectedtrackstree->SetBranchAddress("nseg",  &nseg);       

          const int ntracks = selectedtrackstree->GetEntries();
          int firstsegment = cenv.GetValue("showerrec.firstsegment",0); //take first segment of track
          
          int whichsegment;
          for (int itrk = 0; itrk < ntracks; itrk++){
           selectedtrackstree->GetEntry(itrk);
           if (firstsegment == 0) whichsegment = nseg - 1;
           else whichsegment = 0;       
           EdbSegP *segtest = new EdbSegP(*((EdbSegP*) seg->At(whichsegment)));           
           //set array with inBT (array of initiator base tracks)
           eInBTArray->Add(segtest);                  
          }
          eShowerRec->SetInBTArray(eInBTArray);
          eShowerRec->PrintInitiatorBTs();

          //set edbpvrec
          eShowerRec->SetEdbPVRec(eEdbPVRec);

          cout << " eShowerRec->SetUseAliSub(0)..." << endl;
          eShowerRec->SetUseAliSub(0);

          cout << " eShowerRec->Execute()..." << endl;

          //Start actual reconstruction
          eShowerRec->Execute();

          //Print output
          eShowerRec->PrintRecoShowerArray();
          outputfile->Close();
         }   
}
