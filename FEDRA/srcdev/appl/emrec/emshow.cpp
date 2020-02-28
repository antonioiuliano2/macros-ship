//-- Author :  Frank Meisel   11/November/2010

#include <string.h>
#include <iostream>
#include <TEnv.h>
#include "EdbLog.h"
#include "EdbScanProc.h"
#include "EdbShowerRec.h"

using namespace std;

void print_help_message()
{
    cout<< "\nUsage: \n\t  emshow -set=ID [ -o=DATA_DIRECTORY -v=DEBUG] \n\n";
    cout<< "\t\t  -set=ID   - id of the dataset formed as BRICK.PLATE.MAJOR.MINOR \n";
    cout<< "\t\t  DEBUG     - verbosity level: 0-print nothing, 1-errors only, 2-normal, 3-print all messages\n";
    cout<< "\t\t  -vt    - start showering from linked tracks attached to vertices with IP<250 (not yet included)\n";
    cout<< "\t\t  -lt    - start showering from linked tracks (all; standard)\n";
    cout<< "\t\t  -bt    - start showering from basetracks (all; useful if no linked tracks done yet)\n";
    cout<< "\nExample: \n";
    cout<< "\t  emshow -id=4554.10.1.0 -o/scratch/BRICKS -lt\n";
    cout<< "\n  The data location directory if not explicitly defined will be taken from .rootrc as: \n";
    cout<< "\t  emrec.outdir:      /scratch/BRICKS \n";
    cout<< "\t  emrec.EdbDebugLevel:      1\n";
    cout<<endl;
}

void set_default(TEnv &cenv)
{
    // default parameters for shower reconstruction
    // determined by experimental and simulation studies
    cenv.SetValue("emshow.cpcut","s.eW>13&&eCHI2P<2.5&&s1.eFlag>=0&&s2.eFlag>=0&&eN1==1&&eN2==1");
    cenv.SetValue("emshow.trkcut","nseg>1");
    cenv.SetValue("emshow.ConeRadius", 800);
    cenv.SetValue("emshow.ConeAngle", 0.02);
    cenv.SetValue("emshow.ConnectionDR", 150);
    cenv.SetValue("emshow.ConnectionDT", 0.15);
    cenv.SetValue("emshow.NPropagation", 3);
    cenv.SetValue("emshow.outdir", "..");
    cenv.SetValue("emshow.env", "shower.rootrc");
    cenv.SetValue("emshow.EdbDebugLevel", 1);
}

int main(int argc, char* argv[])
{
    if (argc < 2)   {
        print_help_message();
        return 0;
    }

    TEnv cenv("showerenv");
    set_default(cenv);
    gEDBDEBUGLEVEL        = cenv.GetValue("emshow.EdbDebugLevel", 1);
    const char *env       = cenv.GetValue("emshow.env", "shower.rootrc");
    const char *outdir    = cenv.GetValue("emshow.outdir", "./");

    bool      do_set      = false;
    bool      do_pred     = false;
    bool      do_VSB      = false;
    bool      do_showerfrom_vt      = false;
    bool      do_showerfrom_lt      = false;
    bool      do_showerfrom_bt      = false;

    Int_t     pred_plate  = 0, to_plate=0;
    Int_t     brick=0, plate=0, major=0, minor=0;

    for (int i=1; i<argc; i++ ) {
        char *key  = argv[i];

        if (!strncmp(key,"-set=",5))
        {
            if (strlen(key)>5)	sscanf(key+5,"%d.%d.%d.%d",&brick,&plate,&major,&minor);
            do_set=true;
        }
        else if (!strncmp(key,"-o=",3))
        {
            if (strlen(key)>3)	outdir=key+3;
        }
        else if (!strncmp(key,"-pred=",6))
        {
            if (strlen(key)>6)	{
                pred_plate = atoi(key+6);
                do_pred=true;
            }
        }
        else if (!strncmp(key,"-vt=",3))
        {
            if (strlen(key)>3)	{
                do_showerfrom_vt=true;
            }
        }
        if (!strncmp(key,"-lt=",3))
        {
            if (strlen(key)>3)	{
                do_showerfrom_lt=true;
            }
        }
        else if (!strncmp(key,"-bt=",3))
        {
            if (strlen(key)>3)	{
                do_showerfrom_bt=true;
            }
        }
        else if (!strncmp(key,"-VSB=",5))
        {
            if (strlen(key)>5)  {
                do_VSB=true;
                to_plate = atoi(key+5);
            }
        }
        else if (!strncmp(key,"-v=",3))
        {
            if (strlen(key)>3)	gEDBDEBUGLEVEL = atoi(key+3);
        }
    }

    if (!do_set)   {
        print_help_message();
        return 0;
    }

    cenv.SetValue("emshow.env", env);        
    
    cenv.ReadFile( cenv.GetValue("emshow.env", "shower.rootrc"),kEnvLocal);
    cenv.SetValue("emshow.outdir", outdir);

    cenv.WriteFile("shower.save.rootrc");

    if (do_set) {
        EdbScanProc sproc;
        sproc.eProcDirClient=outdir;
        printf("\n----------------------------------------------------------------------------\n");
        printf("tracking set %d.%d.%d.%d\n", brick,plate, major,minor);
        printf("----------------------------------------------------------------------------\n\n");

        EdbID id(brick,plate,major,minor);
        EdbScanSet *ss = sproc.ReadScanSet(id);
        ss->Brick().SetID(brick);
  
        EdbPVRec * eEdbPVRec = new EdbPVRec();
        TCut c = cenv.GetValue("emshow.cpcut","1");
        TCut trackcut = cenv.GetValue("emshow.trkcut","1");
 
        if (do_showerfrom_lt){
                //loop on plates
          int npl = ss->eIDS.GetEntries();

          for(int i=0; i<npl; i++) {
           EdbID *id = ss->GetID(i);
      
           EdbPlateP *plate = ss->GetPlate(id->ePlate);
           //read pattern information
           EdbPattern *p = new EdbPattern();
           sproc.ReadPatCPnopar(*p,*id, c);
           p->SetZ(plate->Z());
           p->SetSegmentsZ();
           p->SetID(i);
           p->SetPID(i);
           p->SetSegmentsPID();
          //plate->Print();
           p->Transform(    plate->GetAffineXY()   );
           p->TransformShr( plate->Shr() );
           p->TransformA(   plate->GetAffineTXTY() );
           p->SetSegmentsPlate(id->ePlate);
           eEdbPVRec->AddPattern(p);
          } //end of loop on patt
          //get tracks tree and apply selections
          TFile *inputfile = TFile::Open("linked_tracks.root");
          TTree *trackstree = (TTree*) inputfile->Get("tracks");

          TFile *outputfile = new TFile("selected_linkedtracks.root","RECREATE");
          TTree *selectedtrackstree = trackstree->CopyTree(trackcut);
          selectedtrackstree->Write();
          std::cout<<"Selected "<<selectedtrackstree->GetEntries()<<" tracks with cut "<<trackcut<<std::endl;

          float ConeRadius = cenv.GetValue("emshow.ConeRadius", 800);
          float ConeAngle = cenv.GetValue("emshow.ConeAngle", 0.02);
          float ConnectionDR = cenv.GetValue("emshow.ConnectionDR", 150);
          float ConnectionDT = cenv.GetValue("emshow.ConnectionDT", 0.15);
          float NPropagation = cenv.GetValue("emshow.NPropagation", 3);

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
          selectedtrackstree->SetBranchAddress("s",  &seg);

          const int ntracks = selectedtrackstree->GetEntries();


          for (int itrk = 0; itrk < ntracks; itrk++){
           selectedtrackstree->GetEntry(itrk);
    
           EdbSegP *segtest = new EdbSegP(*((EdbSegP*) seg->At(0)));
     
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
        else sproc.TrackSetBT(*ss,cenv);

    }
    if (do_VSB) {
        EdbScanProc sproc;
        sproc.eProcDirClient=outdir;
        printf("\n----------------------------------------------------------------------------\n");
        printf("prepare predictions %d.%d.%d.%d  for the plate %d\n", brick,plate, major,minor, pred_plate);
        printf("----------------------------------------------------------------------------\n\n");

        EdbID id(brick,plate,major,minor);
        EdbScanSet *ss = sproc.ReadScanSet(id);
        ss->Brick().SetID(brick);
        sproc.AssembleScanSet(*ss);

        EdbID predid(id);
        predid.ePlate = pred_plate;

        sproc.FindPredictionsRawSet(predid, *ss, to_plate);
        EdbPVRec ali;
        sproc.ReadFoundTracks(*ss, ali);
        sproc.WriteSBTracks(*(ali.eTracks), id);

    }

    return 1;
}
