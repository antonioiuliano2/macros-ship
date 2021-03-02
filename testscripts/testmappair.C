//test code to count how many segment have from a MCEvent MCTrack pair
//to launch it please do (in a new nusrv9 terminal)
//source ~/newroot_setup.sh
//cd /ship/DESY2019/macros/testgiulianapair
//root -l
//>> .L testmappair.C+
//>> testmappair()
//ROOT libraries
#include "TGraph2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TClonesArray.h"
//FEDRA libraries
#include "/ship/sw_test/fedra_newroot/include/EdbSegP.h"
#include "/ship/sw_test/fedra_newroot/include/EdbPattern.h"
//C++ std libraries
#include <map>
#include <utility>
#include <iostream>
using namespace std;

void testmappair(){
    //reading input tracks file and setting branch
    TFile * inputfile = TFile::Open("linked_tracks.root");
    TTree *inputtree = (TTree*) inputfile->Get("tracks");

    TClonesArray *s = new TClonesArray("EdbSegP",10000);
    int nseg = 0;
    inputtree->SetBranchAddress("s",&s);
    inputtree->SetBranchAddress("nseg",&nseg);

    //declaring variables to be used in the loop;
    map<pair<int,int>, int> frequencyEvent;
    EdbTrackP *trk = NULL;
    EdbSegP *seg = NULL;
    int trk_mcevt, trk_mctrk;

    const int ntracks = inputtree->GetEntries();
    //loop over volume tracks
    for (int itrk = 0; itrk < ntracks; itrk++){
     inputtree->GetEntry(itrk);
     //loop over track segments
     for(Long64_t iseg=0;iseg<nseg; iseg++){
                seg=(EdbSegP*)s->At(iseg);
                trk_mcevt=seg->MCEvt();
                trk_mctrk=seg->MCTrack();
                if (frequencyEvent.find(pair<int,int>(trk_mcevt, trk_mctrk)) == frequencyEvent.end()) //the first time does not exist yet, I create it with value 0
                 frequencyEvent[pair<int,int>(trk_mcevt,trk_mctrk)] = 0;
                frequencyEvent[pair<int,int>(trk_mcevt, trk_mctrk)]++;
            }
    }
    //check if map works, plotting a graph
    TGraph2D *grids = new TGraph2D();
    map<pair<int,int>, int>::iterator it;
    int ipoint = 0;
    cout<<"Starting map loop "<<endl;
    for (it = frequencyEvent.begin(); it != frequencyEvent.end(); it++){
        //the key is a pair (first), the value is the counts (second)
        pair<int,int> ids = it->first;
        int counts = it->second;
        //the key itself, being a pair, has a first (MCEvent) and a second (MCTrack)
        grids->SetPoint(ipoint, ids.first,ids.second,counts);
        ipoint++;
    }
    TCanvas *ctest = new TCanvas();
    grids->SetTitle("How many segments for this track?;MCEvent;MCTrack");
    grids->Draw("COLZ");
}
