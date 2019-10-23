void drawEDAtracks(){
 //which tracks do you want to draw?
 const int ntracks = 1;
 int tracklist[ntracks] = {197};

//getting tracks set from file
 EdbEDA* gEDA = new EdbEDA("linked_tracks.root",-1,"1", kFALSE);
 
 EdbEDATrackSet *set = gEDA->GetTrackSet("TS"); 

//first run without drawing anything
 gEDA->GetTrackSet("TS")->SetDraw(kFALSE);

// reset tracks to be drawn
 set->ClearTracks(); 
 for(int trackID:tracklist){
  EdbTrackP *t = set->GetTrackBase(trackID);
  set->AddTrack(t);
 }
// color selection
 set->SetColorMode(kCOLOR_BY_ID);
 gEDA->GetTrackSet("TS")->SetDraw(kTRUE);
 //gEDA->Redraw();
 gEDA->Run();


} 