void drawEDAtracks(){
 //which tracks do you want to draw?
 const int ntracks = 2;
 int tracklist[ntracks] = {0,1};

//getting tracks set from file
 EdbEDA* gEDA = new EdbEDA("tracks_Giuliana.root",-1,"1", kFALSE);
 
 EdbEDATrackSet *set = gEDA->GetTrackSet("TS"); 

//first run without drawing anything
 gEDA->GetTrackSet("TS")->SetDraw(kFALSE);
 gEDA->Run();
// reset tracks to be drawn
 set->ClearTracks(); 
 for(int trackID:tracklist){
  EdbTrackP *t = set->GetTrackBase(trackID);
  set->AddTrack(t);
 }
// color selection
 set->SetColorMode(kCOLOR_BY_ID);
 gEDA->GetTrackSet("TS")->SetDraw(kTRUE);
 gEDA->Redraw();


} 