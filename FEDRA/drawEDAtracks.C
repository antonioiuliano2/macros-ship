void drawEDAtracks(){
 //which tracks do you want to draw?
 TCut mycut("nseg>20");

//getting tracks set from file
 EdbEDA* gEDA = new EdbEDA("linked_tracks.root",-1,mycut, kFALSE);
 
 EdbEDATrackSet *set = gEDA->GetTrackSet("TS"); 

// color selection
 set->SetColorMode(kCOLOR_BY_ID);
 gEDA->GetTrackSet("TS")->SetDraw(kTRUE);
 //gEDA->Redraw();
 gEDA->Run();


} 