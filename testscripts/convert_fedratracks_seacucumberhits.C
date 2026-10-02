struct MeshData {
    std::string name;
    std::string color;
    int transparency;
    std::vector<double> vertices;
    std::vector<int> indices;

    std::string toJson() const {
        std::ostringstream ss;
        ss << std::fixed << std::setprecision(1);
        ss << "{\"name\":\"" << name << "\","
           << "\"color\":\"" << color << "\","
           << "\"transparency\":" << transparency << ","
           << "\"vertices\":[";
        
        for (size_t i = 0; i < vertices.size(); ++i) {
            ss << vertices[i] << (i + 1 < vertices.size() ? "," : "");
        }
        
        ss << "],\"indices\":[";
        for (size_t i = 0; i < indices.size(); ++i) {
            ss << indices[i] << (i + 1 < indices.size() ? "," : "");
        }
        ss << "]},";
        return ss.str();
    }
};

MeshData createPlateMesh(int ePlate, double eXmin, double eXmax, double eYmin, double eYmax, double eZmin, double eZmax) {
    MeshData mesh;
    mesh.name = "Emulsion" + std::to_string(ePlate) + "#0";
    mesh.color = "#E8CFA0";
    mesh.transparency = 70;

    // 8 Vertices defining a 3D bounding box (or flat plate if eZmin == eZmax)
    // Bottom 4 vertices (Z = eZmin)
    // 0: (Xmin, Ymin, Zmin)
    // 1: (Xmin, Ymax, Zmin)
    // 2: (Xmax, Ymax, Zmin)
    // 3: (Xmax, Ymin, Zmin)
    // Top 4 vertices (Z = eZmax)
    // 4: (Xmin, Ymin, Zmax)
    // 5: (Xmin, Ymax, Zmax)
    // 6: (Xmax, Ymax, Zmax)
    // 7: (Xmax, Ymin, Zmax)
    
    mesh.vertices = {
        eXmin, eYmin, eZmin,  // 0
        eXmin, eYmax, eZmin,  // 1
        eXmax, eYmax, eZmin,  // 2
        eXmax, eYmin, eZmin,  // 3
        eXmin, eYmin, eZmax,  // 4
        eXmin, eYmax, eZmax,  // 5
        eXmax, eYmax, eZmax,  // 6
        eXmax, eYmin, eZmax   // 7
    };

    // 12 Triangles (36 indices) defining the 6 quad faces of the box
    mesh.indices = {
        0, 1, 5,   0, 5, 4,  // Side face - Xmin
        1, 2, 6,   1, 6, 5,  // Side face - Ymax
        2, 3, 7,   2, 7, 6,  // Side face - Xmax
        3, 0, 4,   3, 4, 7,  // Side face - Ymin
        0, 1, 2,   0, 2, 3,  // Bottom face - Zmin
        4, 5, 6,   4, 6, 7   // Top face - Zmax
    };

    return mesh;
}

void convert_fedratracks_seacucumberhits(){
    /*
    using Vec3 = std::array<double, 3>;
    //creating the hit container
    auto model = ROOT::RNTupleModel::Create();
    //auto mcP = model->MakeField<std::vector<SHiP::MCParticle>>("mcParticles");
    auto simH = model->MakeField<std::vector<SHiP::SimHit>>("simHits");
    //auto simR = model->MakeField<SHiP::SimResult>("simResult");
    auto writer = ROOT::RNTupleWriter::Recreate(std::move(model), "events", "demo_display_tracks.root");
    */
    //input tracks
    TFile *inputfile = TFile::Open("/home/utente/Simulations/RUN6_training_30_March_2021/b000006/b000006.0.0.0.trk.root","READ");
    TTree *tracktree = (TTree*) inputfile->Get("tracks");

    EdbDataProc *dproc = new EdbDataProc();
    TCut tracksel("nseg>10");
    dproc->InitVolume(100, tracksel);
    EdbPVRec *ali = new EdbPVRec();
    ali = dproc->PVR();
    ali->FillCell(30,30,0.009,0.009);
    //tracktree.SetAlias("trk","t.") #points create confusion to python
    const std::string outDir = "web/data";
    //setting branches

    EdbEDAAreaSet *areaset = new EdbEDAAreaSet();
    areaset->SetAreas(ali);
    areaset->Print(); 

    for(int itrk = 0; itrk < ali->eTracks->GetEntries(); itrk++){

        const std::string ep = outDir + "/event_" + std::to_string(itrk) + ".json";
        std::ofstream ej(ep);
        if (!ej) {
            std::cerr << "warning: cannot write '" << ep << "'\n";
            continue;
        }
        ej << "{\"event\":" << itrk << ",\"hits\":[";
        //simH->clear(); //for now, each event is a track
        EdbTrackP * tr = (EdbTrackP*) ali->eTracks->At(itrk);
       // temptrack->Copy(EdbTrackP(trk));
        //start loop on segments associated to the track
        int nseg = tr->N();
        cout<<"Track: " << itrk << " with " << nseg << " segments" <<endl;
        for (int i = 0; i< nseg; i++){
            EdbSegP *seg = tr->GetSegmentF(i);
            /*
            SHiP::SimHit h;
            h.detectorId = seg->Plate();
            h.trackId = seg->Track();
            h.pdgCode = seg->Flag(); //if montecarlo, provide pdgcode
            Vec3 vecpos = {seg->X(), seg->Y(), seg->Z()}; 
            Vec3 vecmom = {seg->TX(), seg->TY(), 1};  //if known momentum, can also provide actual momentum magnitude
            h.position = vecpos;
            h.momentum = vecmom;
            h.energyDeposit = seg->W();
            h.time = 1.;
            h.pathLength = 1.;
            simH->push_back(h);
            */

            if (i) ej << ',';
            ej << "{\"x\":" << seg->X() *1e-3 << ",\"y\":" << seg->Y() *1e-3
               << ",\"z\":" << seg->Z() *1e-3 << ",\"e\":" << seg->W()
               << ",\"pdg\":" << seg->Flag() << "}";
        }
        ej << "]";
        ej << "}\n";
        //simR->hits = *simH;  // bundle
        //writer->Fill();
    }
    for (int iplate=0; iplate<=areaset->N(); iplate++){
        EdbEDAArea *platearea = areaset->GetArea(iplate);
        if(platearea){
            MeshData plateMesh = createPlateMesh(platearea->Plate(), platearea->Xmin()*1e-3, platearea->Xmax()*1e-3, platearea->Ymin()*1e-3, platearea->Ymax()*1e-3, platearea->Z()*1e-3, platearea->Z()*1e-3+0.190);
           std::cout << plateMesh.toJson() << std::endl;
        }
    }
}