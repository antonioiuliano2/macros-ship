//simple memento of where I drawn the field map from

void drawfield(){
    TFile *inputfile = TFile::Open("$FAIRSHIP/field/GoliathFieldMap.root");
    TH3F* By = (TH3F*) inputfile->Get("By");
    By->SetTitle("By field;x[cm];y[cm];z[cm]");
    By->Draw("BOX2Z");
}