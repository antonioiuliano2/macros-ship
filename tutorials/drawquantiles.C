//froom ROOT reference
void drawquantiles(TH1D *h){

   const Int_t nq = 20;
   Double_t xq[nq];  // position where to compute the quantiles in [0,1]
   Double_t yq[nq];  // array to contain the quantiles
   for (Int_t i=0;i<nq;i++) xq[i] = Float_t(i+1)/nq;
   h->GetQuantiles(nq,yq,xq);
   //show the original histogram in the top pad
   TCanvas *c1 = new TCanvas("c1","demo quantiles",10,10,700,900);
   c1->Divide(1,2);
   c1->cd(1);
   h->Draw();
 
   // show the quantiles in the bottom pad
   c1->cd(2);
   gPad->SetGrid();
   TGraph *gr = new TGraph(nq,xq,yq);
   gr->SetMarkerStyle(21);
   gr->Draw("alp");

   gr->Print();

}
