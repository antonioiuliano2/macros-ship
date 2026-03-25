vector<double> eff_formula(int found, int total){
  vector<double> efficiency; //value and error
  
  efficiency.push_back((double) found/total);
  double efferr = TMath::Sqrt(efficiency[0] * (1- efficiency[0])/total);
  efficiency.push_back(efferr);
  
  return efficiency;
  
}

void printresults(){
    int ntotal = 2500;
    ROOT::RVec<int> nselected = {2230,1911,1745,1187,584};

    for (auto& ngood: nselected){
        vector<double> efficiency = eff_formula(ngood, ntotal);
        cout<<"We have "<<efficiency[0]<<" pm "<<efficiency[1]<<endl;
    }

}