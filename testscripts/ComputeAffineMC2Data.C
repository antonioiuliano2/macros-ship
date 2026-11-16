//make the histogram of data angles
void ComputeAffineMC2Data(){
    const int npoints_datavsMC = 1;
    float txdata[npoints_datavsMC] = {0.72};
    float tydata[npoints_datavsMC] = {0.18}; 
    float txMC[npoints_datavsMC] = {0.26};
    float tyMC[npoints_datavsMC] = {0.21};

    EdbAffine2D *MC2Data = new EdbAffine2D();
    //MC2Data->Calculate(npoints_datavsMC,txMC,tyMC,txdata,tydata,2); //2 means calculate full
    MC2Data->CalculateShift(npoints_datavsMC,txMC,tyMC,txdata,tydata); 
    MC2Data->Print();
}
