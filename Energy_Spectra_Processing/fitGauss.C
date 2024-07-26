double GEL(double *x, double *par) {  // gauss + erc + lin
    // Seven parameters to fit.
    double fval; double xval; xval = *x;
    fval = TMath::Gaus(xval, par[0], par[1],false) * par[2] + par[3]*(xval-par[0]) + par[4] * TMath::Erfc((xval-par[0])*par[5]) + par[6];

    return fval;

    // The function is:
    // fval = gauss(x, p0,p1) * p2 + p3*x + p4 * Erfx (x-p0)*p5 + p6
}

double fBACK(double *x, double *par) {  // Only the non-gaussian part of the previous function
    double fval; double xval; xval = *x;
    fval = par[3]*(xval-par[0]) + par[4] * TMath::Erfc((xval-par[0])*par[5]) + par[6];
    return fval;
}

// This function, given a histname (that must be inside a Histograms.root file), a xmin and a xmax (ROI), fits a gaussian function (+background)
void fitGauss(TString histName, double xmin, double xmax){

    // we open the file with the Histograms and retreive the second one
    TFile *f = new TFile("Histograms.root");
    TH1D *hist = (TH1D*)f->Get(histName);

    // we draw it (just to see that it worked)
    TCanvas *c1 = new TCanvas();
    hist->Draw();
    hist->GetXaxis()->SetRangeUser(xmin, xmax);
    double ymax = hist->GetBinContent(hist->GetMaximumBin());
    // We try fitting a gaussian with the following parameters
    TF1 *funGEL = new TF1("funGEL", GEL, xmin, xmax, 7); 
        // We defined a TF1 object (function) with name "funGEL"", function GausLin (above), 7 parameters
        funGEL->SetParameters((xmax+xmin)/2, (xmax-xmin)/2, ymax, 0,0,1,0); // We set initial parameters (just a guess)
        funGEL->SetParNames("Centroid","Sigma","PeakVal", "BackgroundSlope","ErcHeight","ErcMult","Baseline");   // We give the parameters a proper name (lineal part ax + b)
        funGEL->SetParLimits(0,xmin, xmax);
        funGEL->SetParLimits(1,0,(xmax-xmin));
        funGEL->SetParLimits(2,0,1.5*ymax);
        funGEL->SetParLimits(4,0,1.5*ymax);

    hist->Fit("funGEL","R+"); //The R parameter restricts the fitting to [xmin, xmax]

    
    TF1 *funBACK = new TF1("funBACK", fBACK,xmin,xmax,7);  // Only the Background
        funBACK->SetParameters(funGEL->GetParameters());
        funBACK->SetLineColor(kGreen);
    funBACK->Draw("SAME");

    double Res = (funGEL->GetParameter(1))*235.5 / (funGEL->GetParameter(0)); // We obtain the resolution
    std::cout<<"Relative resolution (FWHM/ADC): "<< std::to_string(Res) << "%\r\n";
}


