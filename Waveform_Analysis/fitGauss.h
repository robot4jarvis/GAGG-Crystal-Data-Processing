#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <string>

double GEL(double *x, double *par) {  // gauss + erc + lin
    // Seven parameters to fit.
    double fval; double xval; xval = *x;
    fval = TMath::Gaus(xval, par[0], par[1],false) * par[2];

    return fval;
}

// This function, given a histigram a xmin and a xmax (ROI), fits a gaussian function
void fitGauss(TH1D *hist, double xmin, double xmax){

    // we draw it (just to see that it worked)
    TCanvas *c10 = new TCanvas();
    hist->GetXaxis()->SetRangeUser(xmin, xmax);
    double ymax = hist->GetBinContent(hist->GetMaximumBin());
    // We try fitting a gaussian with the following parameters
    TF1 *funGEL = new TF1("funGEL", GEL, xmin, xmax, 3); 
        // We defined a TF1 object (function) with name "funGEL"", function GausLin (above), 7 parameters
        funGEL->SetParameters((xmax+xmin)/2, (xmax-xmin)/2, ymax, 0,0,1,0); // We set initial parameters (just a guess)
        funGEL->SetParNames("Centroid","Sigma","PeakVal", "BackgroundSlope","ErcHeight","ErcMult","Baseline");   // We give the parameters a proper name (lineal part ax + b)
        funGEL->SetParLimits(0,xmin, xmax);
        funGEL->SetParLimits(1,0,(xmax-xmin));
        funGEL->SetParLimits(2,0,2*ymax);

    hist->Fit("gaus","R+"); //The R option restricts the fitting to [xmin, xmax]
    hist->Draw();

    double Res = (funGEL->GetParameter(1))*2.355; // We obtain the resolution
    double ERes = (funGEL->GetParError(1))*2.355;
    std::cout<<"Resolution (FWHM): "<< std::to_string(Res)  << " (" <<std::to_string(ERes)<<") ns\r\n";
}


