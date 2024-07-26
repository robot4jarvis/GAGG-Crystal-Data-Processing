#include "calibrate.h"

double GEL(double *x, double *par) {  // gauss + erc + lin
    // Seven parameters to fit.
    double fval; double xval; xval = *x;
    fval = TMath::Gaus(xval, par[0], par[1],false) * par[2] + par[3]*(xval-par[0]) + par[4] * TMath::Erfc((xval-par[0])*par[5]) + par[6];

    return fval;

    // The function is:
    // fval = gauss(x, p0,p1) * p2 + p3*x + p4 * Erfx (x-p0)*p5 + p6
}

double GEL2(double *x, double *par) { // double gaussian + background. The parameters of the extra gaussian are 7,8 and 9 (10 in total)
    double fval; double xval; xval = *x;
    fval = TMath::Gaus(xval, par[0], par[1],false)*par[2]  + TMath::Gaus(xval, par[7],par[8],false) * par[9] + par[3]*(xval-par[0]) + par[4] * TMath::Erfc((xval-par[0])*par[5]) + par[6];

    return fval;
}

double pol2(double *x, double *par) {  // function y = a*x*x + b*x + c
    double fval; double xval; xval = *x;
    fval = xval *xval* par[0] + xval* par[1] + par[2];
    return fval;
}

// This function, given a histname (that must be inside a Histograms.root file), a xmin and a xmax (ROI), fits a gaussian function (+background) and returns the fitted function
TF1 *fitGEL(TH1D *hist, double xmin, double xmax){

    // we open the file with the Histograms and retreive the second one
    //TFile *f = new TFile("Histograms.root");
    //TH1D *hist = (TH1D*)f->Get(histName);

    // we draw it (just to see that it worked)
    hist->GetXaxis()->SetRange(xmin, xmax);
    double ymax = hist->GetBinContent(hist->GetMaximumBin());
    // We try fitting a gaussian with the following parameters
    TF1 *funGEL = new TF1("funGEL", GEL, xmin, xmax, 7); 
        // We defined a TF1 object (function) with name "funGEL"", function GausLin (above), 7 parameters
        funGEL->SetParameters((xmax+xmin)/2, (xmax-xmin)/2, ymax/2, 0,0,1,0); // We set initial parameters (just a guess)
        funGEL->SetParNames("Centroid","Sigma","PeakVal", "BackgroundSlope","ErcHeight","ErcMult","Baseline");   // We give the parameters a proper name (lineal part ax + b)
        funGEL->SetParLimits(0,xmin, xmax);
        funGEL->SetParLimits(1,0,(xmax-xmin));
        funGEL->SetParLimits(2,0,1.5*ymax);
        funGEL->SetParLimits(4,0,1.5*ymax);

    hist->Fit("funGEL","QR+"); //The R parameter restricts the fitting to [xmin, xmax]

    return funGEL;
}

TF1 *fitGEL2(TH1D *hist, double xmin, double xmid, double xmax){

    // we open the file with the Histograms and retreive the second one
    //TFile *f = new TFile("Histograms.root");
    //TH1D *hist = (TH1D*)f->Get(histName);

    // we draw it (just to see that it worked)
    hist->GetXaxis()->SetRangeUser(xmin, xmax);
    double ymax = hist->GetBinContent(hist->GetMaximumBin());
    // We try fitting a gaussian with the following parameters
    TF1 *funGEL2 = new TF1("funGEL2", GEL2, xmin, xmax, 10); 
        // We defined a TF1 object (function) with name "funGEL"", function GausLin (above), 7 parameters
        funGEL2->SetParameters((xmid+xmin)/2, (xmid-xmin)/2, ymax, 0,0,1,0,(xmax+xmid)/2, (xmax-xmid)/2,ymax); // We set initial parameters (just a guess)
        funGEL2->SetParNames("Centroid1","Sigma1","PeakVal1", "BackgroundSlope","ErcHeight","ErcMult","Baseline","Centroid2","Sigma2","PeakVal2");   // We give the parameters a proper name (lineal part ax + b)
        funGEL2->SetParLimits(0,xmin, xmid);            funGEL2->SetParLimits(7,xmid, xmax);
        funGEL2->SetParLimits(1,0,(xmax-xmin));         funGEL2->SetParLimits(8,0,(xmax-xmin));
        funGEL2->SetParLimits(2,0,1.5*ymax);            funGEL2->SetParLimits(9,0,1.5*ymax);
        funGEL2->SetParLimits(4,0,1.5*ymax);

    hist->Fit("funGEL2","QR+"); //The R parameter restricts the fitting to [xmin, xmax]

    return funGEL2;
}

void macro6(){
    // This macro opens a root file, reads their histograms and fits a gaussian where indicated.
    // Then, using calibrate.h, it calibrates all histograms in the indicated file.
    // Rigth now, it is set up to use the runs 10-13
    
    // We open the Histograms file and extract all data (that we will use for calibration)
    TFile *f = new TFile("Histograms.root");
    TH1D *histNa22 = (TH1D*)f->Get("Na22");
    TH1D *histCs137 = (TH1D*)f->Get("Cs137");
    TH1D *histCo60 = (TH1D*)f->Get("Co60");
    TH1D *histBACK = (TH1D*)f->Get("BACK");

    // We select all the regions of interest
    TCanvas *cNa22 = new TCanvas(); histNa22->Draw();   TF1 *fun1Na22 = fitGEL(histNa22, 1800,2500);
                                                        TF1 *fun2Na22 = fitGEL(histNa22, 4800,5800);
    TCanvas *cCs137 = new TCanvas(); histCs137->Draw(); TF1 *funCs137 = fitGEL(histCs137, 2400,3200);
    TCanvas *cCo60 = new TCanvas(); histCo60->Draw();   TF1 *funCo60 = fitGEL2(histCo60, 4600,5200,6000); // (double gaussian)
                                                        TF1 *funCo60s = fitGEL(histCo60,8400,9000);    
    TCanvas *cBACK = new TCanvas(); histBACK->Draw(); TF1 *funBACK = fitGEL(histBACK, 5600,6600);

    // We have obtained the four fitted Functions -> We now extract the middpoints
    double ADCchann[7] = {fun1Na22->GetParameter(0), fun2Na22->GetParameter(0),funCs137->GetParameter(0),funCo60->GetParameter(0),funCo60->GetParameter(7), funCo60s->GetParameter(0),funBACK->GetParameter(0)};
    double sigma[7] =  {fun1Na22->GetParameter(1), fun2Na22->GetParameter(1),funCs137->GetParameter(1),funCo60->GetParameter(1),funCo60->GetParameter(8), funCo60s->GetParameter(1),funBACK->GetParameter(2)};
    double Res[7];
    for(int ii = 0; ii < 7; ii ++){
        Res[ii] = sigma[ii] * 235.5 / ADCchann[ii];
    }
    //                    Na22    Na 22     Cs137     Co60      Co60   Co60 (sum)   BACK      (keV)
    double Energy[7] = {511.0,  1274.537, 661.657, 1173.228, 1332.492,  2505.72  ,1460.820};

    // We represent the energy vs ADC channel graph
    TGraph *grCal = new TGraph(7,ADCchann, Energy); grCal->SetTitle("Energy calibration");
    grCal->GetXaxis()->SetTitle("ADC channel");     grCal->GetYaxis()->SetTitle("Energy (keV)"); 
    TCanvas *cGr = new TCanvas(); grCal->SetMarkerStyle(5); grCal->Draw("AP"); // We draw the energy-ADC channel graph

    // function to be fitted
    TF1 *funlin0 = new TF1("funlin0",pol2,0.0,16384.0,3);
    funlin0->FixParameter(2,0); funlin0->FixParameter(2,0); // We fix the parameters so that it is a linear function

    grCal->Fit("funlin0",""); // We fit the function

    // Represent the resolution
    TGraph *grRes = new TGraph(6,Energy,Res); grRes->SetTitle("Relative Energy Resolution"); // Plot the resolution
    grRes->GetXaxis()->SetTitle("Energy (keV)");    grRes->GetYaxis()->SetTitle("Relative Energy Resolution (%)");
    TCanvas *cRes = new TCanvas(); grRes->SetMarkerStyle(5); grRes->Draw("AP");

    // We remake the calibrated histograms
    //calibrate("Histo42_45.root","pol.root",funlin0);

}


