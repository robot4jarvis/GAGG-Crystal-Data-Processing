#include "calibrate.h"
#include <TCanvas.h>
#include <TMath.h>
#include <TH1D.h>
#include <TF1.h>
#include <TGraph.h>
#include <TFile.h>
#include <TGraphErrors.h>

#define N 8

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
TF1 *fitGEL(TH1D *hist, double xmin, double xmax, bool erc = true){

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
        if(!erc) {funGEL->FixParameter(4,0);};
    hist->Fit("funGEL","QR+"); //The R parameter restricts the fitting to [xmin, xmax]

    return funGEL;
}

TF1 *fitGEL2(TH1D *hist, double xmin, double xmid, double xmax, bool erc = true){

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

        if(!erc) {funGEL2->FixParameter(4,0);}
    hist->Fit("funGEL2","QR+"); //The R parameter restricts the fitting to [xmin, xmax]

    return funGEL2;
}

void macro59_62(){
    // This macro opens a root file, reads their histograms and fits a gaussian where indicated.
    // Then, using calibrate.h, it calibrates all histograms in the indicated file.
    // Rigth now, it is set up to use the runs 46-49
    TString fnamein = "Histo59_62";
    TString fnameout = fnamein + "_Cal";
    // We open the Histograms file and extract all data (that we will use for calibration)
    TFile *f = new TFile(fnamein + ".root");
    TH1D *histNa22 = (TH1D*)f->Get("run60.root");  histNa22->SetTitle("Na22");
    TH1D *histCs137 = (TH1D*)f->Get("run62.root"); histCs137->SetTitle("Cs137");
    TH1D *histCo60 = (TH1D*)f->Get("run61.root");  histCo60->SetTitle("Co60");
    TH1D *histMn54 = (TH1D*)f->Get("run59.root");  histMn54->SetTitle("Mn54");


    // We select all the regions of interest
    TCanvas *cNa22 = new TCanvas();  histNa22->Draw();   TF1 *fun1Na22 = fitGEL(histNa22, 630,780);
                                                         TF1 *fun2Na22 = fitGEL(histNa22, 1650,1950);
                                                         TF1 *fun3Na22 = fitGEL(histNa22, 2300, 2650); // sum peak
    TCanvas *cCs137 = new TCanvas(); histCs137->Draw();  TF1 *funCs137 = fitGEL(histCs137, 800,1000);
    TCanvas *cCo60 = new TCanvas();  histCo60->Draw();   TF1 *funCo60 = fitGEL2(histCo60, 1500,1800,2100); // (double gaussian)
                                                         TF1 *funCo60s = fitGEL(histCo60,3300,3650, false);   // sum peak 
    TCanvas *cMn54 = new TCanvas();  histMn54->Draw();   TF1 *funMn54 = fitGEL(histMn54, 1050,1300);


    // We have obtained the four fitted Functions -> We now extract the middpoints and sigmas
    double ADCchann[N] = {fun1Na22->GetParameter(0), fun2Na22->GetParameter(0), fun3Na22->GetParameter(0), funCs137->GetParameter(0),funCo60->GetParameter(0),funCo60->GetParameter(7), funCo60s->GetParameter(0),funMn54->GetParameter(0)};
    double EADCchann[N] = {fun1Na22->GetParError(0), fun2Na22->GetParError(0), fun3Na22->GetParError(0), funCs137->GetParError(0),funCo60->GetParError(0),funCo60->GetParError(7), funCo60s->GetParError(0),funMn54->GetParError(0)};
    double sigma[N] =  {fun1Na22->GetParameter(1), fun2Na22->GetParameter(1),fun3Na22->GetParameter(1),funCs137->GetParameter(1),funCo60->GetParameter(1),funCo60->GetParameter(8), funCo60s->GetParameter(1),funMn54->GetParameter(1)};
    double Esigma[N] =  {fun1Na22->GetParError(1), fun2Na22->GetParError(1),fun3Na22->GetParError(1),funCs137->GetParError(1),funCo60->GetParError(1),funCo60->GetParError(8), funCo60s->GetParError(1),funMn54->GetParError(1)};

    double Res[N]; double ERes[N];
    for(int ii = 0; ii < N; ii ++){
        Res[ii] = sigma[ii] * 235.5 / ADCchann[ii];
        ERes[ii] = Res[ii] * (Esigma[ii]/sigma[ii] + EADCchann[ii]/ADCchann[ii]);
    }
    //                    Na22    Na 22     Na22(s), Cs137     Co60      Co60    Co60 (sum)   Mn54      BACK(K40)   (keV)
    double Energy[N] = {511.0,  1274.537, 1785.537, 661.657, 1173.228, 1332.492,  2505.72  , 834.848};
    double EEnergy[N] ={   0.0,     0.007,   0.007,   0.003,    0.003,    0.004,    0.007,     0.003};

    // We represent the energy vs ADC channel graph
    TGraphErrors *grCal = new TGraphErrors(N,ADCchann, Energy, EADCchann, EEnergy); grCal->SetTitle("Energy calibration");
    grCal->GetXaxis()->SetTitle("ADC channel");     grCal->GetYaxis()->SetTitle("Energy (keV)"); 
    TCanvas *cCal = new TCanvas(); grCal->SetMarkerStyle(2); grCal->Draw("AP"); // We draw the energy-ADC channel graph

    // function to be fitted
    TF1 *funlin0 = new TF1("funlin0",pol2,0.0,16384.0,3);
    funlin0->FixParameter(0,0); funlin0->FixParameter(2,0); // We fix the parameters so that it is a linear function

    grCal->Fit("funlin0","W"); // We fit the function

    // Represent the resolution
    TGraphErrors *grRes = new TGraphErrors(N,Energy,Res, EEnergy, ERes); grRes->SetTitle("Relative Energy Resolution"); // Plot the resolution
    grRes->GetXaxis()->SetTitle("Energy (keV))");    grRes->GetYaxis()->SetTitle("Relative Energy Resolution (%)");
    TCanvas *cRes = new TCanvas(); grRes->SetMarkerStyle(2); grRes->Draw("AP");

// Peak to Area ratio:
    double a = 100 / (funlin0->GetParameter(1)); double b = 170 * a;
    double peakArea[4] = {fun1Na22->GetParameter(2) * fun1Na22->GetParameter(1) * TMath::Sqrt(2*3.14159265359),
    funCs137->GetParameter(2) * funCs137->GetParameter(1) * TMath::Sqrt(2*3.14159265359),
    funMn54->GetParameter(2) * funMn54->GetParameter(1) * TMath::Sqrt(2*3.14159265359),
    (funCo60->GetParameter(2) * funCo60->GetParameter(1) + funCo60->GetParameter(9) * funCo60->GetParameter(8) ) *TMath::Sqrt(2*3.14159265359)};
    double totalArea[4] = {histNa22->Integral(a,b), histCs137->Integral(a,b), histMn54->Integral(a,b), histCo60->Integral(a,b)};
    double PTT[4];
    for(int i = 0; i < 4; i ++) PTT[i] = peakArea[i]/totalArea[i];
    double labels[4] = {511.0,661.657, 834.848, 1200};
    TGraph *grPTT = new TGraph(4,labels, PTT); grPTT->SetTitle("Peak to Total ratio (PTT)"); // Plot the resolution
    grPTT->GetXaxis()->SetTitle("Energy (keV))");    grPTT->GetYaxis()->SetTitle("Peak to Total ratio (PTT)");
    TCanvas *cPTT = new TCanvas(); grPTT->SetMarkerStyle(2); grPTT->Draw("AP");


    // We remake the calibrated histograms (REWRITES FILE)
    //calibrate(fnamein + ".root",fnameout + ".root",funlin0);

    TFile *fout = new TFile("Graphs.root", "UPDATE");
    fout->WriteObject(grCal,"CalibrationGraph");
    fout->WriteObject(grRes,"ResolutionGraph");
    fout->WriteObject(grPTT, "PTTGraph");
    fout->Close();

}


