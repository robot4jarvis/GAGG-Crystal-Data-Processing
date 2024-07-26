#define N 6
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TFrame.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TRandom3.h>
#include <TBenchmark.h>
#include <TF1.h>
#include <TGraph.h>
#include <TObject.h>
#include <TGraphErrors.h>
#include <TKey.h>
//#include <TH1D.h>


double GEL(double *x, double *par) {  // gauss + erc + lin
    // Seven parameters to fit.
    double fval; double xval; xval = *x;
    fval = TMath::Gaus(xval, par[0], par[1],false) * par[2] + par[3]*(xval-par[0]) + par[4] * TMath::Erfc((xval-par[0])*par[5]) + par[6];

    return fval;

    // The function is:
    // fval = gauss(x, p0,p1) * p2 + p3*x + p4 * Erfx (x-p0)*p5 + p6
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
        funGEL->FixParameter(4,0);   funGEL->FixParameter(5,0); 


    hist->Fit("funGEL","QR+"); //The R parameter restricts the fitting to [xmin, xmax]

    return funGEL;
}

void checkTimeVariance(){

    //regions of interest
    double xmin[N] =  {  800,   800,   800,    800,    800,   800}; // Double peak: min
    double xmax[N] =  { 1050,  1050,  1050,   1050,   1050,  1050}; //              max
//                        63     64     66      67      68     69
//                      14:00  18:40  21:30   23:40   08:15  9:55
    double time[N] =  { 14.00, 18.66, 21.50,  23.66,  32.25, 33.91}; // hours
    double Etime[N]=  {  0.05,  0.05,  0.05,   0.05,   0.05,  0.05};
    int i = 0;
    double mean[N]; double sigma[N]; double Res[N];
    double Emean[N]; double Esigma[N]; double ERes[N];

    // We open the Histograms file and extract all data (that we will use for calibration)
    TFile *fin = new TFile("Histograms.root");

    for(auto k : *fin->GetListOfKeys()) {  // Loop over all objects in the file
        TKey *key = static_cast<TKey*>(k);
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TH1")) continue;  // Botam aquesta passa si no és histograma
        TH1D *hist = key->ReadObject<TH1D>();
        if(key->GetCycle()==2) continue;

        // We now have the hist "hist"
        TF1 *funfit = fitGEL(hist, xmin[i], xmax[i]); // We fit gaussian

        mean[i]=funfit->GetParameter(0); sigma[i] = funfit->GetParameter(1);
        Emean[i]=funfit->GetParError(0); Esigma[i] = funfit->GetParError(1);

        TCanvas *c1 = new TCanvas(); // We draw the histogram just to see that it worked.

        //hist->GetXaxis()->SetRangeUser(xmin[i], xmax2[i]);
        Res[i]=235.5 * sigma[i]/mean[i];
        ERes[i] = Res[i] * (Esigma[i]/sigma[i] + Emean[i]/mean[i]);

        i++;
   }

    TGraphErrors *grMean = new TGraphErrors(N,time, mean,Etime,Emean); grMean->SetTitle("Position of the Cs137 peak over time");
    grMean->GetXaxis()->SetTitle("time (h)");     grMean->GetYaxis()->SetTitle("Position (ADC channel)"); 
    TCanvas *cGrMean = new TCanvas(); grMean->SetMarkerStyle(2); grMean->Draw("AP"); // We draw the energy-ADC channel graph



    TGraphErrors *grRes = new TGraphErrors(N,time, Res,Etime,ERes); grRes->SetTitle("Energy resolution of the Cs137 peak over time");
    grRes->GetXaxis()->SetTitle("time (h)");     grRes->GetYaxis()->SetTitle("Energy resolution (%)"); 
    TCanvas *cGrRes = new TCanvas(); grRes->SetMarkerStyle(2); grRes->Draw("AP"); // We draw the energy-ADC channel graph

}

