#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TGraph.h>
#include <iostream>
#include <fstream>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TKey.h>
#include <TROOT.h>
#include <TObject.h>
#include <string>
#include "fitGauss.h"
#include <TLine.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompLu.h>
#include <TTree.h>
#include <TH2D.h>

#define noEvents 1000 //15500 //55000
#define recordLength 1024
#define samplingPeriod 0.3125 // 0.3125 ns
#define trigger 0 // 0: CFD, 1: the other method
#define bas 80 // How many points do we use to determine the baseline.
#define CFD 0.5  
#define delayCH0 20 // 20
#define delayCH1 300 // 300
#define DRT 0.5
#define dump 1 // Do we want to dump the waveforms in a dump.root (Recreate) file?
#define output 0 // Do we want to ouput the histograms in a output.root (update)
using namespace std;

double analyticalZero(double xmin, double xmax, double par[4]){// analytical rootfinding for third degree polynomial
    double a = -200; double b = -200; double c = -200;
    bool d = TMath::RootsCubic(par, a, b, c);
    if(d){
        if((a>xmin)&(a<xmax)) return a;
        else{return -100; std::cout<<"ERROR";} 
    }
    else {
        if((a>xmin)&(a<xmax)) return a;
        if((b>xmin)&(b<xmax)) return b;
        if((c>xmin)&(c<xmax)) return c;
        else{return -100; std::cout<<"ERROR";} 
    }

}
double zero(double y[recordLength], int xmin, int xmax){
    // This function finds a zero in the array y, bound by [xmin, xmax]
    // It first finds a guess (xguess), by taking the absolute value of y within [xmin, xmax] and then selecting the minimum value
    // (so, the closest to the real zero). This value is called "xguess". 
    //We select xguess such that the zero is located between (xguess-1) and xguess
    // Then, we use xguess and three adjacent values to interpolate a cubic polynomial.
    // And then, we use "analyticalZero" to find the exact value of the zero.

    double aVal[recordLength]; // Array in which the absolute value of the values will be stored.
    for(int i = 0; i < recordLength; i++) aVal[i] = 10; // We inicialize all the values to 10 (a large number)
    for(int i = xmin; i <= xmax; i++) aVal[i] = abs(y[i]);  // In the ROI, we select the absolute value.

    // The first guess of the zero location is "xguess" -> The position of the point with the minimum absolute value
    int xguess = TMath::LocMin(recordLength, aVal);
    if (y[xguess] <= 0) xguess++;
    
    // We will take the 4 points closer to the guess, and store their values in "vals"
    double vals[4] = {y[xguess - 2], y[xguess-1], y[xguess], y[xguess + 1]};
    
    // third degree polynomial interpolation: a Van der Monde matrix is constructed
    TMatrixD VDM(4,4); TVectorD bvec(4);
    for (int i = 0; i < 4; i ++){
        VDM(i,0) = (xguess-2+i)*(xguess-2+i)*(xguess-2+i);
        VDM(i,1) = (xguess-2+i)*(xguess-2+i);
        VDM(i,2) = (xguess-2+i);
        VDM(i,3) = 1;
        bvec(i) = vals[i];
    }
    // We solve the system to obtain the polynomial coefficients, stored in "par[4]"
    TDecompLU lu(VDM); Bool_t ok;
    TVectorD x = lu.Solve(bvec, ok);
    double par[4] = {x(3), x(2), x(1), x(0)};
    
    return analyticalZero(xguess-2, xguess +1, par); // An analytical method (pol3)
}
double sum(double y[recordLength], int xmin = 0, int xmax = recordLength){
    double res = 0;
    //xmax = TMath::Min(500, int(TMath::LocMax(recordLength, y)));  // To integrate only up to the peak
    for(int i = xmin; i < xmax; i ++){
        res += y[i];
    }
    return res;
}

void filterSpectrum() {

    TString fileName = "Run23";
    //TString title = "Run11 CFD 300 300 0.5 (Pulsar)";
    TString title = "Run23 EDR 0.5";
    // The .root input file will be taken from the "ROOTFILES" folder.
    TFile *fin = new TFile("ROOTFILES/" + fileName + ".root");
    TFile *fout= new TFile(fileName +"_Definitive output.root", "UPDATE");
    TFile *fDump = new TFile("dump.root", "RECREATE");

    TH1D *hist0 = new TH1D("Channel 0", "Amplitude spectrum CH0 " + title,2000, 0,400); 
          hist0->GetXaxis()->SetTitle("Integral (V*ns) "); hist0->GetYaxis()->SetTitle("Counts");
    TH1D *hist1 = new TH1D("Channel 1", "Amplitude spectrum CH1 " + title,2000, 0,400);
          hist1->GetXaxis()->SetTitle("Integral (V*ns)"); hist1->GetYaxis()->SetTitle("Counts");
    TH1D *histDif = new TH1D("Time difference", "Time difference " + title,2048*25, -1024*samplingPeriod,1024*samplingPeriod); 
          histDif->GetXaxis()->SetTitle("Time difference (ns)"); histDif->GetYaxis()->SetTitle("Counts");
    
    TH2D *histCol = new TH2D("ColourGraph", "Two channel graph" + title, 200, 0, 50, 200, 0, 350);
          histCol->GetXaxis()->SetTitle("Channel 0"); histCol->GetYaxis()->SetTitle("Channel 1"); 

    TF1 *fun = new TF1("fun","pol0",0,200);
    TStopwatch t; t.Start();
    double peak0; double pos0; double bas0; double posDif;
    double peak1; double pos1; double bas1;
    int j = 0;

    for(auto k : *fin->GetListOfKeys()) {  // Loop over all objects

        // We get each object in the file
        TKey *key = static_cast<TKey*>(k);
        TClass *cl = gROOT->GetClass(key->GetClassName());
        string className = key->GetClassName();
        if (className.compare("TObjArray")) continue;  // Botam aquesta passa si no és histograma


        // We recover the graphs from the object array
        TObjArray *col = key->ReadObject<TObjArray>();
        TGraph *gr0 = (TGraph*)col->At(0);  // We now have recovered both graphs
        TGraph *gr1 = (TGraph*)col->At(1);  // We now have recovered both graphs

        // We obtain the baselines
        fun->SetRange(0,bas); gr1->Fit(fun,"RQ"); bas1 = fun->GetParameter(0);
        fun->SetRange(0,bas); gr0->Fit(fun,"RQ"); bas0 = fun->GetParameter(0);

        // We recover the arrays
        double val0[recordLength]; double val1[recordLength];
        for(int i = 0; i < recordLength; i++){
            val0[i] =bas0 - gr0->GetPointY(i);
            val1[i] =bas1 - gr1->GetPointY(i);
        }


        //             PEAK DETERMINATION
        peak0 = sum(val0, 0, 400); peak1 = sum(val1); // Integral

        //             FILTER
        //if ((peak0 < 30)||(peak0 > 55)){j++; continue;} // We skip this step if the peak is not in our range
        //if ((peak1 < 110)||(peak1 > 145)){j++; continue;} // We skip this step if the peak is not in our range

        // ========================================================================================
        if(trigger == 0){
        //             CFD 
        for(int i = 0; i < (recordLength-TMath::Max(delayCH0, delayCH1)-10); i++){
            val0[i] = val0[i] - CFD*val0[i + delayCH0];
            val1[i] = val1[i] - CFD*val1[i + delayCH1];
        }
        pos0 = zero(val0, TMath::LocMin(recordLength, val0), TMath::LocMax(recordLength, val0)); 
        pos1 = zero(val1, TMath::LocMin(recordLength, val1), TMath::LocMax(recordLength, val1)); 
        }
        else{
        //      Dynamic Rising Edge Trigger
        double p0 = TMath::MaxElement(recordLength, val0); double p1 = TMath::MaxElement(recordLength, val1);
        for(int i = 0; i < recordLength; i++){
            val0[i] = val0[i] - DRT*p0;
            val1[i] = val1[i] - DRT*p1;
        }
        pos0 = zero(val0, 0, TMath::LocMax(recordLength, val0)); 
        pos1 = zero(val1, 0, TMath::LocMax(recordLength, val1)); 
        }
       // =======================================================================================

        //          ZEROS
        // We determine the zeros for each wave, and fill the histogram
        if((pos0 < 0)||(pos1 < 0)){j++; continue;} // if the rootfinding wasn't succesfull, we filter it out
        posDif = (pos1-pos0);

        // We fill the histograms. We have run all calculations with divisions and not ns, so we multiply by "samplingPeriod".
        hist0->Fill(peak0); hist1->Fill(peak1); 
        histCol->Fill(peak0, peak1); 
        histDif->Fill(posDif*samplingPeriod);
                
        // Dump
        if(dump){
            TMultiGraph *mg = new TMultiGraph(); 
            double xline[2] = {pos0,pos1}; double yline[2] = {0,0}; TGraph *line = new TGraph(2,xline,yline);
            TGraph *grmod0 = new TGraph(recordLength, val0); TGraph *grmod1 = new TGraph(recordLength, val1);
            mg->Add(grmod0); mg->Add(grmod1); mg->Add(line);
            TString name = "waves" + to_string(j); fDump->WriteObject(mg, name);}

        if (j>noEvents) break;
        cout<<100 * j/noEvents<<"%\r"<<flush; j++;
    }

    t.Stop(); t.Print();
    TCanvas *c0 = new TCanvas("hist0","hist0"); hist0->Draw(); 
    TCanvas *c1 = new TCanvas("hist1","hist1"); hist1->Draw();
    TCanvas *c2 = new TCanvas("Combined graph", "CombinedGraph"); histCol->Draw("colz");
    fitGauss(histDif, -10, 10);
    cout<<"Efficienty (events past filter / total events): "<<100 * histDif->GetEntries()/ j <<"%\n";

    TObjArray *arr = new TObjArray();
    arr->Add(hist0); arr->Add(hist1); arr->Add(histDif); arr->Add(histCol);
    if(output) fout->WriteObject(arr, title);

    //fDump->Close()
}


