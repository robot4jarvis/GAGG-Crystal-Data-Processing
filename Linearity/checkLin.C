#define N 12
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
        funGEL->FixParameter(4,0);   funGEL->FixParameter(5,0); 


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

void checkLin(){
    // This macro opens a root file, reads their histograms and fits a gaussian where indicated.
    // Then, it plots the deviation of the three peaks and the chi square against the position of the first peak

    //regions of interest
double xmin[N] =  { 650,  750,  900, 1050, 1300, 1500, 1700, 1900, 2000, 2200, 2400, 2600}; // Double peak: min
double xmax[N] =  { 880, 1050, 1200, 1450, 1700, 2000, 2200, 2500, 2700, 2950, 3200, 3500}; //              max
//                      30    31    32    33    34    35    36    37    38    39    40    41    
double xmin2[N] = {1450, 1700, 2000, 2400, 2800, 3300, 3650, 4100, 4450, 4850, 5200, 5550}; // sum peak min
double xmax2[N] = {1600, 1850, 2200, 2650, 3100, 3550, 3950, 4400, 4800, 5200, 5600, 6000}; //          max
    int i = 0;
    double mean1[N];    double mean2[N];     double mean3[N];  double dev[N]; // The values we need
    double Emean1[N];   double Emean2[N];    double Emean3[N]; double Edev[N]; //their uncertainties
    double chi[N];
    // We open the Histograms file and extract all data (that we will use for calibration)
    TFile *fin = new TFile("HistogramsZoltan.root");

    for(auto k : *fin->GetListOfKeys()) {  // Loop over all objects in the file
        TKey *key = static_cast<TKey*>(k);
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TH1")) continue;  // Botam aquesta passa si no és histograma
        if (key->GetCycle() == 2) continue;
        TH1D *hist = key->ReadObject<TH1D>();

        // We now have the hist "hist"
        TF1 *funfit = fitGEL2(hist, xmin[i],(xmin[i]+xmax[i])/2,xmax[i]); // We fit a double gaussian to the two peaks.
        TF1 *funfitsum = fitGEL(hist, xmin2[i], xmax2[i]); // We fit gaussian to the sum peak.

        mean1[i]=funfit->GetParameter(0);  mean2[i] = funfit->GetParameter(7); mean3[i]=funfitsum->GetParameter(0);// And obtain their means.
        Emean1[i]=funfit->GetParError(0); Emean2[i] = funfit->GetParError(7); Emean3[i]=funfitsum->GetParError(0); // and their errors
        TCanvas *c1 = new TCanvas(); // We draw the histogram just to see that it worked.

        //hist->GetXaxis()->SetRangeUser(xmin[i], xmax2[i]);
        dev[i]=100*(mean1[i]+mean2[i]-mean3[i])/mean1[i];
        Edev[i] = 100*( (mean3[i]-mean2[i])/(mean1[i]*mean1[i])*Emean1[i] + Emean2[i]/mean1[i] +Emean3[i]/mean1[i]); // Error propagation

        //Method 2: comparing the chi square of a linear fit
        double x[3] = {mean1[i],mean2[i],mean3[i]};
        double y[3] = {1173.228, 1332.492,  2505.72};
        TGraph *gr = new TGraph(3, x,y);
        TF1 *funpol2 = new TF1("funpol2",pol2,0,16000,3);
        funpol2->FixParameter(2,0); funpol2->FixParameter(0,0); //We set the non-linear terms of the function to 0.
        gr->Fit("funpol2","Q");
        chi[i]=funpol2->GetChisquare();

        i++;// hist->Draw(); 
   }

    // We plot the deviation of the third peak
    TGraphErrors *grDev = new TGraphErrors(N,mean1, dev,Emean1,Edev); grDev->SetTitle("Deviation of the position of the sum peak");
    grDev->GetXaxis()->SetTitle("ADC channel");     grDev->GetYaxis()->SetTitle("relative dev (%)"); 
    grDev->GetYaxis()->SetRangeUser(0,10); grDev->GetXaxis()->SetRangeUser(500,3000);
    TCanvas *cGrDev = new TCanvas(); grDev->SetMarkerStyle(7); grDev->Draw("AP"); // We draw the energy-ADC channel graph

    // We plot the fitting error (chi square)
    TGraph *grC = new TGraph(N,mean1, chi); grC->SetTitle("Chi square");
    grC->GetXaxis()->SetTitle("ADC channel");     grC->GetYaxis()->SetTitle("Chi square"); 
    TCanvas *cGrC = new TCanvas(); grC->SetMarkerStyle(7); grC->Draw("AP"); // We draw the energy-ADC channel graph */
    cGrC->SetLogy();

}


