void calibrate(TString inputfile, TString outputfile, TF1 *calFun){// Given an imput file with histograms, an output file (which will be overwritten), and a function, recalibrates ALL histograms in imput and puts them in the output.
    // TString inputfile, TString outputfile, TF1 calFun
    // We open the corresponding files
    TFile *fin= new TFile(inputfile);
    TFile *fout= new TFile(outputfile, "RECREATE");

    for(auto k : *fin->GetListOfKeys()) {  // Loop over all objects
        TKey *key = static_cast<TKey*>(k);
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TH1")) continue;  // Botam aquesta passa si no és histograma
        TH1 *hist = key->ReadObject<TH1>();

        // We now have the hist "hist"
        TString oldTitle = hist->GetName();
        TString newTitle = oldTitle + "Cal";
        int N = hist->GetNbinsX();
        TH1D *newhist = new TH1D(newTitle,newTitle,N,0,calFun->Eval(N));

        double xlow[N+1];
        for(int ii = 0; ii < N; ii++){
            xlow[ii] = calFun->Eval(ii-0.5);
            newhist->SetBinContent(ii,hist->GetBinContent(ii));
        }
        xlow[N]=calFun->Eval(N+0.5);
        newhist->GetXaxis()->Set(N,xlow);

        fout->WriteObject(newhist, newTitle);
   }
   fin->Close(); fout->Close();
}