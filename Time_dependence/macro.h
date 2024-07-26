//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jul  4 10:54:05 2024 by ROOT version 6.32.02
// from TTree Data_R/CoMPASS RAW events TTree
// found on file: run01.root
//////////////////////////////////////////////////////////

#ifndef macro_h
#define macro_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector


class macro : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<UShort_t> Channel = {fReader, "Channel"};
   TTreeReaderValue<ULong64_t> Timestamp = {fReader, "Timestamp"};
   TTreeReaderValue<UShort_t> Board = {fReader, "Board"};
   TTreeReaderValue<UShort_t> Energy = {fReader, "Energy"};
   TTreeReaderValue<UInt_t> Flags = {fReader, "Flags"};


   macro(TTree * /*tree*/ =0) { }
   ~macro() override { }
   Int_t   Version() const override { return 2; }
   void    Begin(TTree *tree) override;
   void    SlaveBegin(TTree *tree) override;
   void    Init(TTree *tree) override;
   bool    Notify() override;
   bool    Process(Long64_t entry) override;
   Int_t   GetEntry(Long64_t entry, Int_t getall = 0) override { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   void    SetOption(const char *option) override { fOption = option; }
   void    SetObject(TObject *obj) override { fObject = obj; }
   void    SetInputList(TList *input) override { fInput = input; }
   TList  *GetOutputList() const override { return fOutput; }
   void    SlaveTerminate() override;
   void    Terminate() override;

   ClassDefOverride(macro,0);
   TFile *fin;
   TFile *fout;
   TH1D *HEnergy;


};

#endif

#ifdef macro_cxx
void macro::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

bool macro::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}


#endif // #ifdef macro_cxx
