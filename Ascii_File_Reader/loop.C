#include "asciiReader.h"
#include <TStopwatch.h>
#include <TFile.h>
#include <TSystemDirectory.h>
#include <string>

// This macro loops over all files in a folder (string rootfolder) and reads their data.
// Then, it creates a .root file with all the waveforms in there.
// It requires the "asciiReader.h" file.
void loop()
{ 
	TString rootfolder = "RunNUmber";
	TFile *fout = new TFile(rootfolder + ".root", "RECREATE");

	TStopwatch t; TStopwatch total; total.Start();

  	TSystemDirectory *folder = new TSystemDirectory("folder", "./" + rootfolder);
  	TList *list = folder->GetListOfFiles();
  	int N = list->GetLast(); int n = 0;
  	for(int i = 0; i <= N; i++){
	  	string name = list->At(i)->GetName(); // We obtain the filename
		if (name.size()<5) continue;  // We only consider files long enough to be named "x.dat"
		if (name.find(".dat")!= -1)  //If the file is a .dat file:
		{	
			t.Start();
			ifstream fin(rootfolder + "/" + name, ios::in);  // Open the file
			n = asciiReader(fin, fout, n); // Pass the file to the asciiReader (n keeps track of the number of events)
			fin.close();
			t.Stop();
			std::cout<<"File: "<<(i-2)<<"/"<<(N-2)<<"    #waves: "<<n<<"    "<<t.CpuTime()<<"\n";
			t.Reset();
		}
	}
	cout<<"Total time: "<<total.CpuTime();
	fout->Close();
}
