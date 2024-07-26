#include <string>
void loop()
{ 
	TString rootfolder = "root"; // The name of the folder where the .root files are located
	
  	TSystemDirectory *folder = new TSystemDirectory("folder", "./" + rootfolder);
  	//folder->GetListOfFiles()->Print();
  	TList *list = folder->GetListOfFiles();
  	int N = list->GetLast();

  	for(int i = 0; i <= N; i++){
	  	string name = list->At(i)->GetName();
		if (name.size()<6) continue;  // We only consider files long enough to be named "x.root"
		if (name.compare(name.size()-5,5,".root")==0) //If the file is a .root file:
		{
			TFile *fin = new TFile("./"+ rootfolder + "/" + name);

			for(auto k : *fin->GetListOfKeys()) {  // Loop over all objects
        		TKey *key = static_cast<TKey*>(k);
       			TClass *cl = gROOT->GetClass(key->GetClassName());
        		if (!cl->InheritsFrom("TTree")) continue;  // We skip if is not a TTree
				TTree* t = key->ReadObject<TTree>();
				TString Tname(name);
				t->Process("macro.C", Tname);
   			}
			fin->Close();
		}
	}
}
