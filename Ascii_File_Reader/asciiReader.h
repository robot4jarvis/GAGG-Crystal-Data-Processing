#include <TFile.h>
#include <TObjArray.h>
#include <TH1D.h>
#include <TF1.h>
#include <TGraph.h>
#include <iostream>
#include <fstream>
#include <TMultiGraph.h>
#include <TCanvas.h>
#define noEvents 500 // number of events in each file
#define recordLength 1024
using namespace std;

// This function turns string with coma separeted doubles into a C++ array
std::array<double,recordLength> unStringCSV(string xString){
    std::array<double,recordLength> x{0};
    char c; string val; int j = 0; double xval;
    for(int ii = 0; ii < xString.length(); ii ++){
        c = xString[ii];
        if((c == ' ') || (ii > xString.length()-1)){
            x[j] = stof(val); val = "";
            j++;
        }
        else val = val + c;
    }
    return x;
}

//This function reads a single event from the imputfile, and returns a TObjectArray containing two graphs
TObjArray *twographsArray(ifstream& fin, int j){
        string buffer;
        for(int i = 0; i < 4; i++){ std::getline(fin, buffer);} // We skip three lines and get the forth (they contain a header for each event)
        string line1 = buffer; // We get the line containing CH0 info
        std::getline(fin,buffer); std::getline(fin,buffer); // Skip an empty line and read the next one
        string line2 = buffer;

        std::array<double,recordLength> arr1 = unStringCSV(line1); // Use the previoysly defined function to turn the line into a double array
        std::array<double,recordLength> arr2 = unStringCSV(line2);

        // ROOT is not very compatible with C++ Arrays, so we must turn the c++ array into a C style array
        double x[recordLength]; double val1[recordLength]; double val2[recordLength];
        for(int i = 0; i < recordLength; i++)
            {	val1[i] = arr1[i]; val2[i] = arr2[i];
                x[i] = i*1.0;
            }

        // We build two ROOT graphs with the arrays
        TGraph *gr1 = new TGraph(recordLength, x, val1);
        TGraph *gr2 = new TGraph(recordLength, x, val2);

        // We create a TObjArray and add the two graphs
        TObjArray *col = new TObjArray();
        col->Add(gr1); col->Add(gr2);       
        cout<<100*j/noEvents <<"% \r"<<flush;
        return col;
}

// This program reads a binary file with two channels
int asciiReader(ifstream& fin, TFile *fout, int n) {
    string line;
    for(int i = 0; i < 4; i++){ // Skip the first four lines
        std::getline(fin, line);
    }

    for (int j = 0; j < noEvents; j++){ // For every event, we run the "twographsarray" function
        TString name = "Wave " + to_string(n);
        fout->WriteObject(twographsArray(fin, j),name);
        n++;
    }
    return n;
}

