# GAGG-Crystal-Data-Processing
Hello!
In this repository, I have uploaded all the code (at least, the usable code) that I have used to process the data obtained from a GAGG scintillator detector.
Everything in here was done using C++ and ROOT (root.cern).
## Folders

### Ascii_File_Reader:
This folder serves to turn the ascii files outputed by CAEN's WaveCatcher into much more convenient and light .root files, where each waveform captured by the digitizer corresponds to a TGraph. WaveCatcher usually splits the data into many files, which must be placed in a folder within "Ascii_File_Reader".

When one launches the loop.C macro, it will take all files in said folder (which must be specified inside the script) and use the function asciiReader() (defined in asciiReader.h) to read them and generate a single .root file with all the graphs. In this .root file one can find one TObjArray for every event. This TObjArray contains two TGraph objects, each one corresponding to one channel.

The asciiReader.h function can only read files with two channels. It could be modified to read inputs with more or less channels.

### AS
    
