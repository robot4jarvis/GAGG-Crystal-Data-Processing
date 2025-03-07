# GAGG-Crystal-Data-Processing
Hello!

In this repository, I have uploaded all the code (at least, the usable code) that I have used to process the data obtained from a GAGG scintillator detector. I have removed practically all the input data and most of the output to have a cleaner repository.

Everything in here was done using C++ and ROOT (root.cern).
## Time resolution analysis
### Ascii_File_Reader:
This folder serves to turn the ascii files outputed by CAEN's WaveCatcher into much more convenient and light .root files, where each waveform captured by the digitizer corresponds to a TGraph. WaveCatcher usually splits the data into many files, which must be placed in a folder within "Ascii_File_Reader".

When one launches the loop.C macro, it will take all files in said folder (which must be specified inside the script) and use the function asciiReader() (defined in asciiReader.h) to read them and generate a single .root file with all the graphs. In this .root file one can find one TObjArray for every event. This TObjArray contains two TGraph objects, each one corresponding to one channel.

The asciiReader.h function can only read files with two channels. It could be modified to read inputs with more or less channels.

### Waveform_analysis
In this folder one can find all the necessary macros to read the .root files generated by the previous folder. By placing them in the "ROOTFILES" folder (now empty), and indicating its names inside the filterSpectra.C macro, these waveforms will be read and processed.
**Inside the "Interesting histograms" file there are, as the name implies, some interesting histograms. Those are the relevant runs where the time resolution figures were obtained from. Check the *info.xlsx* excel file (it has two pages) to know what run corresponds to which settings.**

## Energy histogram analysis
### Energy_Histogram_Maker
An event list in .root format was generated using the COMPASS software and a digitizer. To extract a Energy histogram from there, all .root files must be placed in the "root" folder. One energy histogram will be generated for each file (or several, if a file contains more than one cycle). Then, the loop.C script can be used be launched. It needs the macro.C and the macro.h files in the same folder.

The script will loop trough all the .root files in the "root" folder, and build a histogram using the data in the Data_F;1 tree, "Energy" leave.

All the resulting TH1D objects are saved in the "output.root" file. You can then process this "output.root" folder with the following macros.

### Energy_Spectra_Processing
These are the macros used to obtain the energy calibration graphs. Each one of the macros have been tuned to a particular set of runs, however, they contain virtually the same code; the only diffence are the values of some constants.
Therefore, if you need to analyze a different run, you can simply copy and edit one of the macros.

They take, as an input, the files generated in "Energy_Histogram_Maker". You must specify the file name inside the corresponding macro. It extracts the selected histograms and fits a gauss+erc+lin function with the selected ROI, to each peak. Then, it will generate a Resolution (FWHM) - Energy and a Position - Energy ("calibration line"), and fit a line to find the calibration function. It also calculates the "Peak to Total" ratio (PTT). It will save those graphs to a "xxx_Cal" file.

**to be finished... or not**
