This code serves to check the "linearity" of the sensor.
To that extent, we used Co60, which should yield three peaks: two of them close together (1173.228 keV and  1332.492 keV), and their sum
(2505.72 keV). If the sensor is perfectly linear, then the third peak should be at the sum of the two previous peaks. Furthermore, if the sensor
is perfectly linear, then we should expect a very low chi^2 for a linear fit of the three peaks.

So, this folder contains four scripts:
    - loop.C, macro.C and macro.h:
    These scripts take all .root files in the "rootfiles" folder, and convert their data into histograms which are saved in a "Histograms.root" file,
    so that they can be later processed. It's the same thing that the folder "Energy_Histogram_maker" does.
    Right now, the "rootfiles" folder is empty, since I didn't want to surpass GitHub's size limit. But the .root files with the resulting energy histograms
    are available in this same folder.

    - checkLin.C
    This script reads all histograms in "Histograms.root", and fits a double gaussian to the first two peaks and a single gaussian to the third peak.
    Then, it calculates the relative deviation of the third peak and chi^2, and plots them.
    The ROI of the double peak and the third peak must be manually entered in the checkLin.C script.

    The numer of histograms must be specified at the start of the file (#define N 11).


======= ROI for the "Joan" runs:
double xmin[N] = {  700, 1200, 1800, 2200, 2700, 3200, 3600, 1600, 1900, 1950, 2100}; // Double peak: min
double xmax[N] = { 1000, 1800, 2500, 3100, 3700, 4400, 5000, 2200, 2600, 2800, 2900}; //              max
//                  18    19    20    21    22     23   24    25    26    27    29
double xmin2[N] = {1550, 2850, 4000, 5000, 5850, 6700, 7300, 3550, 4200, 4450, 4600}; // sum peak min
double xmax2[N] = {1860, 3100, 4400, 5500, 6300, 7200, 8000, 4000, 4800, 4900, 5200}; //          max

======= ROI for the "Zoltan" runs:

double xmin[N] =  { 650,  750,  900, 1050, 1300, 1500, 1700, 1900, 2000, 2200, 2400, 2600}; // Double peak: min
double xmax[N] =  { 880, 1050, 1200, 1450, 1700, 2000, 2200, 2500, 2700, 2950, 3200, 3500}; //              max
//                      30    31    32    33    34    35    36    37    38    39    40    41    
double xmin2[N] = {1450, 1700, 2000, 2400, 2800, 3300, 3650, 4100, 4450, 4850, 5200, 5550}; // sum peak min
double xmax2[N] = {1600, 1850, 2200, 2650, 3100, 3550, 3950, 4400, 4800, 5200, 5600, 6000}; //          max