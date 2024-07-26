Explanation of all the scripts in this folder.
- fitGauss(histName, xmin, xmax):
	This script opens "Histograms.root", takes the histogram with the indicated name, and then fits a Gaussian+erc+linback in the specified ROI.

- macro6()
	This script opens "Histograms.root", and, from here, extracts the selected histograms and fits a gauss+erc+lin function with the selected ROI, which must be specified inside the script.
	Then, depending on the settings, it fits a parabole or line through the obtained data points to get a calibration curve.
	The, it runs the "calibrate.h" function with the obtained calibration function.
	As per now, it is adapted to the runs 10-13

-macro42_45() and macro46_49(), etc
	They do the same as macro6(), but they are adapted to runs 42-45, 46-49, etc respectively.

- calibrate(input, output, Tfunction)
	This function opens the input file, and loops through all histograms inside it. 
	Then, it used the function inputed to recalibrate the x axis, and saves the updates histograms in the output file.
