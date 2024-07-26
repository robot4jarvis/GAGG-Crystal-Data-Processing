- How it works:
All .root files with the original data must be placed in the "root" folder (if a different folder name is needed, it has to be edited in the loop.C script).
Then, the loop.C script will be launched. It needs the macro.C and the macro.h files in the same folder.
The script will loop trough all the .root files in the "root" folder, and build a histogram using the data in the Data_F;1 tree, "Energy" leave.
All the resulting TH1 objects are saved in the "output.root" file.
