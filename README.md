# SAUSI_paper
All files and code associated with SAUSI paper

z-score_regression_analysis:
Use: run multi-logistic and multi-linear regression on normalized data-sets. Find training accuracy & plot feature importance.
Input Data: a .csv file. each column is a behavior or group. each row is a different mouse/trial. 
The data filled in are the z-score values for each mouse with each behavior. Or "0" or "1" under the groups column to designate "control" or "experimental" for the logistic regression.
The first row contains column headers.

deepSqueak_summary_data:
Use: Calculates total number of USVs, average, and median values for USV attributes identified in deepSqueak for each mouse.
It compiles all of these summary values into a spreadsheet.
Input Data: .xlsx files exported directly from DeepSqueak for each mouse located in a single folder. These contain detailed data for each USV identified.

sleap_analysis:
Use: calculates distances between mice (including nose-nose distance). Also has uses for interpolation and smoothing of data.
Input Data: .csv files exported directly from SLEAP fofr each mouse located in a single folder. These contain the x-y locations body points in every frame for two mice tracked in SLEAP.

behavior_rastor_plot:
Use: Creates a rastor plot overlaying the active behaviors with nose-nose euclidean distance over time for a single mouse/trial. 
Input Data: a .csv file containing a single column with nose-nose distance data for a single mouse/trial (generated by sleap_analysis). 
A .csv file containing frame-by-frame behavior data for a single mouse/trial. First row is column titles (with one column for each behavior). The data in each column are 0 or 1 to indicate whether that behavior took place in that frame (frame indexed by row).
This behavior data can be generated by programs such as EthoVision by setting time bins to every frame and exporting a single trial. Make sure to delete all other descriptions generated in the ethovision file. Note the rate at which ethovision tracked each video - this informatio is used in the code.

motionmapperpy_modifiedForSLEAP:
Use: runs motionmapperpy pipeline for unsupervised machine learning analysis. Modified to work with data generated from SLEAP and other usability.
Input Data: video files (.mp4) and SLEAP-tracked files (.csv) with matching names in a single folder. 
See https://github.com/bermanlabemory/motionmapperpy for original code, use, and instillation instructions.

questions: submit any questions in the issues tab or email jordan.grammer@utah.edu
