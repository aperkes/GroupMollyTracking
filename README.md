# GroupMollyTracking

All the data and code for the paper, "Group-level phenotypes are idiosyncratic yet unpredictable over development in a clonal fish", including python code for dealing with raw tracks.

A brief summary of the relevent scripts: 

processTracks.py - takes the many track CSVs in groupTracks-filtered and calculates summary statistics. 
shuffleAllTracks.py - Runs permutation analysis for figure 2. 
overlayGroupTrack.py - plots track files over video for visualization
parse_labels.py - code to handle size data from sleap labels (requires an installation of SLEAP) 

GroupTracking.R - R Analysis code for analyzing the data
GroupTracking.Rmd - R markdown of R analysis. 
rdsFiles - R objects of various models, etc to save time when running R markdown.
