A Breaf Description of the Scripts:

# 1. wrapup-generate-comp-2nd.R, generate-ccf-2nd.R
These are code to generate comparison data on the server, and extract ccf. Unless  you would like to know how those data is generated, they are never drictly related to other code. There are other operations done in the server, usually data generation for data in organized data, e.g. ccf_nonmatched_km.rds etc. are not included here.

# 2. helper functions.R
Many helper functions, some of which are used from time to time in other scripts if there are functions not defined in some scripts.

# 3. model-full-data-ccf.R
To fit the beta mix models. This is the basis of any following discussion.

# 4. error rate.R 
Corresponding to the theoretical error rate section of the writeup, based on the models/distributions estimated in #3

# 5. application.R
The application section in writeup

# 6. changing_size.R (under updating)
The changing size/sample size/sensitivity section in the writeup

# 7. wrapup-generate-nonmatched-bullets-comp.R is for non-matching bullets comparisons

# 8. generate-ccf-from-knm-bullets-comp.R is to extract CCFs from non-matching bullets comparisons

# bind-ccvs-to-csv.R
To get ccf_data.csv, the final data we will use directly

