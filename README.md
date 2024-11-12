itrSurv is an R package to estimate optimal individualized treatment regimes for various survival endpoints. The authors introduce an multi-utility value function that determine the optimal rule via prioritizing survival while also accounting for a secondary endpoint.

# ** IMPORTANT NOTE: THIS PACKAGE IS COMPLETE FOR "CR" ENDPOINT BUT IS ON-GOING DEVELOPMENT FOR "RE" ENDPOINT, and will continue to be in the midst of updating until this note is removed (estimated - December 2024). **

itrSurv works for multiple treatment options, but is only applicable to the single stage disease setting. Please refer to our papers for technical details.

Currently, itrSurv is only applicable to competing risks with a priority cause of interest (as of May 2024), and recurrent events (on-going, as of July 2024).

This is the doctoral dissertation thesis work of Christina W Zhou under the mentorship of Michael R Kosorok at UNC-CH Department of Biostatistics.


---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Updated July 2024:
Now incorporating a second endpoint, recurrent event. NOTE: data must be properly inputted before running itrSurv function!! Please read the documentation!!

Updated 1/17/24:

To install in cluster (linux):
First, build source .tar.gz file from .RProj on local machine, then copy to cluster location. Then, unzip the file by doing the following:
1) go to the folder that has the R scripts you want to run
2) open R by typing R (make sure you module load r/X.X.X the right version)
3) type: install.packages('[location of .tar.gz file]/itrSurv_0.1.0.tar.gz', repos = NULL, type='source')
4) library(itrSurv)
5) exit and run bash CR_S2run.sh

! If package is ONLY binary and not source, then go to DESCRIPTION and delete the line with "BUILT" and "PACKAGED"
