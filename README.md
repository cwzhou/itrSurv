## itrSurv

**itrSurv** is an R package designed to **estimate optimal individualized treatment regimes for various survival endpoints**. This package provides tools for **survival data involving competing risks or recurrent events**. Please use as appropriate for your data. The authors introduce a joint value function that determine the optimal rule via prioritizing survival while also accounting for a secondary endpoint. itrSurv works for multiple treatment options, but is only applicable to the single stage disease setting. Please refer to our papers for technical details.

This is the doctoral dissertation thesis work of Christina W Zhou under the mentorship of Michael R Kosorok at UNC-CH Department of Biostatistics.

### Installation

You can install the package directly from GitHub:
```r
# install.packages("devtools")
devtools::install_github("cwzhou/itrSurv")
```

If you encounter any issues, please refer to the Installation Guide at the end of this README. Because the package contains C, C++, and Fortran code, appropriate compiler tools are required to successfully build and install it.

---

### Useful Code

- **Simulation Studies**: R code for simulation studies is located at `https://github.com/cwzhou/Analyses/tree/main/Simulations`. See Paper1_CR for Competing Risk simulation studies and Paper3_RE for Recurrent Event simulation studies.
- **University Cohort Application**: R code for applying the methodology to the university cohort (not available for public use) can be found at `https://github.com/cwzhou/Analyses/tree/main/RDA/Paper1_CR`. There is a sensitivity simulation analysis that can be used instead, set revision = 1 for Paper1_CR.

Currently, itrSurv is applicable to settings with either 1) competing risks with a priority cause of interest (as of April 2024), or 2) recurrent events with terminal events (as of January 2025).

***

### Updates

- **January 2025**: Finished updates for developing RCIF and the itrSurv estimator for RE setting.
- **June 2024**: Now incorporating a second endpoint, recurrent events in the recurrent with terminal events setting. NOTE: data must be properly inputted before running itrSurv function!! Please read the documentation!!
- **1/17/24**:

***

### INSTALLATION GUIDE

#### Installing `itrSurv` on macOS
This guide walks you through installing the `itrSurv` R package from GitHub on macOS.  
It includes installing Homebrew, GCC, configuring R to use the correct compiler, and installing the package.

#### 1. Install Homebrew (if not already installed)
Homebrew is the package manager for macOS. Open Terminal and run:
```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```
Verify installation: ```brew --version```

#### 2. Install GCC via Homebrew
itrSurv requires compilation, so you need GCC: ```brew install gcc```

Check the installed GCC version: ```brew list gcc```

You should see something like:
```bash 
/usr/local/Cellar/gcc/15.1.0/bin/gcc-15
/usr/local/Cellar/gcc/15.1.0/bin/g++-15
/usr/local/Cellar/gcc/15.1.0/bin/gfortran-15
...
```

#### 3. Confirm Command Line Tools for Xcode
Ensure Xcode command line tools are installed: ```xcode-select --install```

#### 4. Install itrSurv from GitHub
Finally, open R and install the package:
```r 
# install.packages("remotes")
library(remotes)
remotes::install_github("cwzhou/itrSurv")
library(itrSurv)
```

#### 5. Troubleshooting

If you encounter errors while installing the package, try the tips below (as needed) and then re-run Step 4.

#### Configure R to use Homebrew GCC
Create or edit the R Makevars file:
```bash
mkdir -p ~/.R
nano ~/.R/Makevars
```
Add the following lines (adjust paths if your GCC version is different):
```make 
CC=/usr/local/Cellar/gcc/15.1.0/bin/gcc-15
CXX=/usr/local/Cellar/gcc/15.1.0/bin/g++-15
FC=/usr/local/Cellar/gcc/15.1.0/bin/gfortran-15
F77=/usr/local/Cellar/gcc/15.1.0/bin/gfortran-15
```
Save and exit (Ctrl+O, Enter, Ctrl+X). Check your exact path with: ```which gcc-15```.

NOTE: make sure this path is the same as found in ```brew list gcc```!

#### Verify Build Tools in R
Restart R or RStudio and run: `pkgbuild::check_build_tools(debug = TRUE)`. If everything is configured correctly, it should compile a simple C file without errors.

#### Other troubleshooting
Error: gcc-13: command not found → Make sure ~/.R/Makevars points to your installed GCC version.

Permission issues → Use R as a user, not sudo.

Failed to compile C/C++ files → Verify Xcode command line tools are installed and up to date.

#### Using .tar.gz 
To install in cluster (linux):
First, build source .tar.gz file from .RProj on local machine, then copy to cluster location. Then, unzip the file by doing the following:
1) go to the folder that has the R scripts you want to run
2) open R by typing R (make sure you module load r/4.1.2 the right version)
3) type: install.packages('[location of .tar.gz file]/itrSurv_0.1.0.tar.gz', repos = NULL, type='source')
4) library(itrSurv)
5) exit and run bash CR_S2run.sh

! If package is ONLY binary and not source, then go to DESCRIPTION and delete the line with "BUILT" and "PACKAGED"
