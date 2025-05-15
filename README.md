This repository contains all r code for a stable isotope mixing model termed the Organic Matter Supply Model.
The model description is currently submitted and being considered for publication.

Follow these steps to run the code.
1. Ensure R Studio is properly installed on your machine.
2. Download and unzip the file directory in this repository.
3. Open the R Project file titled "Organic-Matter-Supply-Model.Rproj".
4. To insure all necessary r packages are installed, navigate to the console and run the following command:
      > source("install_packages.R")
      > # The necessary packages to knit any of the .rmd files should now be installed.
5. Open and run/knit any of the included .rmd files.

OMSM_ALOHA_sim-zoops....rmd files contain OMSM model runs using simulated zooplankton data based on AA-CSIA data from particles collected at Station ALOHA. These documents rely on many functions stored in the "Functions" and "Utilities" folders.
- The base version treats Phe as a conservative tracer, though Phe fractionates during zooplankton data simulation.
- The _frac-Phe varsion treats Phe as a fractionating tracer.
- The _no-Thr version excludes Thr as a mixing tracer.

OMSM_TEMPLATE_sim-zoops.rmd provides a template for users to test the OMSM on simulated consumers generated using their own data.
OMSM_Data_Template.xlsx in the Data folder provides a template for users to input their own data for use with the OMSM.

D15N_Regression_Analysis.rmd assesses relationships between d15N values in amino acids with trophic position at a variety of locations.
Source-Separation_ALOHA-All.rmd assess isotopic separation of organic matter sources at Station ALOHA.

AA-CSIA_....xlsx data files are all data from previously published studies.
- AA-CSIA_ALOHA.xlsx, AA-CSIA_ALOHA_Summer.xlsx, and AA-CSIA_ALOHA_Winter.xlsx contain AA-CSIA data from Station ALOHA that was originally published in Hannides et al. (2020). They were originally published as the average d15N values of Ser, Gly, Phe, and Lys. The entire data set is now available here (Zenodo DOI forthcoming), and also as part of Miller et al. (in Revision; DOI forthcoming).
- AA-CSIA_5N.xlsx and AA-CSIA_8N.xlsx contain AA-CSIA data from the Equatorial Pacific that was originally published in Romero-romero et al. (2020), again as the average d15N values of Ser, Gly, Phe, and Lys. The entire data set is now available here (Zenodo DOI forthcoming).
- AA-CSIA_OSP.xlsx contains AA-CSIA data from Ocean Station Papa that was originally published in Shea et al. (2023).

Further model code applying the OMSM to real zooplankton data is on the way soon!
