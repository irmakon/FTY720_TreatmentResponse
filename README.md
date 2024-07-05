# FTY720_TreatmentResponse
## Aim
This README file explains how the code in this repository transforms and analyzes the data retrieved via [Clinical Study Data Request](https://www.clinicalstudydatarequest.com) ([Proposal 11223](https://www.clinicalstudydatarequest.com/Posting.aspx?ID=20600&GroupID=SUMMARIES)) from the sponsor (Novartis) of the [FREEDOMS](https://clinicaltrials.gov/study/NCT00289978) and [FREEDOMS II](https://clinicaltrials.gov/study/NCT00355134) randomized clinical trials (RCT) to generate the [PhD thesis](https://edoc.ub.uni-muenchen.de/32918/) of Begum Irmak On submitted to [LMU Munich](https://www.lmu.de).

## Folder [DataManagement](https://github.com/irmakon/FTY720_TreatmentResponse/tree/main/DataManagement)
This folder contains SAS codes to transform the raw data tables, as shared by the sponsor, into datasets required for the analysis. In order to define the outcomes and choose/define the predictors, the data is described in the files named “basicdesc”.

Then, the selected predictors (“predictors”) and outcomes (“outcome”) are defined for the intention to treat (ITT) population, respectively for the FREEDOMS (“N2301”) and FREEDOMS II (“N2309”) trials.

## Folder [Analysis](https://github.com/irmakon/FTY720_TreatmentResponse/tree/main/Analysis)
This folder contains R codes to describe and analyze the predictors and outcomes tables to meet the objectives of the PhD thesis, namely predicting prognosis and response to fingolimod in relapsing-remitting multiple sclerosis patients. 

To run the analysis, one only needs to run the *000_master.R*, which sources all the other files. The session info can be found in *SessionInfo.txt*. The following .R files are sourced by the master, respectively:

- *00_paths_pckgs_functions* defines all the paths, calls packages, and defines most of the user-defined functions
- *01_studyparam* defines the meta-parameters used in the analysis
- *09_functions_report* defines user-defined functions for reporting
- *10_import_convert* imports the managed datasets and finalizes any data preprocessing, creating the ready-to-analyze datasets
- *20_describe* describes the baseline predictors, outcome and missing data
- *40_select_mdl* imputes the missing data, tunes hyperparameters, runs the benchmark study to select the modeling method and fits the best one to generate the final model for each outcome in the FREEDOMS dataset, and reports internal validation measures alongside the important variables
- *45_mdl_dtls* reports the details of the finally selected models
- *50_evaluate_mdl_overall* externally validates the prognostic aspects of the final models in the FREEDOMS II dataset
- *55_evaluate_mdl_counterfact* creates the counterfactual dataset from FREEDOMS II and evaluates the treatment response predicted by the final models.
- *900_report* reports the findings as Excel tables or colored and black and white figures
- *_Trtsel_debug* corrects a debug identified in functions in the Trtsel package.


