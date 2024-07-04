studies <- c("N2309", "N2301")

for(study in studies){
  source(paste0(getwd(), "/Code/00_paths_pckgs_functions.R"))
  source(paste0(path_code, "10_import_convert.R"))
  source(paste0(path_code, "20_describe.R"))
  
  if(study=="N2301"){
    source(paste0(path_code, "40_select_mdl.R"))
  } else if(study=="N2309"){
    source(paste0(path_code, "45_mdl_dtls.R"))
    source(paste0(path_code, "50_evaluate_mdl.R"))
    source(paste0(path_code, "55_evaluate_mdl_counterfact.R"))
    
  }
  
  source(paste0(path_code, "900_report.R"))
}

