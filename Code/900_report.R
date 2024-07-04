################################################################################
# Prepare results for publication
################################################################################

# Development Results ####

# study <- "N2309"
# source(paste0(getwd(), "/Code/00_paths_pckgs_functions.R"))
source(paste0(getwd(), "/Code/09_functions_report.R"))
load(paste0(path_data, "dict_all_", ana_type, ".RData"))
trtsel_conv_table <- read.csv(paste0(getwd(), "/trtsel_conv_table.csv"))

# DESCRIPTION TABLES ####
load(paste0(path_results, "desc_base_", ana_type,".RData"))

export(list(miss_pred = desc_miss_pred, miss_pat  = desc_miss_pat), 
       paste0(path_results, "Export/", ana_type, "/", study, "_missdesc_", ana_type,".xlsx"))

export(desc_predcat, 
       paste0(path_results, "Export/", ana_type, "/", study, "_categories_", ana_type,".xlsx"))

export(make_num_great(desc_bi), 
       paste0(path_results, "Export/", ana_type, "/", study, "_binaryout_", ana_type,".xlsx"))

export(make_desc_exp(desc_by_drug),
       paste0(path_results, "Export/", ana_type, "/", study, "_basedesc_bydrug_", 
              ana_type,".xlsx"))

export(make_desc_exp(desc_by_miss),
       paste0(path_results, "Export/", ana_type, "/", study, "_basedesc_bymiss_", 
              ana_type,".xlsx"))

export(make_desc_exp(desc_by_censored), 
       paste0(path_results, "Export/", ana_type, "/", study, "_basedesc_bycens_", 
              ana_type,".xlsx"))

export(make_desc_exp(desc_by_outcome),
       paste0(path_results, "Export/", ana_type, "/", study, "_basedesc_byout_", 
              ana_type,".xlsx"))

export(make_desc_exp(desc_predic_all),
       paste0(path_results, "Export/", ana_type, "/", study, "_basedesc_all_", 
              ana_type,".xlsx"))

export(make_epv(desc_by_outcome), 
       paste0(path_results, "Export/", ana_type, "/", study, "_epv_byout_",
              ana_type, ".xlsx"))

# LOOP TABLES ####
if(study=="N2301"){
  load(paste0(path_results, "loops_tab_", ana_type, ".RData"))
  
  export(lapply(tab_auc_all, round, digits = 3), 
         paste0(path_results, "Export/", ana_type, "/", study, "_dev_aucs_",
                ana_type, ".xlsx"), rowNames=T)
  
  export(lapply(tab_brier_all, round, digits = 3), 
         paste0(path_results, "Export/", ana_type, "/", study, "_dev_briers_",
                ana_type, ".xlsx"), rowNames=T)
  
  export(tab_npred, 
         paste0(path_results, "Export/", ana_type, "/", study, "_dev_npreds_",
                ana_type, ".xlsx"))
  
  export(make_varimp(tab_imp_var),
         paste0(path_results, "Export/", ana_type, "/", study, "_dev_varimp_",
                ana_type, ".xlsx"))
} else {
  # VALIDATION TABLES ####
  
  load(paste0(path_results, "mdls_tab_", ana_type, ".RData"))
  load(paste0(path_results, "val_tab_", ana_type, ".RData"))
  load(paste0(path_results, "ctr_tab_", ana_type,".RData"))
  
  export(make_aucbr(tab_auc_val), 
         paste0(path_results, "Export/", ana_type, "/", study, "_val_auc_",
                ana_type, ".xlsx"))
  
  export(tab_probs, 
         paste0(path_results, "Export/", ana_type, "/", study, "_probs_",
                ana_type, ".xlsx"))
  
  export(tab_basehaz, 
         paste0(path_results, "Export/", ana_type, "/", study, "_mdl_bazehaz_",
                ana_type, ".xlsx"))
  
  export(make_aucbr(tab_brier_val), 
         paste0(path_results, "Export/", ana_type, "/", study, "_val_brier_",
                ana_type, ".xlsx"))
  
  export(make_cal_intsl(tab_cal_intsl), 
         paste0(path_results, "Export/", ana_type, "/", study, "_val_cal_",
                ana_type, ".xlsx"))
  
  export(tab_fnls, 
         paste0(path_results, "Export/", ana_type, "/", study, "_fnl_mdls_",
                ana_type, ".xlsx"))
  
  export(lapply(tab_coef, make_coef_great), 
         paste0(path_results, "Export/", ana_type, "/", study, "_mdl_coef_",
                ana_type, ".xlsx"))
  
  export(make_cal_intsl(tab_evalsmry, trtsel.eval = T), 
         paste0(path_results, "Export/", ana_type, "/", study, "_trteval_",
                ana_type, ".xlsx"))
  
  export(lapply(tab_trteffsmry, round, 3), 
         paste0(path_results, "Export/", ana_type, "/", study, "_trtsmry_",
                ana_type, ".xlsx"))

}

# DESCRIPTION FIGURES ####
load(paste0(path_results, "desc_base_fig_", ana_type,".RData"))

if(study=="N2301"){
  load(paste0(path_results, "loops_fig_", ana_type,".RData"))
} else{
  load(paste0(path_results, "val_fig_", ana_type,".RData"))
  load(paste0(path_results, "ctr_fig_", ana_type,".RData"))
}


# bi_outcome <- grid.arrange(fig_desc_bi_all, fig_desc_bi_fnl, nrow = 2)

schemes <- c("", "_gray")

plot_all <- function(ext, devc){
  for(s in schemes){
    surv_outcome <- arrange_ggsurvplots(get(paste0("fig_surv", s)), 
                                        ncol = 2, 
                                        nrow = 3)
    
    ggsave(filename = paste0(study, "_figdesc_surv_", ana_type, s,ext),
           plot = surv_outcome,
           device = devc, 
           path = paste0(path_results, "Export/", ana_type),
           width = 180,
           height = 210,
           units = "mm", 
           dpi = 600)
    
    ggsave(filename = paste0(study, "_figdesc_bi_all_", ana_type, s,ext),
           plot = get(paste0("fig_desc_bi_all", s)),
           device = devc, 
           path = paste0(path_results, "Export/", ana_type),
           width = 180,
           height = 75,
           units = "mm", 
           dpi = 600)
    
    ggsave(filename = paste0(study, "_figdesc_bi_fnl_", ana_type, s,ext),
           plot = get(paste0("fig_desc_bi_fnl", s)),
           device = devc, 
           path = paste0(path_results, "Export/", ana_type),
           width = 180,
           height = 75,
           units = "mm", 
           dpi = 600)
    
    ggsave(filename = paste0(study, "_fig_miss_", ana_type, s,ext),
           plot = get(paste0("fig_miss", s)),
           device = devc, 
           path = paste0(path_results, "Export/", ana_type),
           width = 180,
           height = 225,
           units = "mm", 
           dpi = 600)
    
    # LOOP FIGURES ####
    
    # plot.new()
    
    if(study=="N2301"){
    
    trees <- plot_grid(plotlist = get(paste0("fig_tree", s)), ncol = 2, byrow = F)
    
    ggsave(filename = paste0(study, "_figtrees_", ana_type, s,ext),
           plot = trees,
           device = devc, 
           path = paste0(path_results, "Export/", ana_type),
           width = 180,
           height = 210,
           units = "mm", 
           dpi = 600)
    
    impvar <- plot_grid(plotlist = get(paste0("fig_var_imp", s)), ncol = 2, byrow = F)
    
    ggsave(filename = paste0(study, "_figvarimp_", ana_type, s,ext),
           plot = impvar,
           device = devc,
           path = paste0(path_results, "Export/", ana_type),
           width = 180,
           height = 210,
           units = "mm", 
           dpi = 600)
    } else{
      
      for(i in 1:length(fig_cal_all)){
        attach(fig_cal_all[[i]])
        
        plot.new()
        
        if(devc=="pdf"){
          pdf(file = paste0(path_results, "Export/", ana_type, "/", study, 
                            "_figcal",i , "_", ana_type, s, ext), 
              width = 7.087,
              height = 7.874)
          
        } else if(devc=="tiff"){
          tiff(file = paste0(path_results, "Export/", ana_type, "/", study, 
                             "_figcal",i , "_", ana_type, s, ext), 
              width = 180,
              height = 210,
              units = "mm", 
              res = 600)
        }
        
        replayPlot(get(paste0("fig_cal_roc", s)))
        
        dev.off()
        
        detach(fig_cal_all[[i]])
      }
      
      aucbr_1 <- plot_grid(plotlist = get(paste0("fig_aucbr_val", s))[c(1:3)], nrow = 3)
      aucbr_2 <- plot_grid(plotlist = get(paste0("fig_aucbr_val", s))[c(4:6)], nrow = 3)
      
      ggsave(filename = paste0(study, "_aucbr1_", ana_type, s,ext),
             plot = aucbr_1,
             device = devc,
             path = paste0(path_results, "Export/", ana_type),
             width = 180,
             height = 210,
             units = "mm", 
             dpi = 600)
      
      ggsave(filename = paste0(study, "_aucbr2_", ana_type, s,ext),
             plot = aucbr_2,
             device = devc,
             path = paste0(path_results, "Export/", ana_type),
             width = 180,
             height = 210,
             units = "mm", 
             dpi = 600)
      
      dca_1 <- plot_grid(plotlist = get(paste0("fig_dca_val", s))[c(1:3)], nrow = 3)
      dca_2 <- plot_grid(plotlist = get(paste0("fig_dca_val", s))[c(4:6)], nrow = 3)
      
      ggsave(filename = paste0(study, "_dca1_", ana_type, s,ext),
             plot = dca_1,
             device = devc,
             path = paste0(path_results, "Export/", ana_type),
             width = 180,
             height = 210,
             units = "mm", 
             dpi = 600)
      
      ggsave(filename = paste0(study, "_dca2_", ana_type, s,ext),
             plot = dca_2,
             device = devc,
             path = paste0(path_results, "Export/", ana_type),
             width = 180,
             height = 210,
             units = "mm", 
             dpi = 600)
      
      trtsel_ls <- NULL
      trtsel_ls_extra <- NULL
      my <- 1
      my_extra <- 1
      
      for(i in 1:length(get(paste0("fig_trtsel", s)))){
        pertype <- get(paste0("fig_trtsel", s))[[i]]
        fig_ls <- NULL
        for(f in 1:length(pertype)){
          fig_ls[[f]] <- pertype[[f]]$plot
        }
        trtsel_ls[[my]] <- fig_ls[[1]]
        trtsel_ls[[my+1]] <- fig_ls[[6]]
        my <- length(trtsel_ls)+1
        
        trtsel_ls_extra[[my_extra]] <- fig_ls[[2]]
        trtsel_ls_extra[[my_extra+1]] <- fig_ls[[3]]
        trtsel_ls_extra[[my_extra+2]] <- fig_ls[[4]]
        my_extra <- length(trtsel_ls_extra)+1
      }
      
      trtsel_grd_1 <- plot_grid(plotlist = trtsel_ls[1:6], ncol = 2)
      trtsel_grd_2 <- plot_grid(plotlist = trtsel_ls[7:12], ncol = 2)
      
      ggsave(filename = paste0(study, "_trtsel_1_", ana_type, s,ext),
             plot = trtsel_grd_1,
             device = devc,
             path = paste0(path_results, "Export/", ana_type),
             width = 180,
             height = 210,
             units = "mm", 
             dpi = 600)
      
      ggsave(filename = paste0(study, "_trtsel_2_", ana_type, s,ext),
             plot = trtsel_grd_2,
             device = devc,
             path = paste0(path_results, "Export/", ana_type),
             width = 180,
             height = 210,
             units = "mm", 
             dpi = 600)
      
      trtsel_grd_extra_1 <- plot_grid(plotlist = trtsel_ls_extra[1:9], 
                                    ncol = 3)
      
      trtsel_grd_extra_2 <- plot_grid(plotlist = trtsel_ls_extra[10:18], 
                                      ncol = 3)
      
      ggsave(filename = paste0(study, "_trtsel_extra_1_", ana_type, s,ext),
             plot = trtsel_grd_extra_1,
             device = devc,
             path = paste0(path_results, "Export/", ana_type),
             width = 180,
             height = 210,
             units = "mm", 
             dpi = 600)
      
      ggsave(filename = paste0(study, "_trtsel_extra_2_", ana_type, s,ext),
             plot = trtsel_grd_extra_2,
             device = devc,
             path = paste0(path_results, "Export/", ana_type),
             width = 180,
             height = 210,
             units = "mm", 
             dpi = 600)
      
      ggsave(filename = paste0(study, "_probs_", ana_type, s,ext),
             plot = get(paste0("fig_probs", s)),
             device = devc, 
             path = paste0(path_results, "Export/", ana_type),
             width = 90,
             height = 75,
             units = "mm", 
             dpi = 600)
      
      
    }
  }
}

plot_all(ext = ".pdf", devc = "pdf")
plot_all(ext = ".tif", devc = "tiff")

