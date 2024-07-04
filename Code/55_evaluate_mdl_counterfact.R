################################################################################
# ExternalValidation
################################################################################

# ATTENTION: This script requires debugged TreatmentSelection package according to "_Trtsel_debug.R"

# study <- "N2309"
# source(paste0(getwd(), "/Code/00_paths_pckgs_functions.R"))

load(paste0(path_data, "dict_all_", ana_type, ".RData"))

load(paste0(path_results, "preds_allcntr_", ana_type,".RData"))

targett <- visits[3]

tab_trteffsmry <- NULL
tab_evalsmry <- NULL
fig_trtsel <- NULL
fig_trtsel_gray <- NULL

tic("trtsel evaluate")

for(i in 1:length(preds_allcntr)){
  
  x <- preds_allcntr[[i]]
  
  attach(x)
  
  short_lab <- outcome_cat$short_labels[which(outcome_cat$col_labels==label)]
  
  dat_trtsel <- select(val_dat, pat_id, day, censor, Drug, cntrfctl) %>% 
    inner_join(select(predictions_val, pat_id, cntrfctl, day_720)) %>%
    select(-cntrfctl) %>%
    pivot_wider(names_from = Drug, values_from = day_720) %>%
    inner_join(select(filter(val_dat, cntrfctl==0), pat_id, Drug)) %>%
    mutate(Drug=ifelse(Drug=="FTY720", 1, 0),
           FTY720 = 1-FTY720, 
           Placebo = 1-Placebo) %>%
    data.frame()
  
  # The outcome of below is like a decision curve analysis
  
  obj_trtsel <- trtsel(formula = Surv(day, censor) ~ Drug, data=dat_trtsel, treatment.name = "Drug", 
                       prediction.time = targett, thresh = 0, 
                       fittedrisk.t0 = "Placebo", fittedrisk.t1 = "FTY720")
  
  plot_trtsel <- function(obj, label, scheme = "color"){
    
    if(scheme=="color"){
      colcodes <- drug_colcod
    } else if (scheme=="gray"){
      colcodes <- drug_graycod
    }
    
    fig_trtsel_risk <- plot(obj, plot.type = "risk", ci = "vertical",
                            trt.names = c("FTY720", "Placebo"), 
                            main = paste0(label, " - risk curves"), 
                            clr = colcodes)
    
    fig_trtsel_cal <- calibrate(obj, plot.type="calibration",
                                point.color = colcodes, 
                                ylim = c(0, 1),
                                main = paste0(label, " - calibration plot"),
                                trt.names = c("FTY720", "Placebo"),
                                xlim = c(0,1), groups = 5) #, main = paste0(label, " - impact curve"))
    fig_trtsel_caltr <- calibrate(obj, plot.type="risk.t1", 
                                  line.color = colcodes[1], point.color = colcodes[1], 
                                  main = paste0(label, " - risk curve FTY720"),
                                  ylim = c(0, 1), xlim = c(0, 100), groups = 5)
    fig_trtsel_calpl <- calibrate(obj, plot.type="risk.t0", 
                                  line.color = colcodes[2], point.color = colcodes[2], 
                                  main = paste0(label, " - risk curve Placebo"),
                                  ylim = c(0, 1), xlim = c(0, 100), groups = 5)
    
    # The below option is just like Treatment Effect Distribution
    # plot(obj_trtsel, plot.type = "cdf", 
    #      main = paste0(label, " - treatment effect distribution"))
    
    if(scheme=="color"){
      colcodes <- other_colcod
    } else if (scheme=="gray"){
      colcodes <- other_graycod[c(7:1)]
    }
    
    # what exactly does this display?
    fig_trtsel_selimp <- plot(obj, plot.type = "selection impact", 
                              ylim = c(0, 1), main = paste0(label, " - impact curve"), clr = colcodes[1])
    
    fig_trtsel_trteff <- plot(obj, plot.type = "treatment effect", 
                              ylim = c(min(min(obj$derived.data$trt.effect), -0.05), 
                                       max(obj$derived.data$trt.effect)), 
                              main = paste0(label, " - treatment effect distribution"), clr = colcodes[1])
    
    return(list(fig_trtsel_risk = fig_trtsel_risk, 
                fig_trtsel_cal = fig_trtsel_cal, 
                fig_trtsel_caltr = fig_trtsel_caltr, 
                fig_trtsel_calpl = fig_trtsel_calpl, 
                fig_trtsel_selimp = fig_trtsel_selimp,
                fig_trtsel_trteff = fig_trtsel_trteff))
    
  }
  
  fig_trtsel[[i]] <- plot_trtsel(obj_trtsel,short_lab)
  fig_trtsel_gray[[i]] <- plot_trtsel(obj_trtsel,short_lab, scheme = "gray")
  
  # contintue from here
  eval_obj <- evaluate(obj_trtsel)
  
  tab_evalsmry[[i]] <- data.frame(t(rbind(estimate = data.frame(eval_obj$estimates[,-length(eval_obj$estimates)]), 
                                  data.frame(eval_obj$conf.intervals))))
  
  tab_trteffsmry[[i]] <- summary(obj_trtsel$derived.data$trt.effect)
  
  names(tab_evalsmry)[i]<- names(tab_trteffsmry)[i] <- short_lab

  detach(x)
}

toc()

# save the tables to be reported
save(list = ls()[str_which(ls(), "^tab_")], 
     file = paste0(path_results, "ctr_tab_", ana_type,".RData"))

# save the figures to be reported
save(list = ls()[str_which(ls(), "^fig_")], 
     file = paste0(path_results, "ctr_fig_", ana_type,".RData"))


