################################################################################
# ExternalValidation
################################################################################

# study <- "N2309"
# source(paste0(getwd(), "/Code/00_paths_pckgs_functions.R"))

load(paste0(path_data, "ana_dat_", ana_type,".RData"))
load(paste0(path_data, "dict_all_", ana_type, ".RData"))

load(paste0(path_main, "/N2301/Results/fnl_models_", ana_type,".RData"))

# reorder to the same order as before
fnls <- fnls[c(3, 4, 1, 2, 5, 6)]

n_cores <- detectCores(logical = T)

# to predict for everyday
targett <- 1:765

# to plot AUC/Brier monthly
targett_core <- seq(1:24)*30

cntrfctl_dat <- ana_dat

cntrfctl_dat$Drug <- with(cntrfctl_dat, 
                          factor(ifelse(Drug=="Placebo", 2, 1), 
                                 levels = c(1, 2), labels = c("Placebo", "FTY720")))

predict_val <- function(o,...) {
  
  source(paste0(path_code, "/00_paths_pckgs_functions.R"))
  
  # o <- fnls[[1]]
  
  attach(o)
  
  seed <- 14092022
  
  pattern <- outcome_cat$pattern[which(outcome_cat$col_labels==label)]
  
  test_dat_pre1 <- select(cl_dat, -c(1:2), -starts_with("outcome_"))
  
  test_dat <- data.frame(test_dat_pre1[,str_which(colnames(test_dat_pre1), paste0("_", pattern))], 
                         test_dat_pre1[, c((n_out*2+1):ncol(test_dat_pre1))]) %>%
    rename(c(censor = 1, day = 2)) 
  
  if(type == "forest"){
    
    forest_test <- develop.mdl(dat = test_dat, 
                               type = type, 
                               stump = T,
                               ntree = 1)
    
    predictions <- predict.all(object = model$model, 
                               dat = forest_test$data, 
                               targett = targett, 
                               type = type)
    
    
  } else if (type == "elastic"){
    
    #number of imputed datasets
    m <- 5
    
    elastic_test <- imp.forglm(dat = test_dat, seed = seed+12,
                               m = m)
    
    predictions_pre <- predict.all(object = model, 
                                   dat = elastic_test, 
                                   pred_type = "prob",
                                   targett = targett, 
                                   glm_lambda = glm_lambda, 
                                   type = type)
    
    predictions_pre$lp <- predict.all(object = model, 
                                   dat = elastic_test, 
                                   pred_type = "lp",
                                   glm_lambda = glm_lambda, 
                                   type = type)
    
    n_pred <- nrow(predictions_pre)/m
    
    predictions <- NULL
    
    for(p in 1:(n_pred)) {
      predictions <- rbind(predictions, 
                           apply(predictions_pre, 
                                 2, 
                                 function(x) {mean(x[c(p, p+((1:(m-1))*n_pred))])}))
    }
    
    predictions <- data.frame(predictions)
    
  } else if (type == "grelastic"){
    
    
    if(glm_lambda=="lambda.min"){
      # gr_lambda <- best_hyp[which(names(best_hyp)=="lambda.min")]
      gr_lambda <- model$model$lambda.min
    } else if(glm_lambda=="lambda.1se"){
      # gr_lambda <- best_hyp[which(names(best_hyp)=="lambda.1se")]
      gr_lambda <- model$model$lambda.1se
    }
    
    elastic_test <- imp.forglm(dat = test_dat, seed = seed+13, 
                               m = 1)
    
    predictions <- predict.all(object = model, 
                               dat = elastic_test, 
                               pred_type = "prob",
                               targett = targett, 
                               glm_lambda = gr_lambda,
                               type = type)
    
    predictions$lp <- predict.all(object = model, 
                                  dat = elastic_test, 
                                  pred_type = "lp",
                                  glm_lambda = gr_lambda,
                                  type = type)
    
    predictions <- data.frame(predictions)
    
  } else if (type == "tree"){
    
    tree_test <- develop.mdl(dat = test_dat, 
                             type = type, 
                             stump = T,
                             ntree = 1)
    
    predictions <- predict.all(object = model$model, 
                               dat = tree_test$data, 
                               targett = targett, 
                               type = type)
    
  }
  
  predictions_val <- predictions
  
  preds <- list(predictions_val = predictions, 
                type = type, 
                label = label, 
                val_dat = test_dat)
  
  detach(o)
  
  return(preds)
}

cl <- makeCluster(n_cores-1)
registerDoParallel(cl)

cl_dat <- ana_dat

clusterExport(cl, list("predict_val", "cl_dat", "outcome_cat", 
                       "path_code", "n_out", "study", "targett"))

tic("validation predict")

preds <- c(parLapply(cl, fnls, fun = predict_val))

toc()

stopCluster(cl)

save(preds, file = paste0(path_results, "preds_all_", ana_type,".RData"))

# load(file = paste0(path_results, "preds_all_", ana_type,".RData"))

# Describe predictions ####

# Counterfactual dataset #### 

cl <- makeCluster(n_cores-1)
registerDoParallel(cl)

cl_dat <- cntrfctl_dat

clusterExport(cl, list("predict_val", "cl_dat", "outcome_cat", 
                       "path_code", "n_out", "study", "targett"))

tic("counterfactual predict")

preds_cntrfctl <- c(parLapply(cl, fnls, fun = predict_val))

toc()

stopCluster(cl)

preds_allcntr <- preds

for(al in 1:length(preds)){
  
  preds_allcntr[[al]]$predictions_val$cntrfctl <- 
    preds_allcntr[[al]]$val_dat$cntrfctl <- 0
  
  preds_cntrfctl[[al]]$predictions_val$cntrfctl <- 
    preds_cntrfctl[[al]]$val_dat$cntrfctl <- 1
  
  preds_allcntr[[al]]$predictions_val$pat_id <- 
    preds_allcntr[[al]]$val_dat$pat_id <- 
    preds_cntrfctl[[al]]$predictions_val$pat_id <- 
    preds_cntrfctl[[al]]$val_dat$pat_id <- 1:nrow(ana_dat)
  
  preds_allcntr[[al]]$predictions_val <- rbind(preds_allcntr[[al]]$predictions_val, 
                                               preds_cntrfctl[[al]]$predictions_val)
  
  preds_allcntr[[al]]$val_dat <- rbind(preds_allcntr[[al]]$val_dat, 
                                       preds_cntrfctl[[al]]$val_dat)
  
}

save(preds_allcntr, file = paste0(path_results, "preds_allcntr_", ana_type,".RData"))

### Just to see what's going on
lapply(preds, function(x){summary(x$predictions_val$day_20)})

# View(preds[[1]]$predictions$day_765)
  
# Prob distribution ####
pred_dist <- NULL
for(i in preds){
  probs <- 1-i$predictions_val$day_720
  lbl <- outcome_cat$short_labels[which(outcome_cat$col_labels==i$label)]
  pred_dist <- rbind(pred_dist, 
                     data.frame(lbl = rep(lbl, length(probs)), probs = probs))
}

tab_probs <- NULL
for(i in 1:length(preds)){
  ids <- data.frame(pred = preds[[i]]$predictions_val$day_720, 
                    Drug = preds[[i]]$val_dat$Drug)
  
  tab_probs[[i]] <- data.frame(overall = c(min = min(ids$pred), 
                                              max = max(ids$pred), 
                                              med = median(ids$pred)), 
                                  fty720 = with(subset(ids, Drug=="FTY720"), 
                                                c(min = min(pred), 
                                                  max = max(pred), 
                                                  med = median(pred))),
                                  placebo = with(subset(ids, Drug=="Placebo"), 
                                                 c(min = min(pred), 
                                                   max = max(pred), 
                                                   med = median(pred))))
  names(tab_probs)[i] <- outcome_cat$short_labels[which(outcome_cat$col_labels==preds[[i]]$label)]
}

pred_dist$lbl <-   factor(pred_dist$lbl, levels=c("Relapse", "T2 MRI", "3m CDP", 
                                                  "Safety", "Immune safety", "Composite"))

plot_prob <- function(dat, scheme="color"){
  if(scheme=="color"){
    colcodes <- other_colcod
  } else if (scheme=="gray"){
    colcodes <- other_graycod
  }
  
  trgt <- ggplot(dat, aes(lbl, probs, fill=lbl)) +
    geom_boxplot() +
    scale_fill_manual(values =  colcodes[c(1,7,2,6,3,5)]) +
    labs(x = "Outcome", y = "Predicted risk") +
    coord_cartesian(ylim = c(-0.01, 1.01), xlim = c(0.4, 6.6), expand = F) + 
    theme(text = element_text(size = 7), 
          legend.position = "none", 
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle=45, vjust=1, hjust = 1))
  
  return(trgt)
}

fig_probs <- plot_prob(pred_dist)
fig_probs_gray <- plot_prob(pred_dist, scheme="gray")


calculate_auc_brier <- function(x){
  
  source(paste0(path_code, "/00_paths_pckgs_functions.R"))
  
  # x <- preds[[1]]
  
  attach(x)
  
  short_lab <- outcome_cat$short_labels[which(outcome_cat$col_labels==label)]
  
  val_dat_dca <- data.frame(censor = val_dat$censor, day = val_dat$day,
                            Model = (1-predictions_val$day_720))
  
  obj_dca <- dca(Surv(day, censor)~Model, data = val_dat_dca, time = 720) 
    # plot(test, show_ggplot_code = TRUE)
  
  plot_dca <- function(obj, label, scheme="color"){
    
    if(scheme=="color"){
      colcodes <- other_colcod
    } else if (scheme=="gray"){
      colcodes <- other_graycod
    }
  
  out_1 <- as_tibble(obj) %>%
    dplyr::filter(!is.na(net_benefit)) %>%
    ggplot(aes(x = threshold, y = net_benefit, color = label)) +
    geom_line(size = 1.2) +
    coord_cartesian(ylim = c(-0.2*obj_dca$prevalence, obj_dca$prevalence), expand = F) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) + 
    scale_color_manual(values =  colcodes[c(2, 1, 7)]) + 
    labs(x = "Threshold Probability", y = "Net Benefit", 
         title = paste(short_lab, "- decision curve")) +
    theme(text = element_text(size = 7), 
          legend.key = element_rect(fill = NA),
          legend.background = element_blank(), 
          legend.position = c(0.98, 0.98), 
          legend.justification = c(0.98, 0.98), 
          legend.direction = "vertical", 
          legend.box = "horizontal", 
          legend.margin = margin(0.1, 0.1, 0.1, 0.1, unit = "pt"), 
          legend.title = element_blank())
  
  out_2 <- as_tibble(net_intervention_avoided(obj)) %>%
    dplyr::filter(!is.na(net_intervention_avoided), !(variable %in% 
                                                        c("all", "none"))) %>%
    ggplot(aes(x = threshold, y = net_intervention_avoided, color = label)) +
    geom_line(size = 1.3) +
    coord_cartesian(ylim = c(-20, 100), expand = FALSE) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_color_manual(values =  colcodes[7]) +
    labs(x = "Threshold Probability", y = "Net reduction in interventions\nper 100 patients", 
         color = "",  title = paste(short_lab, "- intervention avoided")) +
    theme(text = element_text(size = 7), 
          legend.position = "none")
  
  out <- grid.arrange(out_1, out_2, ncol = 2)
  
  return(out)
  }
  
  fig_dca <- plot_dca(obj_dca, label = short_lab, scheme="color")
  fig_dca_gray <- plot_dca(obj_dca, label = short_lab, scheme="gray")
  
  # detach(x)
  # calibration (Crowson 2016)
  
  ### error due to absolute probability 1 for MRI up to 28 days
  
  expected <- NULL
  for(i in 1:nrow(predictions_val)){
    pre <- predictions_val[i, val_dat$day[i]]
    pre <- ifelse(pre==1, 1-0.0000000001, pre)
    expected <- c(expected, log(-log(pre)))
  }

  val_dat_sl <- data.frame(censor = val_dat$censor,  
                           expected = expected)
  # calibration in the large
  mdl_int <- glm(censor ~ offset(expected), 
             data = val_dat_sl, family = "poisson")
  cal_int <-summary(mdl_int)$coefficients[,1:2]
  
    
  mdl_sl <- glm(censor ~ expected, 
                data = val_dat_sl, family = "poisson")
  cal_sl <- summary(mdl_sl)$coefficients[2,1:2]
  cal <- rbind(cal_int, cal_sl)
  
  if("lp" %in% colnames(predictions_val)){
      ### this is not taking the baseline hazard into account (unlike poisson)
      val_dat_cox <- data.frame(censor = val_dat$censor, day = val_dat$day,
                               lp = predictions_val$lp)
      mdl_sl_cox <- coxph(Surv(day, censor)~lp, data = val_dat_cox)
      cal_sl_cox <- summary(mdl_sl_cox)$coefficients[1,c(1, 3)]
      cal <- rbind(cal, cal_sl_cox)
    }
    
    cal_intsl <- round(cbind(est = cal[,1], calculate_ci(cal[,1], cal[,2])), 2)
    
    cal_intsl <- rbind(cal_intsl, exp_cal_int = exp(cal_intsl[1,]), 
                       exp_cal_sl = exp(cal_intsl[2,]))

  preds_list <-  preds_list_visits <- preds_list_visits_fty <- preds_list_visits_pl <- NULL
  
  for(vis in targett_core){
    preds_list[[length(preds_list)+1]] <- list(1 - predictions_val[, vis])
  }
  
  for(vis in visits){
    preds_list_visits[[length(preds_list_visits)+1]] <- list(1 - predictions_val[, vis])
    
    preds_list_visits_fty[[length(preds_list_visits_fty)+1]] <- 
      list(1 - predictions_val[which(val_dat$Drug=="FTY720"), vis])
    
    preds_list_visits_pl[[length(preds_list_visits_pl)+1]] <- 
      list(1 - predictions_val[which(val_dat$Drug=="Placebo"), vis])
  }
  
  auc <- NULL
  brier <- NULL
  
  for(nv in 1:length(preds_list)){
    
    #The function expects increasing values for bad outcomes (i.e. risk)
    
    td_auc <- Score(preds_list[[nv]],
                    formula=Surv(day, censor)~1, 
                    metrics=c("auc", "Brier"),
                    data = val_dat, 
                    conf.int = TRUE, 
                    times = targett_core[nv])
    
    auc <- rbind(auc, data.frame(day = targett_core[nv], td_auc$AUC$score))
    
    brier <- rbind(brier, data.frame(day = targett_core[nv], td_auc$Brier$score))
  }
  
  levels(brier$model) <- c("null", "model")

  plot_aucbrier <- function(dat_auc, dat_brier, label, scheme="color"){
    
    if(scheme=="color"){
      colcodes <- other_colcod
    } else if (scheme=="gray"){
      colcodes <- other_graycod
    }
    
    fig_td_auc <- ggplot(dat_auc, aes(x = times, y = AUC)) +
      # geom_hline(yintercept = 0.5, aes(color = "black"), 
      #            size = 1, linetype = 2) +
      geom_ribbon(aes(ymin = lower, ymax = upper, 
                      fill = colcodes[6]), alpha = 0.4) +
      geom_point(aes(color = colcodes[6]), size = 2) + 
      geom_line(aes (color = colcodes[6]), size = 1.2) + 
      labs(x = "Time (days)", y = "Area under the curve (AUC)", 
           title = paste(label, "- AUC(t)")) + 
      scale_fill_manual(name = "95% CI", values =  colcodes[6]) + 
      scale_color_manual(name = "AUC(t)", values =  colcodes[6]) + 
      guides(color = guide_legend(label.theme = element_blank(), 
                                  order = 1), 
             fill = guide_legend(label.theme = element_blank(), 
                                 order = 2)) +
      coord_cartesian(ylim = c(0.5, 1), xlim = c(0, max(dat_auc$times)), expand = F) + 
      scale_x_continuous(breaks = c(0, targett_core[c(1:4)*6])) + 
      theme(text = element_text(size = 7), 
            legend.key = element_rect(fill = NA),
            legend.background = element_blank(), 
            legend.position = c(0.98, 0.98), 
            legend.justification = c(0.98, 0.98), 
            legend.direction = "vertical", 
            legend.box = "horizontal", 
            legend.margin = margin(0.1, 0.1, 0.1, 0.1, unit = "pt"))
    
    fig_td_brier <- ggplot(dat_brier, aes(x = times, y = Brier)) +
      geom_ribbon(aes(ymin = lower, ymax = upper, 
                      fill = model), alpha = 0.3) +
      geom_point(aes(color = model), size = 2, alpha = 0.8) + 
      geom_line(aes (color = model), size = 1.2, alpha = 0.8) + 
      labs(x = "Time (days)", y = "Brier score", 
           title = paste(label, "- Brier(t)")) + 
      scale_fill_manual(name = "95% CI", values =  colcodes[c(1,7)]) +
      scale_color_manual(name = "Brier(t)", values =  colcodes[c(1,7)]) + 
      guides(color = guide_legend(label.position = "left", order = 1), 
             fill = guide_legend(label = F, order = 2)) + 
      coord_cartesian(ylim = c(0, 0.25), xlim = c(0, max(dat_brier$times)), expand = F) + 
      scale_x_continuous(breaks = c(0, targett_core[c(1:4)*6])) + 
      theme(text = element_text(size = 7), 
            legend.key = element_rect(fill = NA),
            legend.background = element_blank(), 
            legend.position = c(0.98, 0.02), 
            legend.justification = c(0.98, 0.02), 
            legend.direction = "vertical", 
            legend.box = "horizontal", 
            legend.margin = margin(0.1, 0.1, 0.1, 0.1, unit = "pt"), 
            legend.title.align = 0.5)
    
    fig_td_aucbrier <- plot_grid(fig_td_auc, fig_td_brier, 
                                 ncol = 2)
    
    return(fig_td_aucbrier)
  }
  
  fig_aucbr <- plot_aucbrier(auc, brier, label = label)
  fig_aucbr_gray <- plot_aucbrier(auc, brier, label = label, scheme="gray")
  
  scores <- list(type = type,
                 label = label, 
                 tab_auc=auc, 
                 tab_brier = brier,
                 cal_intsl = cal_intsl, 
                 fig_dca = fig_dca, 
                 fig_dca_gray = fig_dca_gray, 
                 fig_aucbr = fig_aucbr, 
                 fig_aucbr_gray = fig_aucbr_gray, 
                 preds_list_visits = preds_list_visits,
                 preds_list_visits_fty =  preds_list_visits_fty,
                 preds_list_visits_pl = preds_list_visits_pl)
  
  detach(x)
  return(scores)
}

# Brier score & AUC ####

cl <- makeCluster(n_cores-1)
registerDoParallel(cl)


clusterExport(cl, list("targett_core", "path_code", "study", "calculate_auc_brier", "outcome_cat"))

# test <- NULL
# test[[1]] <- preds[[6]]

tic("model evaluate")

overall_auc_brier <- c(parLapply(cl, preds, fun = calculate_auc_brier))

toc()

stopCluster(cl)


# Calibration Plot ####
fig_cal_all <- list()

scheme_list <- c("color", "gray")

for(mdl in 1:length(overall_auc_brier)){
  
  attach(overall_auc_brier[[mdl]])
  
  short_lab <- outcome_cat$short_labels[which(outcome_cat$col_labels==label)]
  
  val_dat <- preds[[mdl]]$val_dat
  
  for(s in 1:length(scheme_list)){
    
    scheme <- scheme_list[s]
    
    if(scheme=="color"){
      colcodes <- drug_colcod
    } else if (scheme=="gray"){
      colcodes <- drug_graycod
    }
    
    
    plot.new()
    par(mfrow=c(3, 2))
    
    for(nv in 1:length(preds_list_visits)){
      
      td_plots <- Score(preds_list_visits[[nv]],
                        formula=Surv(day, censor)~1, 
                        data = val_dat, 
                        plots = c("Calibration", "ROC"),
                        conf.int = TRUE, 
                        times = visits[nv])
      
      td_plots_fty <- Score(preds_list_visits_fty[[nv]],
                            formula=Surv(day, censor)~1, 
                            data = subset(val_dat, Drug == "FTY720"),
                            plots = c("Calibration", "ROC"),
                            conf.int = TRUE, 
                            times = visits[nv])
      
      td_plots_pl <- Score(preds_list_visits_pl[[nv]],
                           formula=Surv(day, censor)~1, 
                           data = subset(val_dat, Drug == "Placebo"),
                           plots = c("Calibration", "ROC"),
                           conf.int = TRUE, 
                           times = visits[nv])

        plotCalibration(td_plots, legend = T, lty = 6, cens.method = "local", cex = 0.9, 
                        legend.cex = 0.9, legend.legend = "", lwd = 3,
                        legend.x = 0.35, legend.y = 0, legend.xjust = 1, legend.yjust = 0, 
                        plot.main = paste0(short_lab, " Calibration (month", visits[nv]/30, ")"),
                        method = "quantile")
        
        plotCalibration(td_plots_fty, col = colcodes[1], add = T, cens.method = "local", 
                        lwd = 2,
                        method = "quantile")
        plotCalibration(td_plots_pl, col = colcodes[2], add = T, cens.method = "local", 
                        lwd = 2,
                        method = "quantile")
      
      plotROC(td_plots, lty = 6, cex = 0.9, legend = F, lwd = 3, 
              plot.main = paste0(short_lab, " ROC (month", visits[nv]/30, ")"))
      
      plotROC(td_plots_fty, col = colcodes[1], add = T, lwd = 2, legend=F)
      plotROC(td_plots_pl, col = colcodes[2], add = T, lwd = 2, legend=F)
      
      legend(1, 0,  legend = c("Overall", "Drug=FTY720", "Drug=Placebo"), 
             col = c("#000000", colcodes), lty = c(3, 1, 1), lwd = c(3, 2, 2), cex = 0.9,
             xjust = 1, yjust = 0, bty = "n", adj = c(0, 0.5))
      
    }
  
    if(scheme=="color"){
     fig_cal_roc <- recordPlot()
    } else if (scheme=="gray"){
      fig_cal_roc_gray <- recordPlot()
      fig_cal_all[[mdl]] <- list(fig_cal_roc = fig_cal_roc, fig_cal_roc_gray=fig_cal_roc_gray)
    }
    
    par(mfrow=c(1,1))
  }
  
  detach(overall_auc_brier[[mdl]])
  
}

rm(fig_cal_roc, fig_cal_roc_gray)

# Make the outputs reportable ####

tab_auc_val <- NULL
tab_brier_val <- NULL
tab_cal_intsl <- NULL
fig_aucbr_val <- NULL
fig_aucbr_val_gray <- NULL
fig_dca_val <- NULL
fig_dca_val_gray <- NULL

for(ou in 1:n_out){
  
  attach(overall_auc_brier[[ou]])
  
  auc_pre <- subset(tab_auc[,c("day", "AUC", "lower", "upper")], 
                        day %in% visits)
  
  auc_pre2 <- cbind(outcome = rep(label, length(visits)), auc_pre)
  
  tab_auc_val <- rbind(tab_auc_val, auc_pre2)
  
  brier_pre <- subset(tab_brier[,c("model", "day", "Brier", "lower", "upper")], 
                        day %in% visits)
  
  brier_pre2 <- cbind(outcome = rep(label, length(visits)), brier_pre)
  
  tab_brier_val <- rbind(tab_brier_val, brier_pre2)
  
  tab_cal_intsl[[ou]] <- data.frame(cal_intsl)
  fig_aucbr_val[[ou]] <- fig_aucbr
  fig_aucbr_val_gray[[ou]] <- fig_aucbr_gray
  fig_dca_val[[ou]] <- fig_dca
  fig_dca_val_gray[[ou]] <- fig_dca_gray
  
  names(tab_cal_intsl)[ou] <- str_sub(label, end = 10)
  names(fig_aucbr_val)[ou] <- label
  names(fig_aucbr_val_gray)[ou] <- label
  
  detach(overall_auc_brier[[ou]])
}



# save the tables to be reported
save(list = ls()[str_which(ls(), "^tab_")], 
     file = paste0(path_results, "val_tab_", ana_type,".RData"))

# save the figures to be reported
save(list = ls()[str_which(ls(), "^fig_")], 
     file = paste0(path_results, "val_fig_", ana_type,".RData"))

