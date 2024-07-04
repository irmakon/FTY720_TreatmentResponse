################################################################################
# Reporting the final models
################################################################################

# study <- "N2309"
# source(paste0(getwd(), "/Code/00_paths_pckgs_functions.R"))

load(paste0(path_data, "ana_dat_", ana_type,".RData"))
load(paste0(path_data, "dict_all_", ana_type, ".RData"))

load(paste0(path_main, "/N2301/Results/fnl_models_", ana_type,".RData"))

fnls <- fnls[c(3, 4, 1, 2, 5, 6)]

fnls_tbl <- NULL

for(i in 1:length(fnls)){
  best_hyp <- fnls[[i]]$best_hyp
  hyp <- NULL
  if(fnls[[i]]$type %in% c("tree", "forest")){
    for(h in 1:ncol(best_hyp)){
      hyp <- paste0(hyp, colnames(best_hyp)[h], "=", best_hyp[1,h], ", ")
    }
  } else{
    for(h in 1:length(best_hyp)){
      hyp <- paste0(hyp, names(best_hyp)[h], "=", 
                    round(best_hyp[h], 2), ", ")
    }
  }
  
  fnls_tbl <- rbind(fnls_tbl, c(fnls[[i]]$label, fnls[[i]]$type, hyp))
}

tab_fnls <- data.frame(fnls_tbl)
colnames(tab_fnls) <- c("outcome", "method", "hyperparameters")

tab_coef <- NULL
tab_basehaz <- NULL

for(i in 1:length(fnls)){
  if(!(fnls[[i]]$type %in% c("elastic", "grelastic"))){
    tab_coef[[i]] <- tab_basehaz[[i]] <- data.frame(fnls[[i]]$label)
  } else {
    attach(fnls[[i]])
    
    dat_base <- model$data[1,-c(1:3)]
    dat_base <- rbind(dat_base, 0)[2,]
    
    if(fnls[[i]]$type=="elastic"){
      coef_list_pre <- as.matrix(predict(model$model, s = glm_lambda, type = "coefficients"))
      coef_list_pre2 <- data.frame(variable = rownames(coef_list_pre), 
                                   coefficient = as.vector(coef_list_pre))
      
      base_fit <- survfit(model$model, s = glm_lambda, x=as.matrix(model$data[,-c(1:3)]), 
                          y = model$data$outcome, weights = model$data$fit_weight, 
                          newx = as.matrix(dat_base))
      
      test <- base_fit$cumhaz
      
      pre_basehaz <- with(base_fit, (data.frame(time = time, cumhaz = cumhaz)))
      
      pre2_basehaz <- data.frame(time =  1:visits[3]) %>% 
        left_join(pre_basehaz)
      
      for(c in 1:nrow(pre2_basehaz)){
        if(is.na(pre2_basehaz$cumhaz[c]))
          pre2_basehaz$cumhaz[c] <- pre2_basehaz$cumhaz[c-1]
      }
      
      tab_basehaz[[i]] <- pre2_basehaz
      
    } else if(fnls[[i]]$type=="grelastic"){
      coef_list_pre <- coef(model$model, lambda = model$model$lambda.min)
      coef_list_pre2 <- data.frame(variable = names(coef_list_pre), 
                                   coefficient = as.vector(coef_list_pre))
      
      time <- 1:visits[3]
      
      base_fit <- predict(model$model, as.matrix(dat_base), 
                          type="hazard", lambda = model$model$lambda.min)(time)
      
      tab_basehaz[[i]] <- data.frame(time = time, cumhaz = base_fit)
    }
    
    tab_coef[[i]] <- subset(coef_list_pre2, coefficient!=0)
    
    detach(fnls[[i]])
  }
  names(tab_coef)[i] <- names(tab_basehaz)[i] <- str_replace(str_sub(fnls[[i]]$label, end = 10), 
                                                             "/", "")
}

# save the tables to be reported
save(list = ls()[str_which(ls(), "^tab_")], 
     file = paste0(path_results, "mdls_tab_", ana_type,".RData"))

### Check if the basehazard and coefficients give the predictions!

# # load(paste0(path_results, "mdls_tab_", ana_type,".RData"))
# # load(file = paste0(path_results, "preds_all_", ana_type,".RData"))
# 
# # tested for m=c(1, 2, 5, 6) 
# 
# m <- 6
# 
# attach(preds[[m]])
# 
# 
# vars <- str_replace(str_replace(tab_coef[[m]]$variable, "\\.{3}", " - "), "\\.{1}", ":")
# 
# #### Only check for non-missing
# # check lp
# missrow <- which(apply(val_dat, 1, function(x){any(is.na(x))})==T)
# 
# mdl_mtrx <- model.matrix(lm(censor~.*Drug, data=val_dat))
# 
# dat_mtrx <- mdl_mtrx[,vars]
# lp_cal <- dat_mtrx %*% tab_coef[[m]]$coefficient
# 
# all(predictions_val$lp[-missrow]==lp_cal)
# # TRUE
# 
# # check predictions
# surv_cal <- t(exp(-(tab_basehaz[[m]]$cumhaz %*% t(exp(lp_cal)))))
# 
# surv_slct <- predictions_val[-missrow, tab_basehaz[[m]]$time]
# 
# all(round(surv_cal, 10)==round(surv_slct, 10))
# # TRUE for all except FALSE for m=6
# 
# # for m=6, the below is acceptable:
# diff <- surv_cal-surv_slct
# min(diff)
# # 2.482228e-06
# max(diff)
# # 0.0003708107
# 
# detach(preds[[m]])
