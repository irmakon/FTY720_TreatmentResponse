################################################################################
# Fit 1) transformation random forest 2) transformation tree
# 3) elasticnet 4) grouped lasso with ridge penalization
# Internally evaluate them
# All loops all outcomes lasts around 35 hours after parallelization
################################################################################

# study <- "N2301"
# source(paste0(getwd(), "/Code/00_paths_pckgs_functions.R"))

load(paste0(path_data, "ana_dat_", ana_type,".RData"))
load(paste0(path_data, "dict_all_", ana_type, ".RData"))

n_cores <- detectCores(logical = T)

### define the combinations of methods & outcomes
nout <- 1:n_out

# nout <- 1
type <- c( "tree", "forest", "elastic", "grelastic") #  "tree", "forest", "elastic", "grelastic"

loop_dat_pre <- as.matrix(expand.grid(nout = nout, type = type))

loop_dat <- NULL
for(l in 1:(n_out*length(type))){
  loop_dat[[l]] <- loop_dat_pre[l, ]
}

# o <- loop_dat[[5]]

### the function for model choice & development
loop.one <- function(o, ...) {

  source(paste0(path_code, "/00_paths_pckgs_functions.R"))

  # get the necessary parameters for type of outcome
  nout <- as.integer(o[1])
  type <- as.character(o[2])
  
  # target time for prediction and evaluation
  targett <- visits
  
  # create the model dataset based on the outcome
  mdl_dat_pre <- select(ana_dat, -c(1:2), -starts_with("outcome_"))
  
  outcome_nam <- str_split_fixed(colnames(mdl_dat_pre)[2*nout-1], n =2, pattern = "_")[,2]
  
  mdl_dat <- data.frame(mdl_dat_pre[, c((2*nout-1):(2*nout))], 
                        mdl_dat_pre[, c((n_out*2+1):ncol(mdl_dat_pre))]) %>%
    rename(c(censor = 1, day = 2))
  
  lab <- outcome_cat$col_labels[which(outcome_cat$pattern==outcome_nam)]
  
  # create the folds
  set.seed(230622)
  cv_sets <- createFolds(y = mdl_dat$Drug, k = n_folds)
  
  # initialize stuff
  seed <- 2022
  
  # f <- 1
  imp_var_folds <- NULL
  for_avg_varimp <- NULL
  auc <- NULL
  brier <- NULL
  allpredictions <- NULL
  
  # for all the CV-folds + a final fold where the whole dataset is used
  for(f in (1:(n_folds+1))){
    ### n_folds+1 is for the final modeling
      if(f==(n_folds+1)){
        train_dat <- mdl_dat
        final <- T
      } else{
        #create the train and test datasets
        train_dat <- mdl_dat[-cv_sets[[f]],]
        test_dat <- mdl_dat[cv_sets[[f]],]
        final <- F
      }

    if(type == "forest"){
      forest_cox <- develop.mdl(dat = train_dat, 
                                seed = seed+f, 
                                type = type, 
                                final = final) #
      
      if(final==T){
        fnl_mdl <- forest_cox
        break
      } 
      forest_test <- develop.mdl(dat = test_dat, 
                                 type = type, 
                                 stump = T, 
                                 ntree = 1)
      
      predictions <- predict.all(object = forest_cox$model, 
                                 dat = forest_test$data, 
                                 targett = targett, 
                                 type = type)
      
      preds_list <- NULL
      for(vis in 1:ncol(predictions)){
        preds_list[[vis]] <- list(1 - predictions[, vis])
      }
      
      #calculate the variable importances at each loop
      var_imp <- vimp.cv(forest_cox$model, forest_test$data, nperm = 3)
      imp_var_folds <- c(imp_var_folds, names(var_imp[which(var_imp>abs(min(var_imp)))]))
      for_avg_varimp <- cbind(for_avg_varimp, var_imp)
    } else if (type == "elastic"){
      #number of imputed datasets
      m <- 5
      
      list_elastic <- develop.mdl(dat = train_dat, 
                                  seed = seed+f, 
                                  type = type,  
                                  m = m, 
                                  final = final)
      if(final==T){
        fnl_mdl <- list_elastic
        break
      }
      
      elastic_test <- imp.forglm(dat = test_dat, seed = seed+f+1,
                                 m = m)
      
      coef_list <- coef(list_elastic$model, s = glm_lambda)
      
      predictions_pre <- predict.all(object = list_elastic, 
                                 dat = elastic_test, 
                                 pred_type = "prob",
                                 targett = targett, 
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
      
      sel_var <- attr(coef_list, "Dimnames")[[1]][attr(coef_list, "i")+1]
      imp_var_folds <- c(imp_var_folds, sel_var)
      
      preds_list <- NULL
      for(vis in 1:ncol(predictions)){
        preds_list[[vis]] <-  list(1 - predictions[, vis])
      }
    } else if (type == "grelastic"){
      list_grelastic <- develop.mdl(dat = train_dat, 
                                    seed = seed+f, 
                                    type = type,
                                    m = 1, 
                                    final = final)
      
      if(final==T){
        fnl_mdl <- list_grelastic
        break
      }
      
      elastic_test <- imp.forglm(dat = test_dat, seed = seed+f+1, 
                                 m = 1)
      
      if(glm_lambda=="lambda.min"){
        gr_lambda <- list_grelastic$model$lambda.min
      } else if(glm_lambda=="lambda.1se"){
        gr_lambda <- list_grelastic$model$lambda.1se
      }
      
      predictions <- predict.all(object = list_grelastic, 
                                 dat = elastic_test, 
                                 pred_type = "prob",
                                 targett = visits, 
                                 glm_lambda = gr_lambda,
                                 type = type)
      
      predictions <- data.frame(predictions)
      
      preds_list <- NULL
      
      for(vis in 1:ncol(predictions)){
        preds_list[[vis]] <-  list(1 - predictions[, vis])
      }
      
      coef_list <- coef(list_grelastic$model, lambda = gr_lambda)
      
      sel_var <- attr(coef_list, "names")[coef_list!=0]
      
      if(length(sel_var)==0){
        gr_lambda <- list_grelastic$model$lambda.min
        coef_list <- coef(list_grelastic$model, lambda = gr_lambda)
        
        sel_var <- attr(coef_list, "names")[coef_list!=0]
      }
      
      imp_var_folds <- c(imp_var_folds, sel_var)
      
    } else if (type == "tree"){
      
      tree_cox <- develop.mdl(dat = train_dat, seed = seed+f, 
                              type = type, 
                              final = final)
      if(final==T){
        fnl_mdl <- tree_cox
        break
      } 
      
      tree_test <- develop.mdl(dat = test_dat, 
                               type = type, 
                               stump = T, 
                               ntree = 1)
      
      predictions <- predict.all(object = tree_cox$model, 
                                 dat = tree_test$data, 
                                 targett = targett, 
                                 type = type)
      
      preds_list <- NULL
      for(vis in 1:ncol(predictions)){
        preds_list[[vis]] <- list(1 - predictions[, vis])
      }
      
      node_var <- NULL
      if(length(tree_cox)>1){
        for(n in 1:((length(tree_cox)-1)/2)){
          node_var <- c(node_var, names(info_node(node_party(tree_cox[[n]]))$p.value))
        }
        
        sel_var <- unique(node_var)
      } else{
        sel_var <- NA # "none"
      }
      
      imp_var_folds <- c(imp_var_folds, sel_var)
    }
    
    auc_pre <- NULL
    brier_pre <- NULL
    for(nv in 1:length(preds_list)){
      
      #The function expects increasing values for bad outcomes (i.e. risk)
      
      # when the predictions are from trees with no nodes, the Score function gives an error
      # To remedy it change one of the predictions very tiny bit - ERROR
      # if(length(unique(preds_list[[nv]][[1]]))==1){
      #   preds_list[[nv]][[1]][1] <- preds_list[[nv]][[1]][1]-0.000000001
      # }
      # 
      td_auc <- Score(preds_list[[nv]],
                      formula=Surv(day, censor)~1, metrics=c("auc", "Brier"), 
                      data = test_dat, 
                      conf.int = F, 
                      se.fit  = F, 
                      times = visits[nv])
      
      auc_pre[nv] <- td_auc$AUC$score$AUC
      brier_pre[nv] <- td_auc$Brier$score$Brier[2]
      
    }

    auc <- rbind(auc, auc_pre)
    brier <- rbind(brier, brier_pre)
    
    allpredictions <- rbind(allpredictions, predictions)
  }
  
  auc <- data.frame(auc)
  colnames(auc) <- colnames(allpredictions)
  
  brier <- data.frame(brier)
  colnames(brier) <- colnames(allpredictions)
  
  #find which variables were important in at least 2/5 folds
  imp_var <- unique(imp_var_folds[duplicated(imp_var_folds)==T])
  
  if(type %in% c("elastic", "grelastic")){
    imp_var <- imp_var[-which(imp_var=="DrugFTY720")]
  }
  
  output <- list(label = lab, type = type, mark = "folds", 
                 tab_imp_var = imp_var, auc = auc, brier = brier, 
                 mean_brier = apply(brier, 2, mean, na.rm = T),
                 mean_auc = apply(auc, 2, mean, na.rm = T), 
                 allpredictions = allpredictions, 
                 fnl_mdl = fnl_mdl)
  
  #calculate the average of variable importances across the folds
  if(type == "forest"){
    avg_varimp <- apply(for_avg_varimp, 1, mean)
    
    #create a dataframe with the average variable importances, which will be the whole dataset for further figures
    var_imp_dat <- data.frame(importance = unname(avg_varimp), 
                              variable  = names(avg_varimp)) %>%
      arrange(importance) %>%
      mutate(rank = 1:nrow(.))
    
    output$dat_var_imp <- var_imp_dat
  
    } 
  return(output)
}

### compute parallel  ####
# test <- loop.one(loop_dat[[3]])
cl <- makeCluster(n_cores-1)
registerDoParallel(cl)

clusterExport(cl, list("loop.one", "ana_dat", "outcome_cat", "path_code", "n_out", "study"))

tic("model predict loop")

loops_all <- c(parLapply(cl, loop_dat, fun = loop.one))

toc()

stopCluster(cl)

save(loops_all, file = paste0(path_results, "loops_all_", ana_type,".RData"))

# load(paste0(path_results, "loops_all_", ana_type,".RData"))

### check models ####

imp_var_per_out <- NULL
auc_pero_per_out <- NULL
brier_pero_per_out <- NULL
auc_pero_6 <- NULL
auc_pero_12 <- NULL
auc_pero_24 <- NULL
auc_smry_per_out <- NULL

brier_pero_per_out <- NULL
brier_pero_6 <- NULL
brier_pero_12 <- NULL
brier_pero_24 <- NULL
brier_smry_per_out <- NULL

for(ou in 1:n_out){
  
  loops_pero <- NULL
  predictions_pero <- NULL
  auc_pero <- NULL
  brier_pero <- NULL
  all_imp_var <- NULL
  auc_smry <- NULL
  brier_smry <- NULL
  
  for(t in 0:(length(type)-1)){
    loops_pero[[t+1]] <- loops_all[[(t*length(nout))+ou]]
    
    if(is.null(predictions_pero)){
      predictions_pero <- loops_pero[[t+1]]$allpredictions[,3]
    } else{
      predictions_pero <- cbind(predictions_pero, loops_pero[[t+1]]$allpredictions[,3])
    }
    

    predictions_pero <- data.frame(predictions_pero)
    
    colnames(predictions_pero)[t+1] <- paste0(type[t+1], "_day_720")
    
    print(type[t+1])
    print(loops_pero[[t+1]]$auc)
    
    auc_smry_pre <- loops_pero[[t+1]]$mean_auc
    auc_smry_pre <- apply(loops_pero[[t+1]]$auc, 2, mean, na.rm = T)
    
    brier_smry_pre <- loops_pero[[t+1]]$mean_brier
    brier_smry_pre <- apply(loops_pero[[t+1]]$brier, 2, mean, na.rm = T)
    
    auc_pero <- cbind(auc_pero, c(auc_smry_pre, mean = mean(auc_smry_pre)))
    brier_pero <- cbind(brier_pero, c(brier_smry_pre, mean = mean(brier_smry_pre)))
    
    # auc_pero <- c(auc_pero, loops_pero[[t+1]]$mean_auc)
    colnames(auc_pero)[t+1] <- type[t+1]
    
    colnames(brier_pero)[t+1] <- type[t+1]
    
    if(type[t+1]=="grelastic"){
      
      pre_list <- loops_pero[[t+1]]$tab_imp_var[order(loops_pero[[t+1]]$tab_imp_var)]
      all_imp_var[[t+1]] <- pre_list[-str_which(pre_list, 
                                        "^DrugFTY720.")]
    } else{
      all_imp_var[[t+1]] <- loops_pero[[t+1]]$tab_imp_var[order(loops_pero[[t+1]]$tab_imp_var)]
    }
    names(all_imp_var)[t+1] <- type[t+1]
  }
  
  print(summary(predictions_pero))
  
  for(ch in 1:(length(type)-1)){
    cat(type[ch], "-", type[ch+1],
        "\ndifference\n", summary(predictions_pero[,ch]-predictions_pero[,ch+1]), 
        "\ncorrelation\n",cor(predictions_pero[,ch], predictions_pero[,ch+1]), "\n\n")
  }
  
  auc_smry_per_out[[ou]] <- auc_pero
  imp_var_per_out[[ou]] <- all_imp_var
  auc_pero_per_out <- cbind(auc_pero_per_out, auc_pero[length(visits)+1,])
  auc_pero_6 <- cbind(auc_pero_6, auc_pero[1,])
  auc_pero_12 <- cbind(auc_pero_12, auc_pero[2,])
  auc_pero_24 <- cbind(auc_pero_24, auc_pero[3,])
  
  
  names(auc_smry_per_out)[ou] <- loops_all[[ou]]$label
  names(imp_var_per_out)[ou] <- loops_all[[ou]]$label
  
  colnames(auc_pero_per_out)[ou] <- loops_all[[ou]]$label
  
  
  
  brier_smry_per_out[[ou]] <- brier_pero
  brier_pero_per_out <- cbind(brier_pero_per_out, brier_pero[length(visits)+1,])
  brier_pero_6 <- cbind(brier_pero_6, brier_pero[1,])
  brier_pero_12 <- cbind(brier_pero_12, brier_pero[2,])
  brier_pero_24 <- cbind(brier_pero_24, brier_pero[3,])
  
  
  names(brier_smry_per_out)[ou] <- loops_all[[ou]]$label
  
  colnames(brier_pero_per_out)[ou] <- loops_all[[ou]]$label
}

tab_auc <- auc_pero_per_out
tab_imp_var <- imp_var_per_out
tab_auc_all <- list(avg = tab_auc, m6 = auc_pero_6, m12 = auc_pero_12, m24 = auc_pero_24)

tab_brier <- brier_pero_per_out
tab_brier_all <- list(avg = tab_brier, m6 = brier_pero_6, m12 = brier_pero_12, m24 = brier_pero_24)

smry_npred <- NULL

for(i in 1:(length(loops_all))){
  spec_type <- loops_all[[i]]$type
  outcome <- loops_all[[i]]$label
  
  if(spec_type=="grelastic"){
    list_grelastic <- loops_all[[i]]$fnl_mdl$model
    
    if(glm_lambda=="lambda.min"){
      gr_lambda <- list_grelastic$model$lambda.min
    } else if(glm_lambda=="lambda.1se"){
      gr_lambda <- list_grelastic$model$lambda.1se
    }
    
    coef_list <- coef(list_grelastic$model, lambda = gr_lambda)
    sel_var <- attr(coef_list, "names")[coef_list!=0]
  } else if(spec_type=="elastic"){
    list_elastic <- loops_all[[i]]$fnl_mdl$model
    
    coef_list <- coef(list_elastic$model, s = glm_lambda)
    sel_var <- attr(coef_list, "Dimnames")[[1]][attr(coef_list, "i")+1]
  } else if (spec_type=="tree"){
    tree_cox <- loops_all[[i]]$fnl_mdl$model$model
    node_var <- NULL
    
    if(length(tree_cox)>1){

      for(n in 1:((length(tree_cox)-1)/2)){
        node_var <- c(node_var, names(info_node(node_party(tree_cox[[n]]))$p.value))
      }
      
      sel_var <- c(rep(unique(node_var), 2), "Drug")
    } else{
      sel_var <- NA # "none"
    }
    
  }
  
  if(spec_type!="forest"){
    smry_npred <- rbind(smry_npred, c(spec_type, outcome, length(sel_var)))
  }
}

smry_npred <- data.frame(smry_npred)
colnames(smry_npred) <- c("type", "outcome", "n_pred")

tab_npred <- pivot_wider(smry_npred, names_from = "type", values_from = "n_pred", names_prefix = "n_pred_")

#### Just for the figures
var_imp_forest <- loops_all[n_out+c(nout)]

varimp.figures <- function(varimplist, scheme = "color"){
  source(paste0(path_code, "/00_paths_pckgs_functions.R"))
  # mlt etc. depends on survival 3.3.1 as opposed to the newer version 3.4, so:
  # library("survival", .libPaths()[2])
  
  if(scheme=="color"){
    colcodes <- other_colcod
  } else if (scheme=="gray"){
    colcodes <- other_graycod
    colcodes[1] <- "#000000"
  }
  
  imp_var <- varimplist$tab_imp_var
  var_imp_dat <- varimplist$dat_var_imp
  lab <- varimplist$label
  
  importvar <- var_imp_dat$importance>abs(min(var_imp_dat$importance))
  
  fig_varimp <- ggplot(var_imp_dat, aes(x = importance, y = rank)) + 
    geom_segment(aes(yend = rank, x=0, xend = importance, color = "others")) + 
    geom_vline(aes(xintercept = 0), colour = colcodes[2]) +
    geom_text(aes(label = "0", fontface = "bold", family="sans"),
                x = 0, size = 3, y = 0.3, vjust = 0, hjust = 0.5, 
              colour = colcodes[1]) + 
    geom_vline(aes(xintercept = min(importance)), colour = colcodes[1], lty = 2) +
    geom_vline(aes(xintercept = abs(min(importance))), colour = colcodes[1], lty = 2) +
    scale_y_continuous(labels = function(x) max(var_imp_dat$rank)-x+1,
                       breaks = c(1,  max(var_imp_dat$rank)-59, max(var_imp_dat$rank)-34, 
                       max(var_imp_dat$rank)-14, max(var_imp_dat$rank)),
                       trans = squash_axis(0, 
                                           max(var_imp_dat$rank)-length(which(importvar==T)), 
                                           4)) +
    coord_cartesian(ylim = c(0, (max(var_imp_dat$rank)+1)), 
                    xlim = c((min(var_imp_dat$importance)*1.2), 
                             ((max(var_imp_dat$importance)+abs(min(var_imp_dat$importance)))*1.3)), 
                    expand = F)

    if(length(imp_var)!=0){
      fig_varimp <- fig_varimp +
        geom_segment(data = ~ subset(var_imp_dat, (variable %in% imp_var)),
                     aes(yend = rank, x= 0, xend = importance, y = rank, 
                         color = "important (min twice)"),
                     size = 1.5) + 
        geom_text(data = ~ subset(var_imp_dat, (variable %in% imp_var) & importance>=0 & importance<=abs(min(importance))),
                  aes(label = str_replace(variable, "_", " ")),
                  nudge_x = 0.1, size = 2, hjust = 0, color =  colcodes[6]) +
        geom_text(data = ~ subset(var_imp_dat, (variable %in% imp_var) & importance<0 & importance<=abs(min(importance))),
                  aes(label = str_replace(variable, "_", " ")),
                  x = 0.1, size = 2, hjust = 0, color =  colcodes[6])
    } 
  
  if(any(importvar)){
    fig_varimp <- fig_varimp +
      geom_segment(data = ~ subset(var_imp_dat, (importance>abs(min(importance)))),
                   aes(yend = rank, x= 0, xend = importance, y = rank, 
                       color =  "important (by average)"),
                   size = 1.5) +
      scale_colour_manual(name = NULL, 
                          values = c("important (by average)" = colcodes[7], 
                                     "important (min twice)" = colcodes[6],
                                     "others" = colcodes[3]), 
                          guide = guide_legend()) + 
      geom_text(data = ~ subset(var_imp_dat, (importance>abs(min(importance)))),
                aes(label = str_replace(variable, "_", " ")),
                nudge_x = 0.1, size = 2, hjust = 0, color = colcodes[7]) 
  }
  
  fig_varimp <-  fig_varimp +

    labs(title = lab, 
         x = "log-likelihood based permutation importance (average)") +
    theme(legend.position = c(0.99, 0.01), 
          legend.justification = c(1,0),
          legend.spacing.y = unit(0.1, "lines"),
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          text = element_text(family = "sans", size = 7))

  return(fig_varimp)
}

# create parallel computation for figures for all

cl <- makeCluster(n_cores-1)
registerDoParallel(cl)

clusterExport(cl, list("varimp.figures", "path_code", "study"))

tic("figures")

fig_var_imp <- c(parLapply(cl, var_imp_forest, fun = varimp.figures))
fig_var_imp_gray <- c(parLapply(cl, var_imp_forest, fun = varimp.figures, scheme = "gray"))

toc()

stopCluster(cl)

loops_tree <- loops_all[1:6]

plot.new()

tree.figures <- function(loops_out, scheme="color"){

  tree_cox <- loops_out$fnl_mdl$model
  
  if(scheme=="color"){
    colcodes <- drug_colcod
  } else if (scheme=="gray"){
    colcodes <- drug_graycod
  }
  
  plot.new()
    
    plot(tree_cox$model, newdata = ana_dat, 
         tp_args = list(type = "survivor", col = colcodes[c(2, 1)], 
                        gp = gpar(lwd = 1.5, lty = c(1, 2)), fontsize = 5), 
         gp = gpar(family = "mono", fontsize = 6, cex = 1,
        lwd = 1, xlog = F, xaxp = c(0, 720, 5)))
    
    myplot <- grid.grab()
    
    plot.new()
    
    myleg <- legendGrob(labels = c("Drug=FTY720", "Drug=Placebo"), 
                        vgap = unit(1, "scaledpts"), 
                        gp=gpar(lwd = 1.5, 
                                col =colcodes, 
                                fontsize = 6, bty = "n"))
    
    mytit <- textGrob(label = loops_out$label, hjust = 0,
                      x = unit(0.1, "npc"), y = unit(0.5, "npc"),
              gp=gpar(fontsize = 8, fontface = 2))
    
    myfig <- grid.arrange(mytit, myplot, myleg, ncol = 6, 
                 layout_matrix = rbind(c(1,1,1,1,3,3), 
                                       rep(2,6), rep(2,6), rep(2,6), 
                                       rep(2,6), rep(2,6), rep(2,6)))
    return(myfig)

}

fig_tree <- lapply(loops_tree, tree.figures)
fig_tree_gray <- lapply(loops_tree, tree.figures, scheme = "gray")

# choose models (manual)####
best_mdls <- data.frame(best_mdl = 
                          c("elastic", "elastic", "forest", "forest", "elastic", "grelastic"), 
                        outcome = colnames(tab_auc))

best_mdls0 <- str_c(best_mdls$outcome, best_mdls$best_mdl)

fnls <- NULL

for(i in 1:length(loops_all)){
  comb <- paste0(loops_all[[i]]$label,  loops_all[[i]]$type)
  
  if(comb %in% best_mdls0){
    add <- list(model = loops_all[[i]]$fnl_mdl$model,
                best_hyp = loops_all[[i]]$fnl_mdl$best_hyp,
                label = loops_all[[i]]$label,
                type = loops_all[[i]]$type)
    fnls <- c(fnls, 
              list(add))
  }
}

# save the tables to be reported
save(list = ls()[str_which(ls(), "^tab_")], 
                 file = paste0(path_results, "loops_tab_", ana_type,".RData"))

# save the figures to be reported
save(list = ls()[str_which(ls(), "^fig_")], 
     file = paste0(path_results, "loops_fig_", ana_type,".RData"))
