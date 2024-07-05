################################################################################
# Collection of paths, packages, functions - to be sourced
################################################################################

#### PATHS

# The project must be started
path_main <- getwd() ## must be if the project
path_code <- paste0(path_main, "/Code/")

path_data <- paste0(path_main, "/", study, "/Data/")
path_results <- paste0(path_main, "/", study, "/Results/")

source(paste0(path_code, "01_studyparam.R"))

file_form <- ".sas7bdat"

#### PACKAGES
# devtools::install_version("tram", "0.6-4")

.libPaths(c("D:/Users/mse0ibo/Documents/myRlib",
            "D:/Users/mse0ibo/AppData/Local/Programs/R/R-4.2.0/library"))

# install.packages("TreatmentSelection", lib = .libPaths()[1])

pckgs <- c("haven", "survival", "basefun", "survivalROC", "risksetROC", "Hmisc", "gtsummary", "grpreg",
           "tidyverse", "caret", "survminer", "missForest", "tram", "trtf", "parallel", "doParallel", 
           "tictoc", "glmnet", "mice", "riskRegression", "rio", "grid", "gridExtra", "survminer", 
           "scales", "gridGraphics", "ggplotify", "cowplot", "TreatmentSelection", "dcurves")

lapply(pckgs, library, character.only = TRUE)

theme_set(theme_bw())

drug_colcod <- hue_pal()(2)
drug_graycod <- gray.colors(7)[c(1,5)]

other_colcod <-  brewer_pal(palette = "PiYG")(7)
other_graycod <- gray.colors(9)[8:2]
  
### USER_DEFINED FUNCTIONS

#' Describe numeric variables
#' @param x a numeric vector to be described
#' @return a vector of n valid, n missing, prop missing, min, 1stQ, median, 3rdQ, max, IQR, mean, sd, 

desc.num <- function(x){
  if(!is.numeric(x)) {
    print("Error: The vector is not numeric")
  }
  else {
    n_miss <- length(which(is.na(x)))
    n_val <- length(x)-n_miss
    qua <- quantile(x, na.rm = T)
    desc_num <- round(c(n = n_val, n_miss = n_miss, prop_miss = 100*(n_miss/length(x)),
                  min = qua[1], Q1 = qua[2], med = qua[3], Q3 = qua[4], max = qua[5], iqr = qua[4]-qua[2],
                  mean = mean(x, na.rm = T), sd = sd(x, na.rm = T)), 1)
  }
}

#' Describe factor variables
#' @param x a factor vector to be described
#' @return a data frame of n valid (var level), n missing (var level), n category, prop category

desc.cat <- function(x){
  # if(!is.factor(x)) {
  #   as.factor(x)
  # }
    n_miss <- length(which(is.na(x)))
    n_val <- length(x)-n_miss
    cat_lab <- levels(x)
    n <- length(cat_lab)
    
    n_cat <- NULL
    category <- NULL
    for(i in 1:(n-1)){
      n_cat[i] <- length(which(x==cat_lab[i+1]))
      category[i] <- paste(cat_lab[i+1], "(ref.", cat_lab[1], ")")
    }
    
    desc_fac <- data.frame(n = c(n_val, rep(NA, n-2)), n_miss = c(n_miss, rep(NA, n-2)), 
                           prop_miss = c(round(100*(n_miss/length(x))), rep(NA, n-2)), 
                           category = category, n_cat = n_cat, n_prop = round(100*(n_cat/length(x))))
}

#' Summarize all variables, also by subgroup if needed
#' @param dat dataframe to be summarized
#' @param by subgrouping factor variable as character, defaults to NULL
#' @return a list of data summaries - both numeric and factors, by subgroup if provided

smry.all <- function(dat, by = NULL){
  if(is.null(by)){
    out <- smry.per(dat)
  } else{
    fact <- get(by, dat)
    if(!is.factor(fact))
      fact <- factor(fact)
    dat <- dat[,-which(colnames(dat)==by)]
    out <- NULL
    for(i in levels(fact)){
      dat_by <- dat[fact==i & !is.na(fact),]
      out_pre <- smry.per(dat_by)
      out[[length(out)+1]] <- out_pre
    }
    names(out) <- levels(fact)
    
    p_value <- NULL
    for(k in 1:ncol(dat)){
      if(is.numeric(dat[,k])){
        p <- wilcox.test(dat[,k]~fact)$p.value
      }
      else{
        p <- chisq.test(dat[,k], fact)$p.value
      }
      p_value[k] <- p
    }
    names(p_value) <- colnames(dat)
    out$p_wilcox_chisq <- p_value
  }
  return(out)
}

#' Internal function for smry.all to aummarize all variables
#' @param dat dataframe to be summarized
#' @return a list of data summaries - both numeric and factors

smry.per <- function(dat){
  dat_num <- select(dat, where(is.numeric))
  dat_cat <- data.frame(select(dat, where(is.factor)))
  
  smry_num_pre <- apply(dat_num, 2, desc.num)
  
  
  smry_cat_pre <- NULL
  for(ncat in 1:ncol(dat_cat)){
    smry_cat_pre[[ncat]] <- desc.cat(dat_cat[,ncat])
  }
  # smry_cat_pre <- apply(dat_cat, 2, desc.cat)
  
  smry_cat <- NULL
  
  # print(smry_cat_pre)
  
  for(i in 1:length(smry_cat_pre)){
    
    l <- attr(dat_cat[[i]], which = "label")
    label_cat <- ifelse(is.null(l), colnames(dat_cat[i]), l)
    
    smry_cat <- rbind(smry_cat, 
                      cbind(labels = rep(label_cat, nrow(smry_cat_pre[[i]])),
                            var = rep(colnames(dat_cat)[i], nrow(smry_cat_pre[[i]])), 
                            smry_cat_pre[[i]]))
  }
  
  # print(smry_cat)
  
  label_num <- NULL
  for(i in 1:ncol(dat_num)){
    l <- attr(dat_num[[i]], which = "label")
    label_num[i] <- ifelse(is.null(l), colnames(dat_num[i]), l)
  }
  smry_num <- data.frame(labels = label_num, var = colnames(dat_num), t(smry_num_pre))
  
  return(list(smry_num = smry_num, smry_cat = smry_cat))
}

#' log-likelihood permutation variable importances for a traforest in a new dataset
#' @param object a traforest 
#' @param dat dataframe for which variable importances is going to be calculated 
#' @return a vector of log-likelihood permutation variable importances in the (preferably test) dataset

vimp.cv <- function(object, dat, nperm = 1){
  baselik <- logLik(object, newdata = dat)
  tr <- which(names(dat)=="Drug")
  permimps <- sapply(names(dat)[-c(1, tr)], function(x){
    permdat <- dat
    log_lik <- NULL
    for(i in 1:nperm){
      set.seed(35+i)
      permdat[,x] <- sample(permdat[,x])
      log_lik[i] <- logLik(object, newdata = permdat)
    }
    return(mean(log_lik))
  })
  return(baselik-permimps)
}

#' develop the models based on type
#' @param dat dataframe for which to fit the model 
#' @param seed seed for procedures that include randomness
#' @param type the model type. one of "forest", "tree", "elastic", or "grelastic"
#' @param final in order to accomodate the need for more reporting of final models
#' @param minbucket/ntree/stump tuning parameters for the model type = "forest"
#' @param minbucket/test_alpha/inc_test_alpha tuning parameters for the model type = "tree"
#' @param range_alpha/inc_alpha tuning parameters for the model types = "elastic", grelastic
#' @param m number of missing imputations. Currently can only be used with "elastic" - "grelastic does not allow it"
#' @return the final model


develop.mdl <- function(dat, seed=1, type, final = F, 
                        minbucket = c(5, 10, 15),
                        stump = F, ntree = 100,
                        test_alpha = c(0.05, 0.2), inc_test_alpha = 0.05, 
                        range_alpha = c(0.1, 1), inc_alpha = 0.1, m = 1){
  
  if(type %in% c("tree", "forest")){
    dat <- imp.fortree(dat, seed = seed)
    base_model <- Coxph(outcome ~ Drug, data = dat)
    
    if(stump == T){
      if(type =="forest"){
        trgt_model <- traforest(base_model,
                                formula = outcome | Drug ~ .,
                                data = dat, 
                                stump = stump)
      } else{
        trgt_model <- trafotree(base_model,
                                formula = outcome | Drug ~ .,
                                data = dat, 
                                stump = stump)
      }
    } else{
      set.seed(seed)
      partFolds <- createFolds(y = dat$Drug, k = n_folds)
      
      if(type=="forest"){
        mtry <-   c(floor(sqrt(ncol(dat)-2)), ceiling((ncol(dat)-2)/3))
        hyp <- expand.grid(mtry = mtry, minbucket = minbucket)
      } else if(type == "tree"){
        hyp <- expand.grid(alpha = seq(test_alpha[1], 
                                       test_alpha[2], 
                                       by = inc_test_alpha), 
                           minbucket = minbucket)
      }
      
      for(nf in 1:n_folds){
        tstdat <- dat[partFolds[[nf]],]
        trdat <- dat[-partFolds[[nf]],]
        
        tra_base_model <- Coxph(outcome ~ Drug, data = trdat)
        
        if(type == "forest"){
          set.seed(seed+nf)
          hyp[,nf+2] <- apply(hyp[,1:2], 1, function(x){
            tra_model <- traforest(tra_base_model, formula = outcome | Drug ~ ., data = trdat,
                                   ntree = ntree, mtry = x[1], minbucket = x[2])
            lik <- logLik(tra_model, newdata = tstdat)
            
            return(lik)})
          
        } else if (type == "tree"){
          hyp[,nf+2]  <- apply(hyp[,1:2], 1, function(x){
            tra_model <- trafotree(tra_base_model, formula = outcome | Drug ~ ., data = trdat, 
                                   control = ctree_control(alpha = x[1], minbucket = x[2]))
            lik <- logLik(tra_model, newdata = tstdat)
            return(lik)})
        }
      } 
      
      hyp$avg_hyp <- apply(hyp[,3:ncol(hyp)], 1, mean)
      best_hyp <- hyp[which(hyp$avg_hyp == max(hyp$avg_hyp)),][1,1:2]
      
      if(type == "forest"){
        trgt_model_pre <- traforest(base_model, formula = outcome | Drug ~ ., data = dat, ntree = ntree, 
                                    mtry = best_hyp$mtry, minbucket = best_hyp$minbucket)
        } else if (type == "tree"){
        trgt_model_pre <- trafotree(base_model, formula = outcome | Drug ~ ., data = dat,  
                                control = ctree_control(alpha = best_hyp$alpha, 
                                                        minbucket =  best_hyp$minbucket))
        }
      trgt_model <- list(model = trgt_model_pre, data = dat)
      }
    } else if (type %in% c("elastic", "grelastic")){
    mdl_mtrx_dat <- imp.forglm(dat, seed = seed, m = m)
    
    set.seed(seed)
    lassoFolds <- createFolds(y = dat$Drug, k = n_folds)
    
    alphalist <- seq(range_alpha[1], range_alpha[2], by = inc_alpha)
    
    if(m > 1){
      multi_folds <- lapply(lassoFolds, function(x){
        ids <- x
        fnl_ids <- NULL
        for(ind in 0:(m-1)) {
          fnl_ids <- c(fnl_ids, ids+nrow(dat)*ind)}
        return(fnl_ids)})
      
      foldid_pre <- NULL
      for(n_f in 1:n_folds){
        ids  <- multi_folds[[n_f]] 
        foldid_pre <- rbind(foldid_pre, cbind(ids, fold = rep(n_f, length(ids))))
      }
      foldid_pre <- data.frame(foldid_pre)
      foldid_mtrx <- foldid_pre[order(foldid_pre$ids),]
      rownames(foldid_mtrx) <- 1:nrow(foldid_mtrx)
    }
    if(type == "elastic"){
      choose_alpha <- lapply(alphalist, function(a){
        cv.glmnet(x=as.matrix(mdl_mtrx_dat[,-c(1:3)]), 
                  y = mdl_mtrx_dat$outcome, 
                  # the penalty factor is internally rescaled - just doing it manually
                  penalty.factor = c(0, rep((ncol(mdl_mtrx_dat)-3)/(ncol(mdl_mtrx_dat)-4), 
                                     (ncol(mdl_mtrx_dat)-4))),
                  family = "cox", alpha = a, weights = mdl_mtrx_dat$fit_weight, 
                  nfolds = n_folds, foldid = foldid_mtrx$fold)
      })
      
      perf_alpha <- NULL
      for(i in 1:length(alphalist)){
        perf_alpha[i] <- min(choose_alpha[[i]]$cvm)
      }
      
      best_a <- which(min(perf_alpha)==perf_alpha)
      trgt_model_pre <- choose_alpha[[best_a]]
      trgt_model_pre$call$alpha <- alphalist[best_a]
      
      trgt_model <- list(model = trgt_model_pre, data = mdl_mtrx_dat)
      
      best_hyp <- c(alpha = alphalist[best_a], lambda.1se = trgt_model_pre$lambda.1se, 
                    lambda.min = trgt_model_pre$lambda.min)
      
    } else if(type == "grelastic"){
      mdl_matrix_dat_re <- mdl_mtrx_dat
      
      # create the variable groups
      lag_order_i <- (min(str_which(colnames(mdl_matrix_dat_re), "DrugFTY720."))-5)
      
      lag_order <- NULL
      for(ord in 1:lag_order_i){
        lag_order <- c(lag_order, ord, ord+lag_order_i)
      }
       
      dat_order <- c(1:4, 4+lag_order)
      mdl_matrix_dat_re <- mdl_matrix_dat_re[,dat_order]
      
      groups <- c("0", rep(colnames(mdl_matrix_dat_re)[(1:lag_order_i)*2+3], 2))
      groups_order <- c(1, 1+lag_order)
      groups <- groups[groups_order]
      groups[str_which(groups, "Age")] <- "Age"
      groups[str_which(groups, "VisAcuity")] <- "VisAcuity"
      
      all <- data.frame(groups = groups)
      
      id <- data.frame(groups = unique(groups))
      id$trgt <- 1:nrow(id)-1
      
      fnl_groups <- inner_join(all, id)$trgt

      choose_alpha <- lapply(alphalist, function(a){
        cv.grpsurv(X=as.matrix(mdl_matrix_dat_re[1:nrow(dat),-c(1:3)]), 
                  y = mdl_matrix_dat_re$outcome[1:nrow(dat)], 
                  group = fnl_groups, penalty = "grLasso", nfolds = n_folds, alpha = a) 
        ### the folds argument does not work somehow
      })
      track_lambda.1se <- lapply(choose_alpha, function(x){
        max(x$lambda[which(x$cve<x$cvse[which(min(x$cve)==x$cve)]+min(x$cve))])
      })
      
      track_lambda.min <- lapply(choose_alpha, function(x){x$lambda.min})
      
      perf_alpha <- NULL
      for(i in 1:length(alphalist)){
        perf_alpha[i] <- min(choose_alpha[[i]]$cve)
      }
      
      best_a <- which(min(perf_alpha)==perf_alpha)
      
      lambda.1se <- track_lambda.1se[[best_a]]
      lambda.min <- track_lambda.min[[best_a]]
      
      trgt_model_pre <- grpsurv(X=as.matrix(mdl_matrix_dat_re[1:nrow(dat),-c(1:3)]), 
                            y = mdl_matrix_dat_re$outcome[1:nrow(dat)], 
                            group = fnl_groups, penalty = "grLasso", alpha = alphalist[best_a])
      trgt_model_pre$lambda.1se <- lambda.1se
      trgt_model_pre$lambda.min <- lambda.min
      trgt_model <- list(model = trgt_model_pre, data = mdl_matrix_dat_re)
      
      best_hyp <- c(alpha = alphalist[best_a], 
                    lambda.1se = lambda.1se, 
                    lambda.min = lambda.min)
    }
  }
    if(final==T){
      return(list(model = trgt_model, best_hyp = best_hyp))
    } else {
    return(trgt_model)
  }
}

#' impute the dataset and return an m-times stacked dataset ready for analysis/prediction with weights
#' @param dat dataframe for which to impute
#' @param m number of imputed datasets to stack
#' @param seed seed for procedures that include randomness
#' @return imputed stacked datasets

imp.forglm <- function(dat,  m = 5, seed = 1){
  
  y <- rep(with(dat, Surv(time = day, event = censor)), m)
  y_na <- nelsonaalen(dat, day, censor)
  for_weight <- data.frame(select(dat, -day, -censor))
  
  #No missing in drug, so take that into account, also, consider the interactions
  fit_weight <- rep(apply(for_weight, 1, function(x) {
    ((length(which(is.na(x)==F))-1)*2+1)/(((length(x)-1)*2+1)*m)}), 
    m)
  
  for_weight$day <- dat$day
  for_weight$y_na <- y_na
  
  #impute with the outcome
  
  # md.pattern(for_weight, rotate.names = T)
  dat_imp_all <- mice(for_weight, m = m, seed = seed, maxit = 20, printFlag = F)
  
  #stack the datasets
  dat_imp <- NULL
  for(i in 1:m){
    dat_imp <- bind_rows(dat_imp, complete(dat_imp_all, i))
  }
  
  dat_imp <- data.frame(select(dat_imp, -day, -y_na))
  
  ##extract the data matrix
  
  mdl_mtrx_dat <- data.frame(outcome = y , fit_weight = fit_weight,
                             model.matrix(rep(dat$day, m) ~ .*Drug, data = dat_imp))
  return(mdl_mtrx_dat)
}

imp.fortree <- function(dat, seed = seed){
  y <- with(dat, Surv(time = day, event = censor))
  y_na <- nelsonaalen(dat, day, censor)
  
  mdl_dat_pre <- select(dat, -day)
  mdl_dat_pre$y_na <- y_na
  
  set.seed(seed)
  
  imp_dat <- suppressWarnings(missForest(mdl_dat_pre))
  imp_dat_fnl <- select(imp_dat$ximp, -censor, -y_na)
  imp_dat_fnl$outcome <- y
  
  return(imp_dat_fnl)
}

#' perform predictions based on method
#' @param object model for predictions
#' @param dat dataframe for which to perform predictions
#' @param targett the target time to calculate the predictions
#' @param type type of model
#' @glm_lambda which lambda to choose for type = elastic or grelastic. options are "1se" or "min"
#' @pred_type type of predictions for type = elastic or grelastic. options are "prob" or "lp"
#' @return imputed stacked datasets

predict.all <- function(object, dat, targett=1, glm_lambda = NULL, type, pred_type = "prob"){
  
  predictions_df <- data.frame(rep(NA, nrow(dat)))
  co <- 1
  if(type == "forest"){
    predictions_frpre <- predict(object, newdata = dat, mnewdata = dat, type = "survivor", q = targett)
  } else if (type == "tree"){
    predictions_trpre <- predict(object, newdata = dat, type = "survivor", q = targett)
  } else if (type == "grelastic"){
    mdl_mtrx_dat <- object$data
    grelastic <- object$model
    
    dat <- dat %>% select(colnames(mdl_mtrx_dat))
    
    if(pred_type == "lp"){
      predictions <- predict(grelastic, as.matrix(dat[,-c(1:3)]), type = "link", lambda = glm_lambda)
      return(predictions)
    } else if (pred_type == "prob"){
      predictions_grpre <- NULL
      for(pat in 1:(nrow(dat))){
        predictions_grpre <- rbind(predictions_grpre, 
                                   predict(grelastic, as.matrix(dat[pat,-c(1:3)]), 
                                           type="survival", lambda = glm_lambda)(targett))
      }
    }
  }
  
  for(tt in targett){
    if(type == "forest"){
      predictions <- unlist(lapply(predictions_frpre, function(x)x[co,1]))
    } else if (type == "tree"){
      predictions <- predictions_trpre[(0:(nrow(dat)-1))*length(targett)+co]
    } else if (type == "elastic"){
      mdl_mtrx_dat <- object$data
      elastic <- object$model
      
      if(pred_type == "lp"){
        predictions <- predict(elastic, s = glm_lambda, newx = as.matrix(dat[,-c(1:3)]), type = "link")
        return(predictions)
      } else if (pred_type == "prob"){
        all_predictions <- survfit(elastic, s = glm_lambda, x=as.matrix(mdl_mtrx_dat[,-c(1:3)]),
                                   y = mdl_mtrx_dat$outcome, weights = mdl_mtrx_dat$fit_weight,
                                   newx = as.matrix(dat[,-c(1:3)]))
        x <- tt - all_predictions$time
        time_ind <- max(which(x>=0))
        predictions <- all_predictions$surv[time_ind,]
      }
    } else if (type == "grelastic"){
      predictions <- predictions_grpre[,co]
    }
    predictions_df[,co] <- predictions
    colnames(predictions_df)[co] <- paste0("day_", tt)
    co <- co+1
  }
  return(data.frame(predictions_df))
}

calculate_ci <- function (est, se, q = 1.96){
  ci <- cbind(ci_low = est-se*q, ci_hi = est+se*q)
  return(ci)
}

squash_axis <- function(from, to, factor){
  trans <- function(x){
    isq <- x>from & x<to
    ito <- x>=to
    
    x[isq] <- from + (x[isq]-from)/factor
    x[ito] <- from + (to-from)/factor + (x[ito]-to)
    return(x)
  }
  inv <- function(x){
    isq <- x>from & x<from + (to-from)/factor
    ito <- x>=from + (to-from)/factor
    
    x[isq] <- from + (x[isq]-from)*factor
    x[ito] <- to + (x[ito]- (from + (to-from)/factor))
    return(x)
  }
  return(trans_new("squash_axis", trans, inv))
}

### For some reason, when parallelized, R cannot find the following auxillary functions to basefun::Bernstein_basis
### copy-pasted from Github page basefun/R/Bernstein.R

.Bx <- function (x, j, n, deriv = 0L, integrate = FALSE) 
{
  if (deriv == 0L) 
    return(.B(j, n, integrate = integrate)(x))
  return(n * (.Bx(x, j - 1, n - 1, deriv = deriv - 1L, integrate = integrate) - 
                .Bx(x, j, n - 1, deriv = deriv - 1L, integrate = integrate)))
}

.B <- function (j, n, integrate = FALSE) 
{
  function(x) {
    if (j < 0 || (n - j) < 0) 
      return(rep(0, length(x)))
    if (integrate) 
      return(pbeta(x, shape1 = j + 1, shape2 = n - j + 
                     1)/(n + 1))
    return(dbeta(x, shape1 = j + 1, shape2 = n - j + 1)/(n + 
                                                           1))
  }
}

outnames_convert <- function(outcome_cat, col_names){
  fnl_names <- NULL
  for(i in 1:nrow(outcome_cat)){
    col_names <-   str_replace(col_names, paste0("outcome_", outcome_cat$pattern[i]), 
                               outcome_cat$short_labels[i])
  }
  return(col_names)
}
