################################################################################
# Collection of functions relevant for reporting - to be sourced
################################################################################

### for reporting baseline description
make_desc_exp <- function(desc_dat){
  num_header <- NULL
  cat_header <- NULL
  num <- NULL
  cat <- NULL
  pvalue <- NULL
  
  if(names(desc_dat[1])=="smry_num"){
    numerical = make_num_great(desc_dat[[1]])
    categorical = make_cat_great(desc_dat[[2]])
    
  } else if(names(desc_dat[[1]][1])=="smry_num"){
    # do this and that
    for(l1 in 1:(length(desc_dat)-1)){
      if(l1==1){
        numnew <- make_num_great(desc_dat[[l1]][[1]])
        catnew <- make_cat_great(desc_dat[[l1]][[2]])
        
        num <- numnew
        cat <- catnew
      } else{
        numnew <- make_num_great(desc_dat[[l1]][[1]])[,-1]
        catnew <- make_cat_great(desc_dat[[l1]][[2]])[,-1]
        
        num <- cbind(num, numnew)
        cat <- cbind(cat, catnew)
      }
      
      num_header <- c(num_header, rep(names(desc_dat[l1]), ncol(numnew)))
      
      cat_header <- c(cat_header, rep(names(desc_dat[l1]), ncol(catnew)))
    }
    
    pvalue <-  data.frame(variable = names(desc_dat[[length(desc_dat)]]), p = desc_dat[[length(desc_dat)]])
    
    p <- NULL
    for(varname in desc_dat[[1]][[1]]$var){
      p <- c(p, pvalue$p[which(pvalue$variable==varname)])
    }
    
    num$p <- ifelse(p>0.05, "p>0.05", ifelse(p>0.01, "p<0.05", ifelse(p>0.001, "p<0.01", "p<0.001")))
    
    p <- NULL
    for(varname in desc_dat[[1]][[2]]$var){
      p <- c(p, pvalue$p[which(pvalue$variable==varname)])
    }
    
    cat$p <- ifelse(p>0.05, "p>0.05", ifelse(p>0.01, "p<0.05", ifelse(p>0.001, "p<0.01", "p<0.001")))
    
    numerical <- rbind(c(num_header, "NA"), num)
    categorical <- rbind(c(cat_header, "NA"), cat)
    
  } else {
    for(l2 in 1:(length(desc_dat))){
      for(l1 in 1:(length(desc_dat[[l2]])-1)){
        if(l2==1 & l1==1){
          new_num <- make_num_great(desc_dat[[l2]][[l1]][[1]])
          new_cat <- make_cat_great(desc_dat[[l2]][[l1]][[2]])
          
          num <- new_num
          cat <- new_cat
        } else{
          new_num <- make_num_great(desc_dat[[l2]][[l1]][[1]])[,-1]
          new_cat <- make_cat_great(desc_dat[[l2]][[l1]][[2]])[,-1]
          
          num <- cbind(num, new_num)
          
          nr_curr<- nrow(cat)
          nr_new <- nrow(new_cat)
          
          while(nr_curr>nr_new){
            new_cat <- rbind(new_cat, rep(NA, ncol(new_cat)))
            nr_new <- nrow(new_cat)
          }
          cat <- cbind(cat, new_cat)
        }
        cat_header <- c(cat_header, 
                        rep(paste0(names(desc_dat[l2]), "=", names(desc_dat[[l2]][l1])), 
                            ncol(new_cat)))
        
        num_header <- c(num_header, 
                        rep(paste0(names(desc_dat[l2]), "=", names(desc_dat[[l2]][l1])), 
                            ncol(new_num)))
        
        pvalue <-  data.frame(variable = names(desc_dat[[l2]][[length(desc_dat[[l2]])]]), 
                              p = as.numeric(desc_dat[[l2]][[length(desc_dat[[l2]])]]))
        
        pvalue$p <-  ifelse(pvalue$p>0.05, "p>0.05", 
                            ifelse(pvalue$p>0.01, "p<0.05", 
                                   ifelse(pvalue$p>0.001, "p<0.01", 
                                          "p<0.001")))
      }
      p <- NULL
      for(varname in desc_dat[[l2]][[l1]][[1]]$var){
        p <-  c(p, pvalue$p[which(pvalue$variable==varname)])
      }
      num <- cbind(num, p)
      num_header <- c(num_header, NA)
      
      p <- NULL
      for(varname in desc_dat[[l2]][[l1]][[2]]$var){
        p <- c(p, pvalue$p[which(pvalue$variable==varname)])
      }
      cat <- cbind(cat, p)
      cat_header <- c(cat_header, NA)
    }
    numerical <- rbind(num_header, num)
    categorical <- rbind(cat_header, cat)
  }
  
  output <- rbind(numerical, categorical) %>%
    data.frame() %>%
    distinct()  %>%
    left_join(select(pred_namecat, label, order, category), by= c("Characteristic" = "label")) %>%
    mutate(order=ifelse(is.na(order), 0, order)) %>%
    arrange(order) %>%
    select(-order)
  return(output)
}

make_num_great <- function(dat_num) {
  dat_num_fnl <- dat_num %>%
    transmute(Characteristic = labels, 
              "Total_n" = max(n, na.rm = T), 
              "Missing_n_%" = str_c(n_miss, " (", round(as.numeric(prop_miss)), ")"), 
              "Median_IQR_range_Category" = str_c(med.50., "(", Q1.25., "-",Q3.75., ", ", 
                                                  min.0., "-", max.100., ")"), 
              "Mean_sd_n_%" = str_c(mean, " (", sd, ")")) %>%
    data.frame()
  return(dat_num_fnl)
}

make_cat_great <- function(dat_cat) {
  dat_cat_fnl <- dat_cat %>%
    transmute(Characteristic = labels, 
              "Total_n" = max(n, na.rm = T), 
              "Missing_n_%" = str_c(n_miss, " (", prop_miss, ")"), 
              "Median_IQR_range_Category" = category, 
              "Mean_sd_n_%" = str_c(n_cat, " (", n_prop, ")")) %>%
    data.frame()
  return(dat_cat_fnl)
}

make_epv <- function(dat_out){
  dat_out <- desc_by_outcome
  
  out_codes <- NULL
  event_n <- NULL
  
  for(i in 1:length(dat_out)){
    out_codes[i] <- names(dat_out[i])
    event_n[i] <- dat_out[[i]]$`1`$smry_num[1, 3]
  }
  pvar1 <- sum(as.numeric(desc_predcat$Number))+7
  pvar2 <- (pvar1-1)*2+1
  
  epv <- data.frame(outcome = out_codes, n_event = event_n, epv_max = event_n/pvar1,
                    epv_min = event_n/pvar2)
  return(epv)
  
}

make_varimp <- function(list_varimp){
  all_var_list <- NULL
  for(i in 1:length(list_varimp)){
    var_list <- NULL
    new_var_list <- NULL
    new <- 0
    for(a in 1:length(list_varimp[[i]])){
      old <- new
      old_var_list <- new_var_list
      
      if(length(list_varimp[[i]][[a]])==0){
        list_varimp[[i]][[a]] <- "NA"
      }
      new <- length(list_varimp[[i]][[a]])
      var_list <- list_varimp[[i]][[a]]
      
      pred_namecat_use <-  pred_namecat
      
      for(p in 1:nrow(pred_namecat)){
        for(d in 1:nrow(dict)) {
          if(str_detect(pred_namecat_use$name[p], paste0("^", dict$pre[d], "_"))){
            pred_namecat_use$label[p] <- str_c(dict$domain[d], ": ", pred_namecat_use$label[p])
          }
        }
        
        var_list <- str_replace(var_list, paste0(pred_namecat_use$name[p], "(?!_)"), 
                                pred_namecat_use$label[p])
      }
      
      var_list <- str_replace(var_list, "\\.{3}", "-")
      var_list <- str_replace(var_list, "Yes", "=Yes")
      var_list <- str_replace(var_list, "^DrugFTY720", "_DrugFTY720")
      
      var_list <- var_list[order(var_list)]
      
      if(old == 0){
      } else if(old>=new){
        var_list <- c(var_list, rep(NA, old-new))
      } else{
        old_var_list <- apply(old_var_list, 2, function(x) c(x, rep(NA, new-old)))
      }
      new_var_list <- cbind(old_var_list, var_list)
    }
    new_var_list <- data.frame(new_var_list)
    colnames(new_var_list) <- c(names(list_varimp[[i]]))
    all_var_list[[i]] <- new_var_list
    names(all_var_list)[i] <- paste0(str_split_fixed(str_replace(names(list_varimp)[i], 
                                                                 "\\/", "_"), 
                                                     "\\s", n=3)[c(1,2)], collapse = " ")
  }
  return(all_var_list)
}

make_aucbr <- function(dattab){
  refn <- ncol(dattab)
  dattab$day <- dattab$day/30
  colnames(dattab) <- str_replace(colnames(dattab), "day", "month")
  
  dattab2 <- dattab
  dattab2[,(refn-2):refn] <- round(dattab[,(refn-2):refn], 3)
  
  estwci <- apply(dattab2, 1, function(x){
    y <- paste0(x[refn-2], " (", x[refn-1], "-", x[refn], ")")
    return(y)
  })
  
  new_dattab <- data.frame(dattab2[,1:(refn-3)], estwci)
  
  fnl_dat_pre <- pivot_wider(new_dattab, names_from = month, values_from = estwci, 
                         names_prefix = paste0(colnames(dattab)[(refn-2)], " (95% CI) month "))
  
  if(colnames(dattab)[(refn-2)]=="Brier"){
    
    sc_dat_pre <- pivot_wider(dattab[,1:(refn-2)], names_from = model, values_from = Brier)
    sc_dattab_pre <- data.frame(sc_dat_pre[,1:2], 
                                scaled_Brier = round((1-(sc_dat_pre$model/sc_dat_pre$null)), 3))
    
    sc_dattab <- pivot_wider(sc_dattab_pre, names_from = month, values_from = scaled_Brier, 
                             names_prefix = "Scaled Brier month ")
    
    fnl_dat <- list(crude_Br = fnl_dat_pre, scaled_Br = sc_dattab)
  } else{
    fnl_dat <- fnl_dat_pre
    }
  return(fnl_dat)
}

make_cal_intsl <- function(dattab, trtsel.eval = F){
  
  fnl_dat <- lapply(dattab, function(x){
    rnames <- rownames(x)
    out_1 <- apply(x, 1, function(y){
      y <- round(y, 3)
      out_2 <- paste0(y[1], " (", y[2], "-", y[3], ")")
      return(out_2)})
    out_1 <- data.frame(measure = rnames, estimate = out_1)
    if(trtsel.eval == T){
      out_1 <- data.frame(left_join(out_1, trtsel_conv_table, by = c("measure"="var_name")))
    }
    return(out_1)})
  

  
  return(fnl_dat)
}

make_coef_great <- function(dat){
  
  if("variable" %in% colnames(dat)){
    var_list <- dat$variable
    pred_namecat_use <-  pred_namecat
    
    for(p in 1:nrow(pred_namecat)){
      for(d in 1:nrow(dict)) {
        if(str_detect(pred_namecat_use$name[p], paste0("^", dict$pre[d], "_"))){
          pred_namecat_use$label[p] <- str_c(dict$domain[d], ":", pred_namecat_use$label[p])
        }
      }
      
      var_list <- str_replace(var_list, paste0(pred_namecat_use$name[p], "(?!_)"), 
                              pred_namecat_use$label[p])
    }
    
    var_list <- str_replace(var_list, "\\.{3}", "-")
    var_list <- str_replace(var_list, "Yes", "=Yes")
    dat$variable <- str_replace(var_list, "DrugFTY720.", "DrugFTY720*")
  }
    return(dat)
}
