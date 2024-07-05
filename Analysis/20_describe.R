################################################################################
# Describe the datasets
################################################################################

# study <- "N2309"
# source(paste0(getwd(), "/Code/00_paths_pckgs_functions.R"))

load(paste0(path_data,"convdat_", ana_type,".RData"))

### check the correspondence of two studies
# N2301_out <- outcome
# N2309_pred <- desc_predic_all
# 
# all(colnames(N2301_out)==colnames(outcome))
# all(colnames(N2301_pred)==colnames(predictors))

# N2309_desc_predic_all <- desc_predic_all

# summary(N2301_pred)
# 
# desc_predic_all$smry_cat$category
# 
# N2309_pred$smry_num$labels==desc_predic_all$smry_num$labels
# 
# test <- cbind(N2309_desc_predic_all$smry_cat$labels, desc_predic_all$smry_cat$labels)

# N2309_ana_dat <- ana_dat


# For manual reporting
nrow(subset(predictors, EDSS>=6))
nrow(subset(predictors, RlpsDist>24))

# Only for grouping the variables
pred_names <- colnames(predictors)[-1]
pred_cat <- c("Drug", rep("Demographic", 4), rep("Clinical", 8), 
              rep("MRI", 4), rep("Symptoms", 3), rep("MS drug history", 4), 
              rep("Clinical", 2), rep("QoL", 6), rep("Clinical", 3), 
              rep("Comedications", 11), rep("Concomitant diseases", 19), 
              rep("Laboratory", 16))


pred_cat_order <- data.frame(category = c("Drug", "Demographic", "Clinical", "Symptoms", 
                                          "MS drug history", "MRI", "QoL", "Comedications", 
                                          "Concomitant diseases", "Laboratory"), 
                        order = c(1:10))

dict <- data.frame(pre = c("CoMed", "CoDis", "QoL", "Lab", "EDSS", "MSFC"), 
                   domain = c("Comedication", "Concomitant Disease", 
                              "Quality of Life", "Lab", "EDSS", "MSFC"))

pred_labels <- NULL
for(i in 2:ncol(predictors)){
  pred_labels <-  c(pred_labels, attr(predictors[[i]], which = "label"))
}


pred_namecat_pre <- data.frame(name = pred_names, category = pred_cat, label = pred_labels)

pred_namecat <- inner_join(pred_namecat_pre, pred_cat_order) %>%
  arrange(order, name) %>%
  mutate(order = row_number()) %>%
  data.frame()

test2 <- lapply(tapply(pred_namecat$label, pred_namecat$category, paste), 
                function(x){n_pred <- length(x) 
                pred_list <- paste(x, collapse = "; ") 
                return(c(n_pred, pred_list))
  })

desc_predcat <- NULL
for(i in 1:length(test2)){
  desc_predcat <- rbind(desc_predcat, c(names(test2[i]), test2[[i]]))
}

desc_predcat_pre <- data.frame(desc_predcat)
colnames(desc_predcat_pre) <- c("Category", "Number", "List")

desc_predcat <- inner_join(desc_predcat_pre, pred_cat_order, by = c("Category"="category")) %>%
  arrange(order) %>%
  select(-order) %>%
  data.frame()


### OUTCOMES

p <- 2
# create labels for fine presentation
outcome_cat <- data.frame(pattern = unique(str_split_fixed(colnames(outcome[-c(1:p)]), n =2, pattern = "_")[,2]),
                          col_labels = c("Serious Adverse Event (SAE)", "Infections", "Neoplasms", 
                                     "Relapse", "Death", "Discontinue due to AE (Disc AE)", 
                                     "New/enlarging lesions (T2 MRI)",
                                     "Confirmed disability progression (3m CDP)", "Safety", 
                                     "Immunosuppresant safety (Immune safety)", 
                                     "Safety and efficacy (Composite)"),
                          stringsAsFactors = F)

outcome_cat$short_labels <- with(outcome_cat, ifelse(str_detect(col_labels, "\\(")!=T,
                                                     col_labels,
                                                     str_replace(str_split_fixed(col_labels, pattern="\\(", n=2)[,2],
                                                                 "\\)", "")))
c_lab <- NULL
c_attr <- NULL

for(i in 1:nrow(outcome_cat)){
  lab <- rep(outcome_cat$short_labels[i], 3)
  col_att <- rep(outcome_cat$col_labels[i], 3)
  
  names(lab) <- colnames(outcome)[p+1:3+(3*(i-1))]
  names(col_att) <- names(lab)
  
  c_lab <- c(c_lab, lab)
  c_attr <- c(c_attr, col_att)
}

outcome_fnl <- select(outcome, -ends_with(c("sae", "infec", "neopl", "death", "_ae")))
data_pre <- data.frame(outcome_fnl)

for(i in 1:((ncol(outcome_fnl)-p)/3)){
  suff <- str_split_fixed(colnames(outcome_fnl)[i*3], "_", 2)[2]
  surv <- Surv(time = data_pre[,i*3+2], event = data_pre[,i*3+1])
  outcome_fnl[, ncol(outcome_fnl)+1] <- surv
  colnames(outcome_fnl)[ncol(outcome_fnl)] <- paste0("surv_", suff)
}

# if long labels are needed, this one can just be switched to c_attr
# label(outcome) <- as.list(c_lab[match(names(outcome), names(c_lab))])

# extract

outcome_wdrug <- inner_join(outcome, select(predictors, Drug, STYSID1A),
                            by = c("stysid1a"="STYSID1A"))
bi <- data.frame(select (outcome_wdrug, starts_with("outcome_"), Drug))

outcome_wdrug_fnl <- inner_join(outcome_fnl, select(predictors, Drug, STYSID1A),
                                by = c("stysid1a"="STYSID1A"))
bi_fnl <- data.frame(select (outcome_wdrug_fnl, starts_with("outcome_"), Drug))
# cens <- data.frame(select(outcome, starts_with("censor_")))
# days <- data.frame(select(outcome, starts_with("day_")))

## Describe the follow-up
last_vis <- outcome$DayLastVisit

desc_fu <- desc.num(last_vis)

## Describe the binary
desc_bi <- data.frame(cbind(labels = outcome_cat$col_labels, 
                            t(apply(select(bi, -Drug), 2, desc.num))))

plot.desc.bi <- function(dat, scheme = "color"){
  
  if(scheme=="color"){
    colcodes <- other_colcod
  } else if (scheme=="gray"){
    colcodes <- other_graycod
  }
  
  
  dat_fig_desc_bi <- 
    pivot_longer(dat, starts_with("outcome_"), names_to = "type", names_prefix = "outcome_",
                                  values_to = "outcome") %>%
    mutate(outcome = factor(-outcome+2, labels = c("yes", "no")),
           type = factor(type, levels = outcome_cat$pattern[c(6,5,3,2,1,11,10,9,8,7,4)], 
                         labels = outcome_cat$short_labels[c(6,5,3,2,1,11,10,9,8,7,4)]))
  
  fig_desc_bi <- ggplot(dat_fig_desc_bi, aes(y = type, fill = outcome)) + 
    geom_bar(na.rm = T, position = position_stack(reverse = TRUE)) +
    coord_cartesian(expand = F) +
    facet_wrap(~Drug, ncol = 2, scales = "free_x") +
    labs(x = "count", y = "outcome type", 
         title = paste("Outcome by 24 months -", sname)) +
    scale_fill_manual(values = colcodes[c(6,3)], na.value = "black")
    theme_bw() +
    theme(text = element_text(family = "sans", size = 7))
  return(fig_desc_bi)
}

(fig_desc_bi_all <- plot.desc.bi(bi))
(fig_desc_bi_fnl <- plot.desc.bi(bi_fnl))

(fig_desc_bi_all_gray <- plot.desc.bi(bi, scheme = "gray"))
(fig_desc_bi_fnl_gray <- plot.desc.bi(bi_fnl, scheme = "gray"))


### PREDICTORS

desc_predic_all <- smry.per(predictors)

# For manual reporting
nomiss <- length(which(desc_predic_all$smry_num$n_miss==0)) + length(which(desc_predic_all$smry_cat$n_miss==0))
nomiss/ncol(predictors)


### JOINED DATASETS

ana_dat <- inner_join(outcome_fnl, predictors, by = c("stysid1a" = "STYSID1A"))

### Describe by Drug

ana_dat_desc <- select(ana_dat, -starts_with(c("surv_", "day_", "censor_")))

desc_by_drug <- smry.all(ana_dat_desc, "Drug")

### which are different between drug groups? (Did randomization work?) yes!
desc_by_drug$p_wilcox_chisq[which(desc_by_drug$p_wilcox_chisq<0.05)]

### Describe by Outcome

outcome_list <- as.list(colnames(select(ana_dat_desc, starts_with("outcome_"))))

desc_by_outcome <- lapply(outcome_list, function(x){
  smry.all(ana_dat_desc, x)
})

names(desc_by_outcome) <- outnames_convert(outcome_cat, 
                                           colnames(select(ana_dat_desc, starts_with("outcome_"))))

### which are different between outcomes? 
outcomediff <- lapply(desc_by_outcome, function(x){
  x$p_wilcox_chisq[which(x$p_wilcox_chisq<0.05)]
})

### Describe by Missing outcome

ana_dat_miss <- mutate_at(ana_dat_desc, vars(starts_with("outcome_")), ~ifelse(is.na(.x), "miss", "notmiss"))

desc_by_miss <- lapply(outcome_list, function(x){
  smry.all(ana_dat_miss, x)
})

names(desc_by_miss) <- outnames_convert(outcome_cat, 
                                        colnames(select(ana_dat_desc, starts_with("outcome_"))))

### which are different between missing vs not outcomes?
outcomemiss <- lapply(desc_by_miss, function(x){
  x$p_wilcox_chisq[which(x$p_wilcox_chisq<0.05)]
})

### Describe by missing outcome

formiss <- mutate(ana_dat_desc, censored = ifelse(DayLastVisit>675, 0, 1))
desc_by_censored <- smry.all(formiss, "censored")


### Describe missing

desc_miss_pred_pre <- desc.num(apply(predictors[,-1], 2, 
                                     function(x) 100*length(which(is.na(x)==T))/length(x)))

desc_miss_out_pre <- apply(bi[,-12], 2, 
                           function(x) 100*length(which(is.na(x)==T))/length(x))

desc_miss_pred <- data.frame(names = names(desc_miss_pred_pre), pct = desc_miss_pred_pre)
desc_miss_out <- data.frame(names = outnames_convert(outcome_cat, names(desc_miss_out_pre)), 
                            pct = desc_miss_out_pre)
desc_miss_pat <- length(which(apply(predictors[,-1], 1, 
                                    function(x) any(is.na(x)))==T))/nrow(predictors)


fig_miss_dat <- ana_dat_desc %>%
  select(-c(1:2)) %>%
  mutate(id = row_number()) %>%
  gather(-id, key = "key", value = "val") %>%
  mutate(isna = is.na(val)) %>%
  left_join(select(pred_namecat, name, order), by=c("key"="name")) 

for(lbl in 1:nrow(outcome_cat)){
  fig_miss_dat$key <- str_replace(fig_miss_dat$key, paste0("outcome_", outcome_cat$pattern[lbl]), 
                                  paste0("outcome_", outcome_cat$short_labels[lbl]))
}

# for(plb in 1:nrow(pred_namecat)){
#   fig_miss_dat$key <- str_replace(fig_miss_dat$key, pred_namecat$name[plb], 
#                                   paste0(pred_namecat$category[plb], " (",pred_namecat$name[plb], ")"))
# }

fig_miss_dat<- fig_miss_dat %>%
  mutate(key_n = ifelse(is.na(order), 1, order+1)) %>%
  arrange(key_n) %>%
  mutate(key = sub("_", " ", key),
         key = factor(key, levels = rev(unique(key)))) %>%
  select(-order, -key_n)
         
plot.desc.miss <- function(dat, scheme = "color"){
  
  if(scheme=="color"){
    colcodes <- other_colcod
    } else if (scheme=="gray"){
      colcodes <- other_graycod
      }
  
  fig <- dat %>%
    ggplot(aes(key, id, fill = isna)) +
    geom_raster(alpha = 0.8) +
    scale_fill_manual(name = "",
                      values = colcodes[c(7,2)],
                      labels = c("Present", "Missing")) +
    labs(x = "Variable",
         y = "Row Number",
         title = paste("Missing values -", sname)) +
    coord_flip(expand = F) +
    # coord_trans(ytrans = "-")
    theme_bw() +
    theme(text = element_text(family = "sans", size = 7))
  
  return(fig)
}

                   
(fig_miss <- plot.desc.miss(fig_miss_dat))
(fig_miss_gray <- plot.desc.miss(fig_miss_dat, scheme="gray"))

md.pattern(ana_dat_desc, rotate.names = T)

remove(ana_dat_desc, ana_dat_miss, bi, bi_fnl, data_pre, outcome, outcome_fnl, outcome_list, 
       predictors, outcome_wdrug, outcome_wdrug_fnl, fig_miss_dat)

## Describe the survival outcome
### ONCE THE TREATMENT ARRIVES, DO IT BY TREATMENT

ana_dat_surv <- select(ana_dat, starts_with("surv_"), Drug) %>%
  data.frame()

surv_list <- data.frame(names = colnames(select(ana_dat, starts_with("surv_"))))
surv_list$pattern <- str_split_fixed(surv_list$names, n =2, pattern = "_")[,2]
surv_list  <- left_join(surv_list, outcome_cat) 

dat_for_surv <- ana_dat_surv
dat_for_surv$Drug <- relevel(dat_for_surv$Drug, ref = "FTY720")

plot.desc.surv <- function(dat_for_surv, surv_list, scheme = "color"){
  
  if(scheme=="color"){
    colcodes <- drug_colcod
    } else if (scheme=="gray"){
      colcodes <- drug_graycod
      }

  fig_surv_list <- NULL
  
  for(i in 1:nrow(surv_list)){
    dat_fig_desc_surv <- survfit(dat_for_surv[,i]~Drug, data = dat_for_surv)
    fig_surv <- ggsurvplot(dat_fig_desc_surv,
                           risk.table = T,
                           # linetype = "strata",
                           conf.int = T,
                           title = paste(surv_list$col_labels[i], "-", sname),
                           xlim = c(0, maxday),
                           palette = colcodes,
                           tables.y.text = F,
                           tables.col = "strata",
                           risk.table.pos = "in",
                           axes.offset = F,
                           censor.shape = 124,
                           censor.size = 3,
                           # xlab = "Time",
                           break.x.by = 180,
                           ggtheme = theme(text = element_text(family = "sans", size = 7),
                                           legend.key.size = unit(0.2, "pt"),
                                           legend.title = element_blank(),
                                           # axis.title.x = element_blank(),
                                           plot.margin = unit(c(5, 5, 2, 5), "pt")),
                           tables.theme = theme(text = element_text(family = "sans", size = 7,
                                                                    face = "bold"),
                                                plot.margin = unit(c(0, 5, 5, 5), "pt"),
                                                # xlab = "Time"),
                                                axis.title.x = element_blank()),
                           fontsize = 2.2, font.family = "sans",
                           legend = c(0.16, 0.22))
    # text = element_text(family = "sans", size = 7))
    fig_surv_list[[i]] <- fig_surv
  }
  
  names(fig_surv_list) <- surv_list$pattern

  return(fig_surv_list)
  
  # return(dat_fig_desc_surv)
}

fig_surv <- plot.desc.surv(dat_for_surv, surv_list)
fig_surv_gray <- plot.desc.surv(dat_for_surv, surv_list, scheme = "gray")

ana_dat <- select(ana_dat, -starts_with("surv_"))

# save resulting datasets
save(ana_dat, file = paste0(path_data, "ana_dat_", ana_type,".RData"))

save(list=c("outcome_cat", "pred_namecat", "dict"),
     file = paste0(path_data, "dict_all_", ana_type,".RData"))

# save results to be reported
save(list = ls()[str_which(ls(), "^desc_")], 
     file = paste0(path_results, "desc_base_", ana_type,".RData"))
save(list = ls()[str_which(ls(), "^fig_")], 
     file = paste0(path_results, "desc_base_fig_", ana_type,".RData"))

