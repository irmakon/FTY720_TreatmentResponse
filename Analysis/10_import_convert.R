################################################################################
# Import the datasets provided by DM
################################################################################

# study <- "N2301"
# source(paste0(getwd(), "/Code/00_paths_pckgs_functions.R"))

### PREDICTORS

# LAB
labor <- read_sas(paste0(path_data, "labor", file_form))
labor_lbl <- read.csv(paste0(path_main, "/lab_conv_table.csv"))

str(labor)
summary(labor)
# Possible useless values: URPROTST (low variability), PT (missing), ALB (missing)

table(labor$Lab_UPROTST)

#Also eliminate PT because it does not exist in other studies

labor_con <- select(labor, -Lab_UPROTST)

if(study=="N2301"){
  labor_con <- select(labor_con, -Lab_PT)
}

for(nc in 2:ncol(labor_con)){
  cname <- str_sub(colnames(labor_con)[nc], start = 5)
  labval <- labor_lbl$Label[which(labor_lbl$Name==cname)]
  attr(labor_con[[nc]], which = "label") <- paste(labval,  attr(labor_con[[nc]], which = "label"))
}

# CODIS
codis <- read_sas(paste0(path_data, "codis", file_form))

str(codis)
summary(codis)
# Possible useless values: Preg, Hepat, Surg, Inj_P (low variability)

table(codis$CoDis_Preg)
table(codis$CoDis_Hepat)
table(codis$CoDis_Inj_P)
table(codis$CoDis_Surg)
table(codis$CoDis_SocCi)

codis_con <- mutate(codis, 
                    CoDis_Blood = case_when(CoDis_Hepat == 1 | CoDis_Inj_P == 1 
                                            | CoDis_Surg == 1 | CoDis_Preg == 1 
                                            | CoDis_SocCi == 1 | CoDis_Blood == 1 
                                            | CoDis_Card == 1 | CoDis_Ear == 1 ~ 1, 
                                            TRUE ~ 0)) %>%
  mutate_at(-1, factor, levels = c(0, 1), labels = c("No", "Yes")) %>%
  rename(CoDis_Other = CoDis_Blood)

for(i in 1:ncol(codis_con)){
  attr(codis_con[[i]], which="label") <- attr(codis[[i]], which="label")
}

attr(codis_con[[which(colnames(codis_con)=="CoDis_Other")]], 
     which="label") <- 
  paste(attr(codis_con[[which(colnames(codis_con)=="CoDis_Other")]], 
             which="label"), 
             "or", 
        attr(codis_con[[which(colnames(codis_con)=="CoDis_Card")]], 
             which="label"), 
        "or", 
        attr(codis_con[[which(colnames(codis_con)=="CoDis_Ear")]], 
             which="label"), 
        "or", 
        attr(codis_con[[which(colnames(codis_con)=="CoDis_Hepat")]], 
             which="label"), 
        "or", 
        attr(codis_con[[which(colnames(codis_con)=="CoDis_Inj_P")]], 
             which="label"), 
        "or", 
        attr(codis_con[[which(colnames(codis_con)=="CoDis_Preg")]], 
             which="label"), 
        "or", 
        attr(codis_con[[which(colnames(codis_con)=="CoDis_SocCi")]], 
             which="label"), 
        "or", 
        attr(codis_con[[which(colnames(codis_con)=="CoDis_Surg")]], 
             which="label"))

codis_con <- 
  select(codis_con, -CoDis_Hepat, -CoDis_Inj_P, -CoDis_Surg, 
         -CoDis_Preg, -CoDis_SocCi, -CoDis_Card, -CoDis_Ear)

summary(codis_con)

# COMED
comed <- read_sas(paste0(path_data, "comed", file_form)) 

str(comed)
summary(comed)
# Possible useless values: CoMed_P, CoMed_S, CoMed_L

table(comed$CoMed_P)
table(comed$CoMed_S)
table(comed$CoMed_L)

comed_con <- mutate(comed, CoMed_J = case_when(CoMed_S == 1 | CoMed_L == 1 
                                               | CoMed_P == 1 | CoMed_J == 1 ~ 1,
                                               TRUE ~ 0)) %>%
  mutate_at(-1, factor, levels = c(0, 1), labels = c("No", "Yes")) %>%
  rename(CoMed_Other = CoMed_J)

for(i in 1:ncol(comed_con)){
  attr(comed_con[[i]], which="label") <- attr(comed[[i]], which="label")
}

attr(comed_con[[which(colnames(comed_con)=="CoMed_Other")]], 
     which="label") <- 
  paste(attr(comed_con[[which(colnames(comed_con)=="CoMed_Other")]], 
             which="label"), 
        "or", 
        attr(comed_con[[which(colnames(comed_con)=="CoMed_L")]], 
             which="label"), 
        "or", 
        attr(comed_con[[which(colnames(comed_con)=="CoMed_P")]], 
             which="label"), 
        "or", 
        attr(comed_con[[which(colnames(comed_con)=="CoMed_S")]], 
             which="label"))

comed_con <- 
  select(comed_con, -CoMed_S, -CoMed_P, -CoMed_L)

summary(comed_con)

# atc_miss <- read_sas(paste0(path_data, "atccodemissing", file_form))
# nothing to be done for these 12 observations. Ignored

# PREDIC
predic <- read_sas(paste0(path_data, "predic", file_form))

str(predic)
summary(predic)
# Still needs some work in terms of factor codings


#check
# table(predic$VIS1N)

if(study=="N2309"){
  predic <- rename(predic,
                   TGPDSC1A = tgpdsc1a)
}

with(predic, table(TGPDSC1A, Drug))
table(predic$Age)

predic_con <- mutate(predic, 
                     Drug = factor(4-Drug, labels = c("Placebo", "FTY720 0.5mg", "FTY720 1.25mg")), 
                     Race = as.factor(ifelse(Race == 1, "Caucasian", 
                                             "non-Caucasian")), 
                     Sex = factor(Sex, labels = c("Male", "Female")), 
                     PriorMSDr = factor(PriorMSDr, levels = c(0, 1), labels = c("No", "Yes")),
                     Age = factor(Age, levels = c("16 - 20", "21 - 25", "26 - 30", "31 - 35", "36 - 40",
                                                      "41 - 45", "46 - 50", "51 - 55")),
                     PriorMSDr_IB = factor(PriorMSDr_IB, levels = c(0, 1), labels = c("No", "Yes")),
                     PriorMSDr_GA = factor(PriorMSDr_GA, levels = c(0, 1), labels = c("No", "Yes")),
                     PriorMSDr_Nat = factor(case_when(PriorMSDr_Nat == 1 | 
                                                              PriorMSDr_Other == 1 ~ 1,
                                                            TRUE ~ 0), 
                                                  levels = c(0, 1), labels = c("No", "Yes"))) %>%
  rename(PriorMSDr_Nat_Other = PriorMSDr_Nat) %>%
  select(-TRTREG1C, -TGPDSC1A, -PriorMSDr_Other, -PriorMSDr) %>%
  select(STYSID1A, Drug, everything())

summary(predic_con)

varname <- c("EDSS", "MRI_NumGdT1", 
             "MRI_VolGdT1", "MRI_VolT2", 
             "MRI_VolHypT1", "RlpsDist", 
             "QoL_Vas", "MSFC_25FW", "MSFC_9HPT", 
             "MSFC_PASAT", "PriorMSDr_IB", 
             "PriorMSDr_GA", "PriorMSDr_Nat_Other",  
             "PriorMSDr_Count", "Drug", "Sex", "Age", "Race", 
             "BMI", "STYSID1A")
labelname <- c("EDSS score (total)", "Number of Gd-enhanced T1 lesions",
               "Total volume of Gd-enhanced T1 lesions", "Total volume of T2 lesions", 
               "Total volume of T1 hypointense lesions", "Number of months since recent relapse",
               "Visual anlaog scale", "Mean of timed 25-foot walk", "Mean of 9-hole peg test", 
               "Paced auditory serial addition test", "Prior Interferon beta use",
               "Prior Glatiramer acetate use", "Prior Natalizumab or other MS treatment use", 
               "Number of prior MS treatments", "Drug", "Sex", "Age", "Race", 
               "Body Mass Index (kg/m^2)", "Study and subject ID")

for(cnt in 1:length(varname)){
  attr(predic_con[[which(colnames(predic_con)==varname[cnt])]], 
       which="label") <- labelname[cnt]
}




## RlpsDist problematic: there are 4 patients > 24 months (contrary to exclusion criteria)

## OUTCOMES

outcome <- read_sas(paste0(path_data, "outcome", file_form))

# check edss & ae (row) tables for decisions

#finalized definition of edss by looking into it
# edss <- read_sas(paste0(path_data, "edss", file_form))
#filtering !is.na(SOC_TXT) is harmless because they don't contain any SOC_ABB, too
# ae <- read_sas(paste0(path_data, "a_aeinf", file_form))

str(outcome)

summary(outcome)
outcome_fact <- data.frame(sapply(data.frame(outcome), as.factor))
summary(outcome_fact)

#run regular data quality checks
check_day <- select(outcome, starts_with("day_", ignore.case = F))
all(check_day<=765 & check_day>=1)
summary(check_day)

check_cens <- select(outcome_fact, starts_with("censor_", ignore.case = F),
                     starts_with("outcome_", ignore.case = F))

summary(check_cens)


sub_outcome <- select(outcome, starts_with(c("outcome_", "censor_", "day_"), ignore.case = F))
ind <- which(colnames(outcome)=="outcome_sae")

# for(i in 0:11){
#   bi <- data.frame(outcome)[, ind+(3*i)]
#   cens <- data.frame(outcome)[, ind+1+(3*i)]
#   days <- data.frame(outcome)[, ind+2+(3*i)]
#   
#   print(colnames(outcome[ind+(3*i)]))
#   
#   if(all(cens[which(bi==0 | is.na(bi))]==0) & all(cens[which(bi==0 | is.na(bi))]==0))
#     print("cens-outcome match")
#   else
#     print("cens-outcome DO NOT match")
#   
#   if(all(days[which(bi==0)]<=765 & days[which(bi==0)]>=676) & all(days[which(is.na(bi))]<676))
#     print("days-bi match")
#   else
#     print("days-bi DO NOT match")
# }

for(i in 1:n_outcome){
  bi <- data.frame(sub_outcome)[, i]
  cens <- data.frame(sub_outcome)[, i+n_outcome]
  days <- data.frame(sub_outcome)[, i+2*n_outcome]
  
  print(colnames(sub_outcome[i]))
  
  if(all(cens[which(bi==0 | is.na(bi))]==0) & all(cens[which(bi==0 | is.na(bi))]==0))
    print("cens-outcome match")
  else
    print("cens-outcome DO NOT match")
  
  if(all(days[which(bi==0)]<=maxday & days[which(bi==0)]>=minday) & all(days[which(is.na(bi))]<minday))
    print("days-bi match")
  else
    print("days-bi DO NOT match")
}


tapply(outcome$day_edss, outcome$outcome_edss, summary)

# tapply(outcome$TCDPG3N1, outcome$FCDPG3N1, summary)


###EDSS check

# edssvis_daycheck <- select(outcome, DayLastVisit, DayLastEdss) %>%
#   mutate(diff = DayLastVisit-DayLastEdss) %>%
#   filter(DayLastVisit<=765)
# 
# summary(edssvis_daycheck)
# 
# hist(edssvis_daycheck$diff)
# boxplot(edssvis_daycheck$diff)
# plot(edssvis_daycheck$DayLastVisit, edssvis_daycheck$DayLastEdss)
# hist(edssvis_daycheck$DayLastEdss)

# lastvisday <- select(outcome, STYSID1A, DayLastVisit, censor_edss, TCDPG3N1, day_edss) %>% 
#   subset(censor_edss==0 & TCDPG3N1<=765 & TCDPG3N1!=day_edss) %>%
#   mutate(dif = day_edss-TCDPG3N1) %>%
#   subset(dif<0)
# 
# summary(lastvisday)

# write.csv(lastvisday, file = paste0(path_results, "lastvisday.csv"))

predictors <- inner_join(predic_con, comed_con) %>%
  inner_join(codis_con) %>% 
  inner_join(labor_con)

labname <- NULL
for(cnt in 1:ncol(predictors)){
  labname[cnt] <- attr(predictors[[cnt]], 
                       which="label")
}


if(ana_type == "main"){
  predictors <- mutate(predictors, 
                       Drug = fct_collapse(Drug, FTY720 = c("FTY720 1.25mg", "FTY720 0.5mg")))
} else if (ana_type == "s_drug05"){
  ids_125 <- predictors$STYSID1A[which(predictors$Drug == "FTY720 1.25mg")]
  predictors <- mutate(predictors, 
                       Drug = fct_collapse(Drug, FTY720 = c("FTY720 1.25mg", "FTY720 0.5mg"))) %>% 
    subset(!(STYSID1A %in% ids_125))
  outcome <- subset(outcome, !(stysid1a %in% ids_125))
}

save(outcome, predictors, file = paste0(path_data, "convdat_", ana_type,".RData"))


# outcome_str <- outcome[1,]
# outcome_str[1,] <- NA
# write.csv(outcome_str, file = paste0(path_results, "outcome_str.csvcsv"))