################################################################################
# Study parameters - to be sourced
################################################################################

###DEFINE INDIVIDUALLY FOR EACH STUDY


#to differentiate between main analyses versus sensitivity
# ana_type <- "main"
ana_type <- "s_drug05"

# number of CV folds
n_folds=5

#strategy for lasso
glm_lambda <- "lambda.min"
# glm_lambda <- "lambda.1se"

if(study=="N2301"){
  n_outcome <- 11
  minday <- 676
  maxday <- 765
  
  visits <- c(180, 360, 720)
  
  sname <- "development"
  
  #number of outcomes
  n_out <- 6

} else if(study=="N2309"){
  n_outcome <- 11
  minday <- 676
  maxday <- 765
  
  visits <- c(180, 360, 720)
  
  sname <- "validation"
  
  #number of outcomes
  n_out <- 6
}

# if (ana_type == "s_drug05"){
#   n_folds=3
# }

