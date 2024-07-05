# to be copy-pasted within the following function:
# fixInNamespace("evaluate.trtsel", pos = "package:TreatmentSelection", editor = "internal")

function (x, ..., bias.correct = TRUE, bootstraps = 1000, alpha = 0.05) 
{
  if (!is.trtsel(x)) 
    stop("x must be an object of class 'trtsel' created by using the function 'trtsel' see ?trtsel for more help")
  if (alpha < 0 | alpha > 1) 
    stop("Error: alpha should be between 0 and 1")
  if (bootstraps == 0) 
    print("bootstrap confidence intervals will not be calculated")
  if (bootstraps == 1) 
    warning("Number of bootstraps must be greater than 1, bootstrap confidence intervals will not be computed")
  stopifnot(is.logical(bias.correct))
  if (x$model.fit$family$family == "risks_provided") {
    bias.correct = FALSE
  }
  if (missing(bias.correct)) {
    if (x$model.fit$family$family == "risks_provided") {
      bias.correct = FALSE
    }
    else if (bootstraps > 1) {
      bias.correct = TRUE
      message("Bootstrap bias-correction will be implemented to correct for over-optimism bias in estimation.")
    }
    else {
      bias.correct = FALSE
    }
  }
  if (bias.correct & x$model.fit$family$family == "risks_provided") {
    message("risks are already provided from a fitted model, bias-correction is not available. \n Reported measures will not be bias-corrected.")
    bias.correct = FALSE
  }
  data <- x$derived.data
  study.design <- x$model.fit$study.design
  rho <- x$model.fit$cohort.attributes
  boot.sample <- x$functions$boot.sample
  get.summary.measures <- x$functions$get.summary.measures
  test.Null.val <- NA
  event.name = as.character(x$formula[[2]])
  treatment.name = x$treatment.name
  family <- x$model.fit$family
  if (x$model.fit$family$family == "risks_provided") 
    provided_risk <- cbind(x$derived.data$fittedrisk.t0, 
                           x$derived.data$fittedrisk.t1)
  else provided_risk = NULL
  data$prediction.time = x$prediction.time
  if (bootstraps > 1) {
    boot.data <- replicate(bootstraps, one.boot.eval(data = data, 
                                                     formula = x$formula, treatment.name = x$treatment.name, 
                                                     rho = rho, d = x$model.fit$thresh, study.design = study.design, 
                                                     obe.boot.sample = boot.sample, obe.get.summary.measures = get.summary.measures, 
                                                     family = family, disc.rec.no.trt = x$model.fit$disc.rec.no.trt, 
                                                     provided_risk = provided_risk, prediction.time = x$prediction.time, 
                                                     bbc = bias.correct))
    quantile <- NULL
    conf.intervals <- apply(boot.data[-c(1:4), ], 1, quantile, 
                            probs = c(alpha/2, 1 - alpha/2), type = 1, na.rm = TRUE)
    row.names(conf.intervals) <- c("lower", "upper")
    if (bias.correct) {
      bias <- apply(boot.data[c(5:36), ], 1, mean, na.rm = TRUE)
      bias <- bias[1:16] - bias[17:32]
      conf.intervals <- conf.intervals[, 1:16] - rbind(bias, 
                                                       bias)
    }
    else {
      bias = rep(0, dim(conf.intervals)[[2]]/2)
    }
  }
  else {
    conf.intervals = NULL
  }
  if (x$model.fit$outcome == "time-to-event") 
    event.name = x$formula[[2]]
  summary.measures <- data.frame(get.summary.measures(data, 
                                                      event.name, treatment.name, rho))
  if (any(data$rec.no.trt == 0) & any(data$rec.no.trt == 1) & 
      is.null(x$model.fit$disc.rec.no.trt) & family$family != 
      "risks_provided") {
    if (length(x$model.fit$marker.names) == 1) {
      marker = data[[x$model.fit$marker.names]]
      summary.measures$Marker.Thresh <- ifelse(data$trt.effect[which.min(marker)] < 
                                                 0, max(marker[data$rec.no.trt == 1]), min(marker[data$rec.no.trt == 
                                                                                                    1]))
    }
  }
  else if (any(data$rec.trt == 1) & any(data$rec.trt == 0) & 
           is.null(x$model.fit$disc.rec.no.trt) & family$family != 
           "risks_provided") {
    if (length(x$model.fit$marker.names) == 1) {
      marker = data[[x$model.fit$marker.names]]
      summary.measures$Marker.Thresh <- ifelse(data$trt.effect[which.min(marker)] < 
                                                 0, max(marker[data$rec.trt == 0]), min(marker[data$rec.trt == 
                                                                                                 0]))
    }
  }
  if (is.null(summary.measures$Marker.Thresh)) 
    summary.measures$Marker.Thresh <- NA
  result <- list(estimates = summary.measures, conf.intervals = conf.intervals, 
                 bias = bias, bias.correct = bias.correct)
  if (!is.null(x$model.fit$disc.rec.no.trt)) 
    result$discrete.marker = TRUE
  class(result) <- "eval.trtsel"
  return(result)
}
