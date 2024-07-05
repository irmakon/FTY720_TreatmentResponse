# to be copy-pasted within the following function:
# fixInNamespace("one.boot.plot", pos = "package:TreatmentSelection", editor = "internal")

function (x, ci, fixed.values, fix.ind, out.ind) 
{
  if (x$model.fit$outcome == "time-to-event") {
    mysurv <- with(x$derived.data, eval(x$formula[[2]]))
    event <- mysurv[, 2]
  }
  else {
    event <- x$derived.data[[as.character(x$formula[[2]])]]
  }
  myboot.sample <- x$functions$boot.sample(event, x$derived.data[[x$treatment.name]], 
                                           rho = x$model.fit$cohort.attributes)
  if (substr(ci, 1, 1) == "h") 
    addind = 0
  else addind = 1
  rho.b <- myboot.sample[1:7]
  ind <- myboot.sample[-c(1:7)]
  if (x$model.fit$outcome == "time-to-event") {
    event.b <- 0
  }
  else {
    event.b <- x$derived.data[[as.character(x$formula[[2]])]][ind]
  }
  trt.b <- x$derived.data[[x$treatment.name]][ind]
  if (x$model.fit$family$family == "risks_provided") {
    obsrisk.t0.b <- x$derived.data$fittedrisk.t0[ind]
    obsrisk.t1.b <- x$derived.data$fittedrisk.t1[ind]
    linkinvfun <- NULL
    marker.b <- obsrisk.t0.b - obsrisk.t1.b
  }
  else if (x$model.fit$family$family == "time-to-event") {
    coxfit <- do.call(coxph, list(x$formula, x$derived.data[ind, 
                                                            ]))
    obsrisk.t0.b <- get.risk.t_coxph(coxfit, x$treatment.name, 
                                     x$derived.data[ind, ], x$prediction.time, t = 0)
    obsrisk.t1.b <- get.risk.t_coxph(coxfit, x$treatment.name, 
                                     x$derived.data[ind, ], x$prediction.time, t = 1)
  }
  else {
    coef <- unname(get.coef(x$formula, x$treatment.name, 
                            x$derived.data[ind, ], x$model.fit$study.design, 
                            rho.b, family = x$model.fit$family)[, 1])
    linkinvfun <- x$model.fit$family$linkinv
    obsrisk.t0.b <- get.risk.t(coef, x$formula, x$treatment.name, 
                               data = x$derived.data[ind, ], linkinvfun, t = 0)
    obsrisk.t1.b <- get.risk.t(coef, x$formula, x$treatment.name, 
                               data = x$derived.data[ind, ], linkinvfun, t = 1)
    wi = 0
  }
  obsdelta.b <- obsrisk.t0.b - obsrisk.t1.b
  if (length(x$model.fit$marker.names) == 1) {
    marker.b <- x$derived.data[ind, x$model.fit$marker.names]
  }
  else {
    marker.b <- NULL
  }
  if (is.null(marker.b)) 
    marker.b <- obsrisk.t0.b - obsrisk.t1.b
  F.Y <- x$functions$get.F(marker.b, event.b, trt.b, rho.b) * 
    100
  F.D <- x$functions$get.F(obsdelta.b, event.b, trt.b, rho.b) * 
    100
  theta.c <- EventRateVec(obsrisk.t0.b, obsrisk.t1.b, F.D, 
                          rho.b, event.b, trt.b)
  all <- cbind(F.Y, obsrisk.t0.b, obsrisk.t1.b, F.D, obsdelta.b, 
               theta.c)
  all <- unique(all)
  if (length(fix.ind) > 1) {
    myorder <- apply(all[, fix.ind], 2, order)
    out <- matrix(0, ncol = length(fixed.values), nrow = length(fix.ind))
    for (i in 1:length(fix.ind)) {
      if (is.element(fix.ind[i], c(2, 3)) & !(all.equal(order(all[, 
                                                                  fix.ind[i]]), order(all[, 1])) == TRUE)) {
        addind = 1
      }
      tmpind <- sum.I(fixed.values, ">=", all[myorder[, 
                                                      i], fix.ind[i]]) + addind
      tmpind[tmpind <= 0] <- NA
      tmpall <- all[myorder[, i], out.ind[i]]
      out[i, ] <- tmpall[tmpind]
    }
  }
  else {
    myorder <- order(all[, fix.ind])
    out <- numeric(length(fixed.values))
    tmpind <- sum.I(fixed.values, ">=", all[myorder, fix.ind]) + 
      addind
    tmpind[tmpind == 0] <- NA
    tmpall <- all[myorder, out.ind]
    out <- tmpall[tmpind]
  }
  out
}
