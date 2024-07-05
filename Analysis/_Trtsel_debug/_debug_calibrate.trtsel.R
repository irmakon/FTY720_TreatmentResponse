# to be copy-pasted within the following function:
# fixInNamespace("calibrate.trtsel", pos = "package:TreatmentSelection", editor = "internal")

function (x, ..., groups = 10, plot.type = "calibration", trt.names = c("Treatment", 
                                                                        "No Treatment"), line.color = "black", point.color = "grey10", 
          ci = 0.95, main = NULL, ylim = NULL, xlim = NULL, ylab = NULL, 
          xlab = NULL) 
{
  if (!is.trtsel(x)) 
    stop("x must be an object of class 'trtsel' created by using the function 'trtsel' see ?trtsel for more help")
  if (!is.null(x$model.fit$disc.rec.no.trt)) 
    stop("Calibration not supported for a discrete marker")
  event.name = as.character(x$formula[[2]])
  lower <- upper <- NULL
  if (x$model.fit$outcome == "time-to-event") {
    event.name = x$formula[[2]]
    mysurv <- with(x$derived.data, eval(event.name))
    event <- mysurv[, 2]
    stime <- mysurv[, 1]
  }
  else {
    event.name <- as.character(x$formula[[2]])
    event <- x$derived.data[[event.name]]
  }
  trt <- x$derived.data[[x$treatment.name]]
  n <- length(trt)
  if (!is.numeric(groups)) 
    stop("groups must be an integer")
  if (groups < 2) 
    stop("Must have more than 2 groups!")
  if (!is.element(plot.type, c("calibration", "risk.t0", "risk.t1", 
                               "treatment effect", NA, "none"))) {
    stop("available plot.type options are \"calibration\", \"risk.t0\", \"risk.t1\", \"treatment effect\", \"none\" or NA")
  }
  fittedrisk.t0 <- x$derived.data$fittedrisk.t0
  fittedrisk.t1 <- x$derived.data$fittedrisk.t1
  fitteddelta <- x$derived.data$trt.effect
  fittedrisk.c.t0 <- x$derived.data$fittedrisk.t0[trt == 0]
  fittedrisk.c.t1 <- x$derived.data$fittedrisk.t1[trt == 1]
  D.t0 <- event[trt == 0]
  D.t1 <- event[trt == 1]
  trt.t0 <- trt[trt == 0]
  n.t0 <- sum(trt == 0)
  trt.t1 <- trt[trt == 1]
  n.t1 <- sum(trt == 1)
  rho <- x$model.fit$cohort.attributes
  study.design <- x$model.fit$study.design
  if (study.design == "RCT") {
    F.risk.t0 <- get.F.cohort(fittedrisk.c.t0, D.t0, trt.t0, 
                              rho, return.fun = TRUE)
    F.risk.t1 <- get.F.cohort(fittedrisk.c.t1, D.t1, trt.t1, 
                              rho, return.fun = TRUE)
    F.delta <- get.F.cohort(fitteddelta, event, trt, rho, 
                            return.fun = TRUE)
  }
  else if (substr(study.design, 1, 4) == "nest") {
    F.risk.t0 <- get.F.case.control(fittedrisk.c.t0, D.t0, 
                                    trt.t0, rho, return.fun = TRUE)
    F.risk.t1 <- get.F.case.control(fittedrisk.c.t1, D.t1, 
                                    trt.t1, rho, return.fun = TRUE)
    F.delta <- get.F.case.control(fitteddelta, event, trt, 
                                  rho, return.fun = TRUE)
  }
  else if (substr(study.design, 1, 5) == "strat") {
    get.F.tmp <- function(marker, event, trt, rho, tmp.trt) {
      tmp.index <- ifelse(tmp.trt == 0, 1, 2)
      if (tmp.trt == 0) {
        Pr.D1.trt <- (rho[3])/(rho[3] + rho[2])
        Pr.D0.trt <- 1 - Pr.D1.trt
      }
      else if (tmp.trt == 1) {
        Pr.D1.trt <- (rho[5])/(rho[5] + rho[4])
        Pr.D0.trt <- 1 - Pr.D1.trt
      }
      marker.D1.trt <- marker[event == 1 & trt == tmp.trt]
      marker.D0.trt <- marker[event == 0 & trt == tmp.trt]
      FY.D1.trt <- ecdf(marker.D1.trt)
      FY.D0.trt <- ecdf(marker.D0.trt)
      function(x) FY.D1.trt(x) * (Pr.D1.trt) + FY.D0.trt(x) * 
        (Pr.D0.trt)
    }
    F.risk.t0 <- get.F.tmp(fittedrisk.c.t0, D.t0, trt.t0, 
                           rho, tmp.trt = 0)
    F.risk.t1 <- get.F.tmp(fittedrisk.c.t1, D.t1, trt.t1, 
                           rho, tmp.trt = 1)
    F.delta <- get.F.stratified.case.control(fitteddelta, 
                                             event, trt, rho, return.fun = TRUE)
  }
  else {
    stop("study.design not specified correctly")
  }
  breaks.t0 <- sort(fittedrisk.c.t0)[sum.I(seq(0, 1, 1/groups), 
                                           "<", F.risk.t0(fittedrisk.c.t0))]
  breaks.t1 <- sort(fittedrisk.c.t1)[sum.I(seq(0, 1, 1/groups), 
                                           "<", F.risk.t1(fittedrisk.c.t1))]
  breaks.delta <- sort(fitteddelta)[sum.I(seq(0, 1, 1/groups), 
                                          "<", F.delta(fitteddelta))]
  if (!(length(unique(round(breaks.t0, 5))) == (groups + 1) & 
        length(unique(round(breaks.t1, 5))) == (groups + 1) & 
        length(unique(round(breaks.delta, 5))) == (groups + 1))) {
    stop("Error: Too many groups, cut points are not unique. Please reduce number of groups")
  }
  cut.t0 <- cut(fittedrisk.c.t0, breaks = breaks.t0, include.lowest = TRUE)
  cut.t1 <- cut(fittedrisk.c.t1, breaks = breaks.t1, include.lowest = TRUE)
  cut.delta <- cut(fitteddelta, breaks = breaks.delta, include.lowest = TRUE)
  if (x$model.fit$outcome == "time-to-event") {
    sfit.t0 <- summary(survfit(Surv(stime[trt == 0], event[trt == 
                                                             0] == 1) ~ cut.t0, se.fit = TRUE), extend = TRUE, 
                       times = x$prediction.time)
    sfit.t1 <- summary(survfit(Surv(stime[trt == 1], event[trt == 
                                                             1] == 1) ~ cut.t1, se.fit = TRUE), extend = TRUE, 
                       times = x$prediction.time)
    obs.risk.t0 <- 1 - sfit.t0$surv
    obs.risk.t1 <- 1 - sfit.t1$surv
    obs.sd.t0 <- sfit.t0$std.err
    obs.sd.t1 <- sfit.t1$std.err
    sfit.t0.delta <- summary(survfit(Surv(stime[trt == 0], 
                                          event[trt == 0] == 1) ~ cut.delta[trt == 0], se.fit = TRUE), 
                             times = x$prediction.time, extend = TRUE)
    sfit.t1.delta <- summary(survfit(Surv(stime[trt == 1], 
                                          event[trt == 1] == 1) ~ cut.delta[trt == 1], se.fit = TRUE), 
                             times = x$prediction.time, extend = TRUE)
    obs.risk.t0.tmp <- 1 - sfit.t0.delta$surv
    obs.risk.t1.tmp <- 1 - sfit.t1.delta$surv
  }
  else {
    obs.risk.t0 <- aggregate(D.t0, by = list(cut.t0), FUN = "mean")$x
    obs.risk.t1 <- aggregate(D.t1, by = list(cut.t1), FUN = "mean")$x
    obs.sd.t0 <- aggregate(D.t0, by = list(cut.t0), FUN = "sd")$x
    obs.sd.t1 <- aggregate(D.t1, by = list(cut.t1), FUN = "sd")$x
    obs.risk.t1.tmp <- aggregate(D.t1, by = list(cut.delta[trt == 
                                                             1]), FUN = "mean")$x
    obs.risk.t0.tmp <- aggregate(D.t0, by = list(cut.delta[trt == 
                                                             0]), FUN = "mean")$x
    obs.var.t1.tmp <- aggregate(D.t1, by = list(cut.delta[trt == 
                                                            1]), FUN = "var")$x
    obs.var.t0.tmp <- aggregate(D.t0, by = list(cut.delta[trt == 
                                                            0]), FUN = "var")$x
  }
  exp.risk.t0 <- aggregate(fittedrisk.c.t0, by = list(cut.t0), 
                           FUN = "mean")$x
  ng.t0 <- as.numeric(unlist(table(cut.t0)))
  if (any(ng.t0 < 5)) 
    warning(paste("For Treatment = 0,", sum(ng.t0 < 5), "groups have less than 5 observations."))
  exp.risk.t1 <- aggregate(fittedrisk.c.t1, by = list(cut.t1), 
                           FUN = "mean")$x
  ng.t1 <- as.numeric(unlist(table(cut.t1)))
  if (any(ng.t1 < 5)) 
    warning(paste("For Treatment = 1,", sum(ng.t1 < 5), "groups have less than 5 observations."))
  exp.delta <- aggregate(fitteddelta, by = list(cut.delta), 
                         FUN = "mean")$x
  if (!(length(obs.risk.t0.tmp) == length(obs.risk.t1.tmp))) 
    stop("Failure to observe at least one individual per treatment arm in each group. Please reduce the number of groups")
  obs.delta <- obs.risk.t0.tmp - obs.risk.t1.tmp
  ng.delta <- as.numeric(unlist(table(cut.delta)))
  if (any(ng.delta < 5)) 
    warning(paste("For observed treatment effects,", sum(ng.delta < 
                                                           5), "groups have less than 5 observations."))
  if (substr(study.design, 1, 4) == "nest") {
    obs.risk.t0 = expit(logit(obs.risk.t0) + logit(rho[3]) - 
                          logit(mean(event)))
    obs.risk.t1 = expit(logit(obs.risk.t1) + logit(rho[3]) - 
                          logit(mean(event)))
    obs.delta = expit(logit(obs.risk.t1.tmp) + logit(rho[3]) - 
                        logit(mean(event))) - expit(logit(obs.risk.t0.tmp) + 
                                                      logit(rho[3]) - logit(mean(event)))
  }
  else if (substr(study.design, 1, 5) == "strat") {
    Pr.D1.givT0 <- rho[3]/(rho[2] + rho[3])
    Pr.D1.givT1 <- rho[5]/(rho[4] + rho[5])
    obs.risk.t0 = expit(logit(obs.risk.t0) + logit(Pr.D1.givT0) - 
                          logit(mean(event[trt == 0])))
    obs.risk.t1 = expit(logit(obs.risk.t1) + logit(Pr.D1.givT1) - 
                          logit(mean(event[trt == 1])))
    obs.delta = expit(logit(obs.risk.t1.tmp) + logit(Pr.D1.givT1) - 
                        logit(mean(event[trt == 1]))) - expit(logit(obs.risk.t0.tmp) + 
                                                                logit(Pr.D1.givT0) - logit(mean(event[trt == 0])))
  }
  if (study.design == "RCT") {
    hl.t0 <- sum(ng.t0 * (obs.risk.t0 - exp.risk.t0)^2/(exp.risk.t0 * 
                                                          (1 - exp.risk.t0)))
    hl.t1 <- sum(ng.t1 * (obs.risk.t1 - exp.risk.t1)^2/(exp.risk.t1 * 
                                                          (1 - exp.risk.t1)))
  }
  else {
    if (x$model.fit$family$family == "risks_provided") 
      stop("cannot calculate Hosmer Lemeshow statistic when fitted risks are provided and study design is not RCT")
    risk.naive.all <- fitted(glm(x$formula, data = x$derived.data, 
                                 family = x$model.fit$family))
    risk.naive.t0.all <- risk.naive.all[trt == 0]
    risk.naive.t1.all <- risk.naive.all[trt == 1]
    exp.risk.naive.t0 <- aggregate(risk.naive.t0.all, by = list(cut.t0), 
                                   FUN = "mean")$x
    exp.risk.naive.t1 <- aggregate(risk.naive.t1.all, by = list(cut.t1), 
                                   FUN = "mean")$x
    hl.t0 <- sum(ng.t0 * (obs.risk.t0 - exp.risk.t0)^2/((exp.risk.t0^2 * 
                                                           (1 - exp.risk.t0)^2)/(exp.risk.naive.t0 * (1 - exp.risk.naive.t0))))
    hl.t1 <- sum(ng.t1 * (obs.risk.t1 - exp.risk.t1)^2/((exp.risk.t1^2 * 
                                                           (1 - exp.risk.t1)^2)/(exp.risk.naive.t1 * (1 - exp.risk.naive.t1))))
  }
  Df <- groups - 2
  pval.t0 <- 1 - pchisq(hl.t0, Df)
  pval.t1 <- 1 - pchisq(hl.t1, Df)
  if (is.element(plot.type, c("calibration", "risk.t0", "risk.t1", 
                              "treatment effect"))) {
    min.risk <- min(c(fittedrisk.c.t0, fittedrisk.c.t1))
    max.risk <- max(c(fittedrisk.c.t0, fittedrisk.c.t1))
    cen <- mean(c(min.risk, max.risk))
    ran <- max.risk - min.risk
    ran <- ran * 1.1
    mylim <- c(cen - ran/2, cen + ran/2)
  }
  observedRisk <- expectedRisk <- F.risk <- risk <- y <- NULL
  if (is.element(plot.type, "calibration")) {
    if (!is.null(xlim)) {
      if (any(xlim < 0)) 
        stop("Parameters of xlim must be > 0 due to log scaling")
    }
    if (!is.null(ylim)) {
      if (any(ylim < 0)) 
        stop("Parameters of ylim must be > 0 due to log scaling")
    }
    if (is.null(xlab)) 
      xlab <- "Observed risk"
    if (is.null(ylab)) 
      ylab <- "Predicted risk"
    if (is.null(main)) 
      main <- "Calibration plot"
    data <- data.frame(observedRisk = c(obs.risk.t0, obs.risk.t1), 
                       expectedRisk = c(exp.risk.t0, exp.risk.t1), trt = rep(c(0, 
                                                                               1), c(length(obs.risk.t0), length(obs.risk.t1))))
    data <- subset(data, observedRisk > 0)
    data <- subset(data, expectedRisk > 0)
    p <- ggplot(data = data, aes(x = observedRisk, y = expectedRisk, 
                                 color = factor(trt, levels = c(1, 0), labels = trt.names))) +
      geom_point(size = 2) + 
      scale_color_manual(values = point.color, labels = trt.names)
    p <- p + # coord_trans(x = "log", y = "log") + 
      # scale_shape_discrete("", labels = trt.names) + 
      ylab(ylab) + xlab(xlab) + 
      ggtitle(main) + 
      geom_abline(intercept = 0, slope = 1, 
                  colour = "black", linetype = 3) +
      geom_smooth(aes(x = observedRisk, y = expectedRisk), 
                colour = "grey50", size = 1.1, se = F) + 
      coord_cartesian(expand = F) + 
      theme(legend.title = element_blank(), 
            legend.background = element_blank(), 
            legend.key = element_rect(fill = NA),
            legend.key.size = unit(1.5, "lines"), 
            text = element_text(size = 7),
            legend.position = c(0.98, 0.02),
            legend.justification = c(0.98, 0.02),
            legend.direction = "horizontal",
            legend.box = "horizontal",
            legend.margin = margin(0.1, 0.1, 0.1, 0.1, unit = "pt"))

    if (!is.null(xlim)) {
      # if (xlim[1] == 0) {
      #   warning("due to log scaling, the lower bound of xlim must be > 0, changing xlim[1] <- .01")
      #   xlim[1] <- 0.01
      # }
      p <- p + scale_x_continuous(limits = xlim)
    }
    if (!is.null(ylim)) {
      # if (ylim[1] == 0) {
      #   warning("due to log scaling, the lower bound of ylim must be > 0, changing ylim[1] <- .01")
      #   ylim[1] <- 0.01
      # }
      p <- p + scale_y_continuous(limits = ylim)
    }
    print(p)
  }
  if (is.element(plot.type, "risk.t0")) {
    if (is.null(xlab)) 
      xlab <- "% population below risk"
    if (is.null(ylab)) 
      ylab <- "Risk"
    if (is.null(xlim)) 
      xlim <- c(0, 100)
    if (is.null(main)) 
      main <- "Risk curve for non treated individuals"
    data = data.frame(F.risk = F.risk.t0(sort(fittedrisk.c.t0)) * 
                        100, risk = sort(fittedrisk.c.t0))
    p <- ggplot() + geom_step(data = data, direction = "vh", size = 1.1, 
                              color = line.color, aes(x = F.risk, y = risk))
    obsdata <- data.frame(x = (1:groups/groups - 1/(2 * groups)) * 
                            100, y = obs.risk.t0, sd = obs.sd.t0)
    if (x$model.fit$outcome == "time-to-event") {
      obsdata$upper <- 1 - sfit.t0$lower
      obsdata$lower <- 1 - sfit.t0$upper
    }
    else if (x$model.fit$outcome == "binary") {
      obsdata$upper <- binom.confint(obsdata$y * ng.t0, 
                                     ng.t0, methods = "wilson", conf.level = ci)$upper
      obsdata$lower <- binom.confint(obsdata$y * ng.t0, 
                                     ng.t0, methods = "wilson", conf.level = ci)$lower
    }
    else if (x$model.fit$outcome == "continuous") {
      obsdata$upper = obsdata$y + qnorm(1 - (1 - ci)/2) * 
        obsdata$sd/sqrt(ng.t0)
      obsdata$lower = obsdata$y - qnorm(1 - (1 - ci)/2) * 
        obsdata$sd/sqrt(ng.t0)
    }
    p <- p + geom_errorbar(data = obsdata, color = point.color, 
                           aes(ymin = lower, ymax = upper, x = x), width = 2) + 
      geom_point(data = obsdata, aes(x = x, y = y), color = point.color, size = 2) 
    p <- p + ylab(ylab) + xlab(xlab) + ggtitle(main) + 
      coord_cartesian(expand = F) +
      theme(text = element_text(size = 7))
    if (!is.null(xlim)) 
      p <- p + xlim(xlim)
    if (!is.null(ylim)) 
      p <- p + ylim(ylim)
    print(p)
  }
  if (is.element(plot.type, "risk.t1")) {
    if (is.null(xlab)) 
      xlab <- "% population below risk"
    if (is.null(ylab)) 
      ylab <- "Risk"
    if (is.null(xlim)) 
      xlim <- c(0, 100)
    if (is.null(main)) 
      main <- "Risk curve for treated individuals"
    data = data.frame(F.risk = F.risk.t1(sort(fittedrisk.c.t1)) * 
                        100, risk = sort(fittedrisk.c.t1))
    p <- ggplot() + geom_step(data = data, direction = "vh", size = 1.1, 
                              color = line.color, aes(x = F.risk, y = risk))
    obsdata <- data.frame(x = (1:groups/groups - 1/(2 * groups)) * 
                            100, y = obs.risk.t1, sd = obs.sd.t1)
    if (x$model.fit$outcome == "time-to-event") {
      obsdata$upper <- 1 - sfit.t1$lower
      obsdata$lower <- 1 - sfit.t1$upper
    }
    else if (x$model.fit$outcome == "binary") {
      obsdata$upper <- binom.confint(obsdata$y * ng.t1, 
                                     ng.t1, methods = "wilson")$upper
      obsdata$lower <- binom.confint(obsdata$y * ng.t1, 
                                     ng.t1, methods = "wilson")$lower
    }
    else if (x$model.fit$outcome == "continuous") {
      obsdata$upper = obsdata$y + qnorm(1 - (1 - ci)/2) * 
        obsdata$sd/sqrt(ng.t0)
      obsdata$lower = obsdata$y - qnorm(1 - (1 - ci)/2) * 
        obsdata$sd/sqrt(ng.t0)
    }
    p <- p + geom_errorbar(data = obsdata, color = point.color, 
                           aes(ymin = lower, ymax = upper, x = x), width = 2) + 
      geom_point(data = obsdata, aes(x = x, y = y), color = point.color, size = 2)
    p <- p + ylab(ylab) + xlab(xlab) + ggtitle(main) + 
      coord_cartesian(expand = F) +
      theme(text = element_text(size = 7))
    if (!is.null(xlim)) 
      p <- p + xlim(xlim)
    if (!is.null(ylim)) 
      p <- p + ylim(ylim)
    print(p)
  }
  if (is.element("treatment effect", plot.type)) {
    if (is.null(xlab)) 
      xlab <- "% population below treatment effect"
    if (is.null(ylab)) 
      ylab <- "Treatment effect"
    if (is.null(xlim)) 
      xlim <- c(0, 100)
    if (is.null(main)) 
      main <- "Treatment effect distribution"
    data = data.frame(F.risk = F.delta(sort(fitteddelta)) * 
                        100, risk = sort(fitteddelta))
    p <- p <- ggplot() + geom_step(data = data, direction = "vh", 
                                   color = line.color, aes(x = F.risk, y = risk))
    obsdata <- data.frame(x = (1:groups/groups - 1/(2 * groups)) * 
                            100, y = obs.delta)
    if (x$model.fit$outcome == "time-to-event") {
      obsdata$var <- sfit.t0.delta$std.err^2 + sfit.t1.delta$std.err^2
    }
    else if (x$model.fit$outcome == "binary") {
      obsdata$var <- obs.risk.t1.tmp * (1 - obs.risk.t1.tmp)/ng.t1 + 
        obs.risk.t0.tmp * (1 - obs.risk.t0.tmp)/ng.t0
    }
    else if (x$model.fit$outcome == "continuous") {
      obsdata$var <- obs.var.t1.tmp/ng.t1 + obs.var.t0.tmp/ng.t0
    }
    obsdata$upper <- obsdata$y + qnorm(1 - (1 - ci)/2) * 
      sqrt(obsdata$var)
    obsdata$lower <- obsdata$y + qnorm((1 - ci)/2) * sqrt(obsdata$var)
    p <- p + geom_errorbar(data = obsdata, color = point.color, 
                           aes(ymin = lower, ymax = upper, x = x), width = 2) + 
      geom_point(data = obsdata, aes(x = x, y = y), color = point.color, size = 1.2)
    p <- p + ylab(ylab) + xlab(xlab) + ggtitle(main)
    if (!is.null(xlim)) 
      p <- p + xlim(xlim)
    if (!is.null(ylim)) 
      p <- p + ylim(ylim)
    print(p)
  }
  if (!is.element(plot.type, c("calibration", "risk.t0", "risk.t1", 
                               "treatment effect"))) 
    p = NULL
  res <- list(HL.TestStat = c(trt0 = hl.t0, trt1 = hl.t1), 
              p.value = c(trt0 = pval.t0, trt1 = pval.t1), Df = c(Df), 
              plot = p)
  class(res) = "calibrate.trtsel"
  return(res)
}

