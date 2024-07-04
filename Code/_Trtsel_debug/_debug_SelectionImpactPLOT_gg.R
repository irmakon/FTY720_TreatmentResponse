# to be copy-pasted within the following function:
# fixInNamespace("SelectionImpactPLOT_gg", pos = "package:TreatmentSelection", editor = "internal")

function (x, ci, ci.bounds, get.F, fixed.values, conf.bands, 
          rho, xlab, ylab, xlim, ylim, main, clr = "gray50") 
{
  risk.t0 <- x$derived.data$fittedrisk.t0
  risk.t1 <- x$derived.data$fittedrisk.t1
  trt.effect <- x$derived.data$trt.effect
  if (x$model.fit$outcome == "time-to-event") {
    event = rep(0, nrow(x$derived.data))
    event.name = x$treatment.name
  }
  else {
    event <- x$derived.data[[as.character(x$formula[[2]])]]
    event.name = as.character(x$formula[[2]])
  }
  trt <- x$derived.data[[x$treatment.name]]
  n = length(trt)
  F.D <- get.F(trt.effect, event, trt, rho = rho) * 100
  theta.curve <- EventRateVec(risk.t0, risk.t1, F.D, rho, event, 
                              trt)
  lty = 1
  mydata = data.frame(theta.curve, F.D, lty)
  mydata = mydata[with(mydata, order(F.D)), ]
  allMeasures <- x$functions$get.summary.measures(x$derived.data, 
                                                  event.name = event.name, treatment.name = x$treatment.name, 
                                                  rho = rho, d = x$model.fit$thresh)
  avglines <- cbind(allMeasures$ER.trt0.mod, sort(F.D), 4)
  avglines <- rbind(avglines, cbind(allMeasures$ER.trt1.mod, 
                                    sort(F.D), 3))
  avglines = data.frame(avglines)
  names(avglines) = names(mydata)
  mydata <- rbind(mydata, avglines)
  cen <- mean(c(min(theta.curve, na.rm = TRUE), max(theta.curve, 
                                                    na.rm = TRUE)))
  ran <- max(theta.curve, na.rm = TRUE) - min(theta.curve, 
                                              na.rm = TRUE)
  if (substr(ci, 1, 1) == "v") {
    cen <- mean(c(min(c(theta.curve, ci.bounds), na.rm = TRUE), 
                  max(c(theta.curve, ci.bounds), na.rm = TRUE)))
    ran <- max(c(theta.curve, ci.bounds), na.rm = TRUE) - 
      min(c(theta.curve, ci.bounds), na.rm = TRUE)
  }
  ran <- ran * 1.1
  mylim <- c(cen - ran/2, cen + ran/2)
  if (is.null(xlab)) 
    xlab <- "% population recommended treatment"
  if (is.null(ylab)) 
    ylab <- "Rate of outcome under treatment policy"
  if (is.null(xlim)) 
    xlim <- c(0, 100)
  if (is.null(ylim)) 
    ylim <- mylim
  if (is.null(main)) 
    main <- "Selection impact curve"
  p <- ggplot()
  p <- p + xlab(xlab) + ylab(ylab) + ylim(ylim[1], ylim[2]) + 
    ggtitle(main)
  p <- p + scale_x_continuous(limits = xlim)
  p <- p + coord_cartesian(expand = F)
  if (!is.null(ci.bounds)) {
    ci.bounds <- matrix(ci.bounds, ncol = length(fixed.values), 
                        nrow = 2)
    if (substr(ci, 1, 1) == "v") {
      index.fix <- (fixed.values <= max(F.D) & fixed.values >= 
                      min(F.D))
      width = 5
    }
    else {
      width = 0.05
      index.fix <- (fixed.values <= max(theta.curve) & 
                      fixed.values >= min(theta.curve))
    }
    p <- shade_gg(p, ci.bounds[, index.fix], fixed.values[index.fix], 
                  type = substr(ci, 1, 1), bands = conf.bands, lty, 
                  width = width, clr = clr[1])
  }
  p <- p + geom_step(data = mydata[(1:(n)), ], 
                     aes(x = F.D, y = theta.curve), direction = "vh", 
                     color = clr[1], size = 1.1)
  
  p <- p + geom_line(data = mydata[-c(1:(n)), ], 
                     aes(x = F.D, y = theta.curve, linetype = factor(lty)), 
                     size = 1)
  p <- p + scale_linetype_manual(name = "Event Rate", breaks = c("3", 
                                                                 "4"), values = c(3, 4), labels = c("treat all", "treat none"))
  p <- p + theme(legend.title = element_blank(), 
                 legend.key.size = unit(1.5, "lines"), 
                 legend.background = element_blank(), 
                 legend.key = element_rect(fill = NA),
                 text = element_text(size = 7),
                 legend.position = c(0.98, 0.02),
                 legend.justification = c(0.98, 0.02),
                 legend.direction = "horizontal",
                 legend.box = "horizontal",
                 legend.margin = margin(0.1, 0.1, 0.1, 0.1, unit = "pt"))
  print(p)
  return(list(p = p))
}
