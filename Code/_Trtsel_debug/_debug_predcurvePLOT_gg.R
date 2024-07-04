# to be copy-pasted within the following function:
# fixInNamespace("predcurvePLOT_gg", pos = "package:TreatmentSelection", editor = "internal")

function (x, ci, ci.bounds, get.F, fixed.values, conf.bands, 
          rho, trt.names, xlab, ylab, xlim, ylim, main, show.marker.axis, 
          offset = 0.01, clr) 
{ 
  fittedrisk.t0 <- x$derived.data$fittedrisk.t0
  fittedrisk.t1 <- x$derived.data$fittedrisk.t1
  if (length(x$model.fit$marker.names) == 1) {
    marker <- x$derived.data[, x$model.fit$marker]
    if (is.null(xlab)) 
      xlab <- "% population below marker value"
    if (is.null(ylab)) 
      ylab <- "Risk given marker"
  }
  else {
    marker <- x$derived.data$trt.effect
    if (is.null(xlab)) 
      xlab <- "% population below treatment effect"
    if (is.null(ylab)) 
      ylab <- "Risk given treatment effect"
  }
  if (x$model.fit$outcome == "time-to-event") {
    event = 0
  }
  else {
    event <- x$derived.data[[as.character(x$formula[[2]])]]
  }
  trt <- x$derived.data[[x$treatment.name]]
  F.Y <- get.F(marker, event, trt, rho = rho) * 100
  n = length(fittedrisk.t0)
  mydata <- data.frame(risk = c(fittedrisk.t0, fittedrisk.t1), 
                       trt = c(rep(1, n), rep(0, n)), Fy = rep(F.Y, 2))
  mydata <- mydata[!duplicated(mydata), ]
  Fy <- risk <- NULL
  if (!is.null(ci.bounds)) {
    ci.bounds <- matrix(ci.bounds, ncol = length(fixed.values), 
                        nrow = 4)
    if (substr(ci, 1, 1) == "h") {
      width = 0.05
      index.fix.t0 <- (fixed.values <= max(fittedrisk.t0[trt == 
                                                           0]) & fixed.values >= min(fittedrisk.t0[trt == 
                                                                                                     0]))
      index.fix.t1 <- (fixed.values <= max(fittedrisk.t1[trt == 
                                                           1]) & fixed.values >= min(fittedrisk.t1[trt == 
                                                                                                     1]))
    }
    else {
      width = 5
      index.fix.t0 <- (fixed.values <= max(F.Y[trt == 0]) & 
                         fixed.values >= min(F.Y[trt == 0]))
      index.fix.t1 <- (fixed.values <= max(F.Y[trt == 1]) & 
                         fixed.values >= min(F.Y[trt == 1]))
    }
    p <- shade_gg(ggplot(mydata), ci.bounds[1:2, index.fix.t0],
                  fixed.values[index.fix.t0], type = substr(ci, 1,
                                                            1), bands = conf.bands, clr = clr[2], 
                  width = width)
    p <- shade_gg(p, ci.bounds[3:4, index.fix.t1], fixed.values[index.fix.t1] +
                    offset, type = substr(ci, 1, 1), bands = conf.bands, clr = clr[1],
                  width = width)
  }
  else {
    p <- ggplot(mydata)
  }
  if (is.null(xlim)) 
    xlim <- c(0, 100)
  if (is.null(ylim)) 
    ylim <- c(0, 1)
  if (is.null(main)) 
    main <- "Risk curves by treatment"
  
  breaks = seq(xlim[1], xlim[2], length.out = 5)
  p <- p + geom_smooth(data = mydata, aes(x = Fy, y = risk,
                                        color = factor(trt)), #linetype = factor(trt)),
                     direction = "vh", size = 1.1, se = F)
  p <- p + scale_color_manual(values = clr, labels = trt.names) 
  p <- p + xlab(xlab) + ylab(ylab) + ylim(ylim[1], ylim[2]) + 
    ggtitle(main) +
  # p <- p + scale_linetype_manual(values = c(2, 1), labels = trt.names) + 
    theme(legend.title = element_blank(),
          legend.background = element_blank(), 
          legend.key = element_rect(fill = NA),
          legend.key.size = unit(1.5, "lines"), text = element_text(size = 7),
          legend.position = c(0.98, 0.02),
          legend.justification = c(0.98, 0.02),
          legend.direction = "horizontal",
          legend.box = "horizontal",
          legend.margin = margin(0.1, 0.1, 0.1, 0.1, unit = "pt"))
  p <- p + scale_x_continuous(breaks = breaks, limits = xlim)
  p <- p + coord_cartesian(expand = F)
  if (show.marker.axis) {
    p <- p + theme(plot.margin = unit(c(1, 1, 4, 1), "lines"))
    p <- p + annotation_custom(grob = xaxisGrob(at = breaks, 
                                                label = round(quantile(marker, prob = breaks/100), 
                                                              1)), xmin = 0, xmax = 1, ymin = ylim[1] - diff(ylim) * 
                                 0.25, ymax = ylim[1] - diff(ylim) * 0.25)
    p <- p + annotation_custom(grob = textGrob(label = "marker value"), 
                               xmin = mean(xlim), xmax = mean(xlim), ymin = ylim[1] - 
                                 diff(ylim) * 0.4, ymax = ylim[1] - diff(ylim) * 
                                 0.4)
    plot.new()
    gt <- ggplotGrob((p))
    gt$layout$clip[gt$layout$name == "panel"] <- "off"
    grid.draw(gt)
  }
  else {
    print(p)
  }
  list(p = p)
}

