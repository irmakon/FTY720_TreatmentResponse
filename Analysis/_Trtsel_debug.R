# install.packages("TreatmentSelection")
# library("TreatmentSelection")
      # data("tsdata")
      # 
      # mymod <- glm(event~trt*(Y2), data= tsdata, family = binomial("logit"))
      # 
      # tsdata$fitted.t0 <- predict(mymod, newdata=data.frame(trt = 0,Y1 = tsdata$Y1, Y2 = tsdata$Y2), type = "response")
      # tsdata$fitted.t1 <- predict(mymod, newdata=data.frame(trt = 1,Y1 = tsdata$Y1, Y2 = tsdata$Y2), type = "response")
      # 
      # tsdata_new <- tsdata
      # tsdata_new$days <- sample(1:750, 1000, replace = T)

fixInNamespace("plot.trtsel", pos = "package:TreatmentSelection", editor = "internal")

fixInNamespace("get.summary.measures.cohort.survival", pos = "package:TreatmentSelection", editor = "internal")
# line 7 to 
# tmp <- data[,1:2]


      # myfitted.trtsel <- trtsel( Surv(days, event)~trt, treatment.name = "trt", 
      #                            data = tsdata_new,
      #                            fittedrisk.t0 = "fitted.t0",
      #                            fittedrisk.t1 = "fitted.t1",
      #                            study.design = "RCT", 
      #                            default.trt = "trt all",
      #                            prediction.time = 500)

      # The below throws in error if there are no fixes
      # plot(myfitted.trtsel, plot.type = "risk")

fixInNamespace("one.boot.plot", pos = "package:TreatmentSelection", editor = "internal")
# lines 3 and 18 to 
# if (x$model.fit$outcome == "time-to-event") {

fixInNamespace("shade_gg", pos = "package:TreatmentSelection", editor = "internal")
# function (p, bounds, fixed, type, lty = 1, bands, width = 5, clr = "gray50")

# p <- p + geom_polygon(data = mydat, alpha = 0.4, 
#                       aes(x = bounds, y = fixed), fill = clr)

fixInNamespace("predcurvePLOT_gg", pos = "package:TreatmentSelection", editor = "internal")
# line 21 to 
# if (x$model.fit$outcome == "time-to-event") {
# add line
# p <- p + coord_cartesian(expand = F)
# revise the first theme to:
# theme(legend.title = element_blank(), legend.key.size = unit(1.5,
#                                                              "lines"), text = element_text(size = 8),
#       legend.position = c(0.98, 0.02),
#       legend.justification = c(0.98, 0.02),
#       legend.direction = "horizontal",
#       legend.box = "horizontal",
#       legend.margin = margin(0.1, 0.1, 0.1, 0.1, unit = "pt"))

# p <- p + scale_color_manual(values = colcodes, labels = trt.names) + 
# p <- p + geom_step(data = mydata, aes(x = Fy, y = risk,
#                                       color = factor(trt)), #linetype = factor(trt)),
#                    direction = "vh", size = 1.1)

# p <- shade_gg(ggplot(mydata), ci.bounds[1:2, index.fix.t0],
#               fixed.values[index.fix.t0], type = substr(ci, 1,
#                                                         1), bands = conf.bands, width = width, clr = colcodes[2])
# p <- shade_gg(p, ci.bounds[3:4, index.fix.t1], fixed.values[index.fix.t1] +
#                 offset, type = substr(ci, 1, 1), bands = conf.bands, clr = colcodes[1],
#               width = width)

      # plot(myfitted.trtsel, plot.type = "risk")

# No need to apply the below: CDF plot is not used
# fixInNamespace("CDFdeltaPLOT_gg", pos = "package:TreatmentSelection", editor = "internal")
# line 5 to 
# if (x$model.fit$outcome == "time-to-event") {
# add line
# p <- p + coord_cartesian(expand = F)

      # plot(myfitted.trtsel, plot.type = "cdf")

fixInNamespace("SelectionImpactPLOT_gg", pos = "package:TreatmentSelection", editor = "internal")
# line 7 to 
# if (x$model.fit$outcome == "time-to-event") {
# add line
# p <- p + coord_cartesian(expand = F)

# p <- p + geom_line(data = mydata[-c(1:(n)), ], aes(x = F.D, 
#                                                    y = theta.curve, linetype = factor(lty)), size = 1)
# p <- p + geom_step(data = mydata[(1:(n)), ], aes(x = F.D, 
#                                                  y = theta.curve), direction = "vh", color = colcodes[1], size = 1.1)

# width = width, clr = colcodes[1])

# theme(legend.title = element_blank(), legend.key.size = unit(1.5,
#                                                              "lines"), text = element_text(size = 8),
#       legend.position = c(0.98, 0.02),
#       legend.justification = c(0.98, 0.02),
#       legend.direction = "horizontal",
#       legend.box = "horizontal",
#       legend.margin = margin(0.1, 0.1, 0.1, 0.1, unit = "pt"))

# Note that get.summary.measures need to be corrected before creation of the trtsel object

      # plot(myfitted.trtsel, plot.type = "selection impact")

fixInNamespace("trteffectPLOT_gg", pos = "package:TreatmentSelection", editor = "internal")
# line 6 to 
# if (x$model.fit$outcome == "time-to-event") {
# theme(legend.title = element_blank(), legend.key.size = unit(1.5,
#                                                              "lines"), text = element_text(size = 8),
#       legend.position = c(0.98, 0.02),
#       legend.justification = c(0.98, 0.02),
#       legend.direction = "horizontal",
#       legend.box = "horizontal",
#       legend.margin = margin(0.1, 0.1, 0.1, 0.1, unit = "pt"))

      # plot(myfitted.trtsel, plot.type = "treatment effect")

fixInNamespace("calibrate.trtsel", pos = "package:TreatmentSelection", editor = "internal")
# lines 12, 103(104), 309(311), 350(353) to 
# if (x$model.fit$outcome == "time-to-event") {

      # calibrate(myfitted.trtsel)
      # calibrate(myfitted.trtsel, plot.type="risk.t1")
      # calibrate(myfitted.trtsel, plot.type = "risk.t0")

fixInNamespace("evaluate.trtsel", pos = "package:TreatmentSelection", editor = "internal")
# NO NEED TO add at 37 
#  if(!(x$model.fit$outcome == "time-to-event"))
# line 70 (72) to 
# if (x$model.fit$outcome == "time-to-event")
# comment out line 77
# #summary.measures <- summary.measures - bias

fixInNamespace("one.boot.eval", pos = "package:TreatmentSelection", editor = "internal")
# lines 5, NO NEED TO 26(25) to 
#  if (!is.null(prediction.time)) {

      # evaluate(myfitted.trtsel, bias.correct = F)

fixInNamespace("trtsel.boot", pos = "package:TreatmentSelection", editor = "internal")
# add at 8 
# coef <- NULL
