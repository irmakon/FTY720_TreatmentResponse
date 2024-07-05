# to be copy-pasted within the following function:
# fixInNamespace("get.summary.measures.cohort.survival", pos = "package:TreatmentSelection", editor = "internal")

function (data, event.name, treatment.name, rho, d = 0) 
{
  t0 <- data$prediction.time[1]
  trt = data[[treatment.name]]
  trt.effect <- data$trt.effect
  wi <- data$censoring.weights
  tmp <- data[,1:2]
  stime = tmp[, 1]
  status = tmp[, 2]
  if (is.null(data[["rec.trt"]])) {
    neg <- data$rec.no.trt
    pos <- 1 - neg
  }
  else {
    pos <- data$rec.trt
    neg <- 1 - pos
  }
  p.rec.no.trt <- mean(neg)
  p.rec.trt <- mean(pos)
  num <- (wi * I(stime < t0) * I(neg))
  den <- (wi * I(neg))
  if (sum(den[trt == 0]) == 0 | sum(den[trt == 1]) == 0) {
    B.neg.emp <- 0
    ER.neg.emp <- 0
  }
  else {
    B.neg.emp <- sum(num[trt == 1])/sum(den[trt == 1]) - 
      sum(num[trt == 0])/sum(den[trt == 0])
    ER.neg.emp <- sum(num[trt == 0])/sum(den[trt == 0])
  }
  B.neg.mod <- ifelse(sum(neg) > 0, -mean(trt.effect[neg == 
                                                       1]), 0)
  ER.neg.mod <- ifelse(sum(neg) > 0, mean(data$fittedrisk.t0[neg == 
                                                               1]), 0)
  num <- (wi * I(stime < t0) * pos)
  den <- (wi * pos)
  if (sum(den[trt == 0]) == 0 | sum(den[trt == 1]) == 0) {
    B.pos.emp <- 0
    ER.pos.emp <- 0
  }
  else {
    B.pos.emp <- sum(num[trt == 0])/sum(den[trt == 0]) - 
      sum(num[trt == 1])/sum(den[trt == 1])
    ER.pos.emp <- sum(num[trt == 1])/sum(den[trt == 1])
  }
  B.pos.mod <- ifelse(sum(pos) > 0, mean(trt.effect[pos == 
                                                      1]), 0)
  ER.pos.mod <- ifelse(sum(pos) > 0, mean(data$fittedrisk.t1[pos == 
                                                               1]), 0)
  if (is.null(data[["rec.trt"]])) {
    Theta.emp <- B.neg.emp * p.rec.no.trt
    Theta.mod <- B.neg.mod * p.rec.no.trt
  }
  else {
    Theta.emp <- B.pos.emp * p.rec.trt
    Theta.mod <- B.pos.mod * p.rec.trt
  }
  p0.hat <- sum(wi * I(stime < t0) * (1 - trt))/sum(wi * (1 - 
                                                            trt))
  p1.hat <- sum(wi * I(stime < t0) * trt)/sum(wi * trt)
  Var.Delta <- mean((trt.effect - (p0.hat - p1.hat))^2)
  event = 0
  delta.F <- get.F.cohort(trt.effect, event, trt, rho = NULL)
  ooo <- order(trt.effect)
  s.delta.F <- delta.F[ooo]
  TG <- sum(diff(c(0, s.delta.F)) * abs(sort(trt.effect) - 
                                          (p0.hat - p1.hat)))
  ER.trt0.emp = p0.hat
  ER.trt0.mod = mean(data$fittedrisk.t0)
  ER.trt1.emp = p1.hat
  ER.trt1.mod = mean(data$fittedrisk.t1)
  ER.mkrbased.emp = ER.pos.emp * p.rec.trt + ER.neg.emp * p.rec.no.trt
  ER.mkrbased.mod = ER.pos.mod * p.rec.trt + ER.neg.mod * p.rec.no.trt
  list(p.rec.no.trt = p.rec.no.trt, p.rec.trt = p.rec.trt, 
       B.neg.emp = B.neg.emp, B.neg.mod = B.neg.mod, B.pos.emp = B.pos.emp, 
       B.pos.mod = B.pos.mod, Theta.emp = Theta.emp, Theta.mod = Theta.mod, 
       Var.Delta = Var.Delta, TG = TG, ER.trt0.emp = ER.trt0.emp, 
       ER.trt0.mod = ER.trt0.mod, ER.trt1.emp = ER.trt1.emp, 
       ER.trt1.mod = ER.trt1.mod, ER.mkrbased.emp = ER.mkrbased.emp, 
       ER.mkrbased.mod = ER.mkrbased.mod)
}
