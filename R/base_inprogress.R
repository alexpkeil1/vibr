

# .Haw <- function(gn, ga, gb){
#   # I(A<u(w))g(a-delta|w)/g(a|w) + I(a>=u(w)-delta)
#   Haw <- ifelse(gn>0, ga/gn, 0) + as.numeric(gb == 0)
#   Haw
# }

.OneStepTmleContBounded <- function(Y,
                             Qinit,
                             Qdawinit,
                             Haw,
                             Hdaw,
                             isbin=FALSE,
                             weighted=FALSE,
                             wt=wt,
                             .link=.identity,
                             .ilink=.invidentity
){
  # Q_update(A,W) = Q(A,W) + eps*H(A,W)
  # Q_update(A,W) = expit(logit(Q(A,W)) + eps*H(A,W))
  # 1. estimate epsilon
  if(.link(.3)==.ilink(.3)) fam <- gaussian()
  if(.link(.3)!=.ilink(.3)) fam <- binomial()
  if(weighted){
    epsk <- as.numeric(glm(Y~ offset(.link(Qinit)), weights = Haw*wt, family=fam)$coefficients[1])
  } else{
    epsk <- as.numeric(glm(Y~ -1 + offset(.link(Qinit)) + Haw, weights=wt, family=fam)$coefficients[1])
  }
  # 2. update Qk
  Qk1 <- .ilink(.link(Qinit) + epsk*Haw)
  Qdawk1 <- .ilink(.link(Qdawinit) + epsk*Hdaw)
  cbind(Qk1, Qdawk1)
}

.plotfun <- function(){
  data(metals, package="qgcomp")
  XYlist = list(X=metals[,1:23], Y=metals$y)
  set.seed(123123)
  (vimp <- varimp(data.frame(XYlist$X),XYlist$Y, delta=.1, Y_learners=.default_continuous_learners_big(),
                  Xdensity_learners=.default_density_learners_big(), Xbinary_learners=.default_binary_learners_big(),
                  verbose=FALSE, estimator="TMLE", estimand="diff", updatetype = "unweighted"))
  (vimp2 <- varimp(data.frame(XYlist$X),XYlist$Y, delta=.1, Y_learners=.default_continuous_learners_big(),
                  Xdensity_learners=.default_density_learners_big(), Xbinary_learners=.default_binary_learners_big(),
                  verbose=FALSE, estimator="AIPW", estimand="diff", updatetype = "unweighted"))

  gf <-  vimp$gfits[[2]]
  target <- gf$training_task$Y
  gf$training_task$column_names
  gf$training_task$nodes$outcome
  pred <- gf$predict()[[1]]
  plot(target,pred, cex=.3, pch=19)
  points(target+0.1,pred, cex=.3, pch=19, col="red")

  gf2 <-  vimp$gfits[[2]]
  target2 <- gf2$training_task$Y
  pred2 <- gf2$predict()[[1]]
  # should be a note that p is uniform
  hist(pred2)
  points(target+0.1,pred, cex=.3, pch=19, col="red")

}


