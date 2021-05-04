

# .Haw <- function(gn, ga, gb){
#   # I(A<u(w))g(a-delta|w)/g(a|w) + I(a>=u(w)-delta)
#   Haw <- ifelse(gn>0, ga/gn, 0) + as.numeric(gb == 0)
#   Haw
# }

.OneStepTmleContBounded <- function(Y,Qinit,Qdawinit,Haw,Hdaw,isbin=FALSE){
  # IN PROGRESS
  # will scale outcome to unit length and use logistic regression to keep bounded
  yb = .bound_zero_one(Y)
  Ybound = yb[[1]]
  epsk <- as.numeric(glm(Y~ -1 + offset(.logit(Qinit)) + Haw, family=binomial(link="logit"))$coefficients[1])
  # 2. update Qk
  Qk1 <- .expit(.logit(Qinit) + epsk*Haw)
  Qdawk1 <- .expit(.logit(Qdawinit) + epsk*Hdaw)
  cbind(Qk1, Qdawk1)

}
