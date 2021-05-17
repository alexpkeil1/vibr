library("vibr")
library("ggplot2")
pcbs1 <- paste0("lbx", sprintf("%03d",c(74, 99, 105, 118,138,153,156,157,167,170,180,187,189,194,196)))
pcbs2 <- tolower(c(paste0("lbd", sprintf("%03d",c(199))), "LBXPCB", "LBXHXC"))
pcbs1la <- paste0(paste0("lbx", sprintf("%03d",c(74, 99, 105, 118,138,153,156,157,167,170,180,187,189,194,196, 199))), "la")
pcbs2la <- paste0(tolower(c("LBXPCB", "LBXHXC")), "la")
pcbs <- c(pcbs1,pcbs2)
pcbsla <- c(pcbs1la,pcbs2la)
dioxins <- tolower(c("LBXTCD", "LBXD01", "LBXD02", "LBXD03", "LBXD04", "LBXD05", "LBXD07"))
dioxinsla <- paste0(dioxins, "la")
furans <- tolower(paste0(c("LBXF0"), 1:9))
furansla <- paste0(furans, "la")

mixturela = c(pcbsla,dioxinsla, furansla)

dat <- read.csv("/Users/akeil/EpiProjects/NHANES/pcbs_furans_dioxins/mitro_data_cc.csv")
dim(dat)

Y = dat$telomean
Xcr = dat[,c(mixturela)]
X = dat[,c(mixturela, "ridageyr")]
V =dat[,c("wtspo2yr"),drop=FALSE]
V$wtspo2yr = V$wtspo2yr/(mean(V$wtspo2yr))
#X = data.frame(cbind(log(dat[,c(mixturela)]), ridageyr=dat$ridageyr))
summary(vibr:::.scale_continuous(X[,2, drop=FALSE]))

Xi = X[,,drop=FALSE]
glm.fit(cbind(rep(1, length(Y)),Xi),Y,weights=V$wtspo2yr)$coef
glm.fit(cbind(rep(1, length(Y)),Xi),Y)$coef

pass <- function(){
  tsk = sl3_Task$new(data = data.frame(cbind(vibr:::.scale_continuous(Xi),V)),
                     covariates = names(Xi)[1:2],
                     outcome=names(Xi)[3],
                     outcome_type="continuous",
                     weights="wtspo2yr")
  head(tsk$X)
  head(tsk$Y)
  mean(V$wtspo2yr)
  lnr = Lrnr_polspline$new(name="polymars", classify=FALSE)
  trn <- lnr$train(tsk)
  typeof(trn$get_outcome_type(tsk)$format(tsk$Y))
  trn$predict()
  trn$fit_object$call
}

(ncores <- future::availableCores())
future::plan("multisession", workers=ncores/2)


(viNULL <- vibr::varimp(X=Xi,Y=Y, V=V,delta=0.01, weights="wtspo2yr",
                    Y_learners = .default_continuous_learners()[1:4],
                    Xdensity_learners = .default_density_learners_big()[c(1)],
                    Xbinary_learners = list(Lrnr_stepwise$new()), estimator="IPW"))

vi0 <- vibr::varimp(X=Xi,Y=Y, V=V,delta=0.01, weights="wtspo2yr",
                    Y_learners = .default_continuous_learners_big()[1:5],
                    Xdensity_learners = .default_density_learners_big()[c(1:3,8)],
                    Xbinary_learners = list(Lrnr_stepwise$new()), estimator="TMLE")
vi1 <- vibr::varimp_refit(vi0,Xi,Y,delta=0.01, estimator="AIPW")
vi2 <- vibr::varimp_refit(vi0,Xi,Y,delta=0.01, estimator="GCOMP")
vi3 <- vibr::varimp_refit(vi0,Xi,Y,delta=0.01, estimator="IPW")
vi0
vi1
vi2
vi3


round(cor(X), 2)
X = dat[,c(mixturela, "ridageyr")]
plotshift_scatter(vi0, Acol=12, Bcol=2, delta=.01)
plotshift_dens(vi0, X, Acol=c(12), delta=.01)
plotshift_dens(vi0, X, Acol=c(12,1:11,13:35), delta=.01)


plotshift_scatter(vi0, Acol=13, Bcol=2)
plotshift_dens(vi0, X, Acol=1, delta=.01)
plotshift_wt(vi0, X, Acol=c(12), delta=.01)
plotshift_wt(vi0, X, Acol=c(12,1:11), delta=0.00001)
dx_dens(vi0, X, Acol=c(1), delta=.0001)



ggplot(data=data.frame(x=vi0$rank, y=vi1$rank), aes(x=x,y=y)) + geom_smooth(method = "lm", se=FALSE) + geom_point() # tmle, aipw
ggplot(data=data.frame(x=vi0$rank, y=vi2$rank), aes(x=x,y=y)) + geom_smooth(method = "lm", se=FALSE) + geom_point() # tmle, gcomp
ggplot(data=data.frame(x=vi0$rank, y=vi3$rank), aes(x=x,y=y)) + geom_smooth(method = "lm", se=FALSE) + geom_point() # tmle, ipw
ggplot(data=data.frame(x=vi1$rank, y=vi2$rank), aes(x=x,y=y)) + geom_smooth(method = "lm", se=FALSE) + geom_point() # aipw, gcomp
ggplot(data=data.frame(x=vi1$rank, y=vi3$rank), aes(x=x,y=y)) + geom_smooth(method = "lm", se=FALSE) + geom_point() # aipw, ipw
ggplot(data=data.frame(x=vi2$rank, y=vi3$rank), aes(x=x,y=y)) + geom_smooth(method = "lm", se=FALSE) + geom_point() # gcomp, ipw

ggplot(data=data.frame(x=vi0$res$est, y=vi1$res$est), aes(x=x,y=y)) + geom_smooth(method = "lm", se=FALSE) + geom_point() # tmle, aipw
ggplot(data=data.frame(x=vi0$res$est, y=vi2$res$est), aes(x=x,y=y)) + geom_smooth(method = "lm", se=FALSE) + geom_point() # tmle, gcomp
ggplot(data=data.frame(x=vi0$res$est, y=vi3$res$est), aes(x=x,y=y)) + geom_smooth(method = "lm", se=FALSE) + geom_point() # tmle, ipw
ggplot(data=data.frame(x=vi1$res$est, y=vi2$res$est), aes(x=x,y=y)) + geom_smooth(method = "lm", se=FALSE) + geom_point() # aipw, gcomp
ggplot(data=data.frame(x=vi1$res$est, y=vi3$res$est), aes(x=x,y=y)) + geom_smooth(method = "lm", se=FALSE) + geom_point() # aipw, ipw
ggplot(data=data.frame(x=vi2$res$est, y=vi3$res$est), aes(x=x,y=y)) + geom_smooth(method = "lm", se=FALSE) + geom_point() # gcomp, ipw


# same effect direction
top10 = which(vi0$rank %in% 1:10)
cor(as.matrix(cbind(tmle=sign(vi0$res$est),
                    aipw=sign(vi1$res$est),
                    gcomp=sign(vi2$res$est),
                    ipw=sign(vi3$res$est))))

cor(as.matrix(cbind(tmle=sign(vi0$res$est[top10]),
                    aipw=sign(vi1$res$est[top10]),
                    gcomp=sign(vi2$res$est[top10]),
                    ipw=sign(vi3$res$est[top10]))))


cor(as.matrix(cbind(tmle=vi0$rank, aipw=vi1$rank, gcomp=vi2$rank, ipw=vi3$rank)))
cor(as.matrix(cbind(tmle=vi0$res$est, aipw=vi1$res$est, gcomp=vi2$res$est, ipw=vi3$res$est)))


(qgcft <- qgcomp::qgcomp(telomean~., expnms = mixturela, data=dat[,c("telomean", mixturela)]))
(qgcft <- qgcomp::qgcomp(telomean~., expnms = pcbsla,  data=dat[,c("telomean", mixturela)]))
(qgcft <- qgcomp::qgcomp(telomean~., expnms = furansla, data=dat[,c("telomean", mixturela)]))
(qgcft <- qgcomp::qgcomp(telomean~., expnms = dioxinsla, data=dat[,c("telomean", mixturela)]))

(qgcftadj <- qgcomp::qgcomp.noboot(telomean~.+I(ridageyr^2), expnms = mixturela, data=dat[,c("telomean", mixturela, "ridageyr")]))
(qgcftadj <- qgcomp::qgcomp.noboot(telomean~.+I(ridageyr^2), expnms = pcbsla, data=dat[,c("telomean", mixturela, "ridageyr")]))
(qgcftadj <- qgcomp::qgcomp.noboot(telomean~.+I(ridageyr^2), expnms = furansla, data=dat[,c("telomean", mixturela, "ridageyr")]))
(qgcftadj <- qgcomp::qgcomp.noboot(telomean~.+I(ridageyr^2), expnms = dioxinsla, data=dat[,c("telomean", mixturela, "ridageyr")]))
qgcftadj$fit


wqft <- gWQS::gwqs(telomean~wqs+ridageex+I(ridageex^2), rs=TRUE, n_vars=6, mix_name=mixture, data=dat[,c("telomean", mixture, "ridageex")])
summary(wqft)
head(wqft$final_weights)

wqft2<- gWQS::gwqs(telomean~wqs, b1_pos = FALSE, rs=TRUE, n_vars=6, mix_name=mixture, data=dat[,c("telomean", mixture)])
summary(wqft2)



