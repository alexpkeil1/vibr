library("vibr")
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


Y = dat$telomean
Xcr = dat[,c(mixturela)]
X = dat[,c(mixturela, "ridageyr")]
X = data.frame(cbind(log(dat[,c(mixturela)]), ridageyr=dat$ridageyr))

lnr  = Lrnr_density_discretize$new(name=nm, categorical_learner = Lrnr_multinom$new(trace=FALSE), n_bins = 5, bin_method="equal.mass")

vi0 <- vibr::varimp(X,Y,delta=0.1,
                   Y_learners = .default_continuous_learners(),
                   Xdensity_learners = .default_density_learners(),
                   Xbinary_learners = list(Lrnr_stepwise$new()), estimator="TMLE")
vi1 <- vibr::varimp_refit(vi0,X,Y,delta=0.1, estimator="AIPW")
vi2 <- vibr::varimp_refit(vi0,X,Y,delta=0.1, estimator="GCOMP")
vi3 <- vibr::varimp_refit(vi0,X,Y,delta=0.1, estimator="IPW")
vi0
vi1
vi2
vi3

ridx <- which.max(vi$res$est)
vi$res$est[ridx]
vi$res[ridx,]
predmod <- vi$gfits[[ridx]]
hist(X[,ridx])
plot(X[,ridx], predmod$predict()[[1]])
crdat = data.frame(telomean=dat$telomean, X[,ridx])
names(crdat)[2] <- names(X)[ridx]
summary(glm(telomean~ ., data=crdat))

ridx <- which.min(vi$res$est)
vi$res$est[ridx]
rbind(tmle=vi$res[ridx,], aipw=vi1$res[ridx,], gcomp=vi2$res[ridx,], ipw=vi3$res[ridx,])
predmod <- vi$gfits[[ridx]]
predmod$learner_fits$Stack
hist(X[,ridx])
plot(X[,ridx], predmod$predict()[[1]])
crdat = data.frame(telomean=dat$telomean, X[,ridx])
names(crdat)[2] <- names(X)[ridx]
summary(glm(telomean~ ., data=crdat))


(lm(telomean~lbxf09la, data=dat[,c("telomean", mixture)]))
(lm(telomean~lbxf01la, data=dat[,c("telomean", mixture)]))
(lm(telomean~lbxf06la, data=dat[,c("telomean", mixture)]))
(lm(telomean~lbxf02la, data=dat[,c("telomean", mixture)]))
(lm(telomean~lbx153la, data=dat[,c("telomean", mixture)]))
(lm(telomean~lbx189la, data=dat[,c("telomean", mixture)]))


(qgcft <- qgcomp::qgcomp(telomean~., expnms = mixture, data=dat[,c("telomean", mixture)]))
(qgcft <- qgcomp::qgcomp(telomean~., expnms = pcbsla,  data=dat[,c("telomean", mixture)]))
(qgcft <- qgcomp::qgcomp(telomean~., expnms = furansla, data=dat[,c("telomean", mixture)]))
(qgcft <- qgcomp::qgcomp(telomean~., expnms = dioxinsla, data=dat[,c("telomean", mixture)]))

(qgcftadj <- qgcomp::qgcomp.noboot(telomean~.+I(ridageyr^2), expnms = mixture, data=dat[,c("telomean", mixture, "ridageyr")]))
(qgcftadj <- qgcomp::qgcomp.noboot(telomean~.+I(ridageyr^2), expnms = pcbsla, data=dat[,c("telomean", mixture, "ridageyr")]))
(qgcftadj <- qgcomp::qgcomp.noboot(telomean~.+I(ridageyr^2), expnms = furansla, data=dat[,c("telomean", mixture, "ridageyr")]))
(qgcftadj <- qgcomp::qgcomp.noboot(telomean~.+I(ridageyr^2), expnms = dioxinsla, data=dat[,c("telomean", mixture, "ridageyr")]))
qgcftadj$fit


wqft <- gWQS::gwqs(telomean~wqs+ridageex+I(ridageex^2), rs=TRUE, n_vars=6, mix_name=mixture, data=dat[,c("telomean", mixture, "ridageex")])
summary(wqft)
head(wqft$final_weights)

wqft2<- gWQS::gwqs(telomean~wqs, b1_pos = FALSE, rs=TRUE, n_vars=6, mix_name=mixture, data=dat[,c("telomean", mixture)])
summary(wqft2)

