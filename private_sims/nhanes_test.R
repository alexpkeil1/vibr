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

  vijoint <- vibr:::.varimp_gcomp_joint(X=Xi,Y=Y, V=V,
                                         expnms=c(mixturela),
                                        delta=.01,
                                        weights="wtspo2yr",
                                        estimand="diff",
                          Y_learners = .default_continuous_learners_big()[1:4],
                          Xdensity_learners = NULL,
                          Xbinary_learners = NULL)
  vijoint <- vibr:::.attach_misc(vijoint, scale_continuous=TRUE, delta=vijoint$delta)
  vijoint


  vijointb <- vibr:::.varimp_gcomp_joint_boot(X=Xi,Y=Y, V=V,
                                             expnms=c(mixturela),
                                             delta=.01,
                                             #weights="wtspo2yr",
                                             estimand="diff",
                                             Y_learners = .default_continuous_learners_big()[1],
                                             Xdensity_learners = NULL,
                                             Xbinary_learners = NULL,
                                             B=200, verbose=FALSE)
  vijointb$est <- vibr:::.attach_misc(vijointb$est, scale_continuous=TRUE, delta=vijoint$delta)

  vijoint
  vijointb
  vijointb$est$qfit$coefficients
}

(ncores <- future::availableCores())
future::plan("multisession", workers=ncores)

gcompxfit1 <- vibr::varimp(X=Xi,Y=Y, V=V,delta=0.01, weights="wtspo2yr",
                    Y_learners = .default_continuous_learners_big()[1:5],
                    Xdensity_learners = .default_density_learners_big()[c(1:3,8)],
                    Xbinary_learners = list(Lrnr_stepwise$new()), estimator="GCOMPX",
                    xfitfolds = 1,foldrepeats = 50)

ipwxfit1 <- vibr::varimp(X=Xi,Y=Y, V=V,delta=0.01, weights="wtspo2yr",
                    Y_learners = .default_continuous_learners_big()[1:5],
                    Xdensity_learners = .default_density_learners_big()[c(1:3,8)],
                    Xbinary_learners = list(Lrnr_stepwise$new()), estimator="IPWX",
                    xfitfolds = 1,foldrepeats = 50)

gcompxfit2 <- vibr::varimp(X=Xi,Y=Y, V=V,delta=0.01, weights="wtspo2yr",
                    Y_learners = .default_continuous_learners_big()[1:5],
                    Xdensity_learners = .default_density_learners_big()[c(1:3,8)],
                    Xbinary_learners = list(Lrnr_stepwise$new()), estimator="GCOMPX",
                    xfitfolds = 2,foldrepeats = 50)

ipwxfit2 <- vibr::varimp(X=Xi,Y=Y, V=V,delta=0.01, weights="wtspo2yr",
                    Y_learners = .default_continuous_learners_big()[1:5],
                    Xdensity_learners = .default_density_learners_big()[c(1:3,8)],
                    Xbinary_learners = list(Lrnr_stepwise$new()), estimator="IPWX",
                    xfitfolds = 2,foldrepeats = 50)



save(file="/Users/akeil/temp/nhanes.Rdata", list=c( "gcompxfit1", "ipwxfit1",  "gcompxfit2", "ipwxfit2"))

tmlefit <- vibr::varimp(X=Xi,Y=Y, V=V,delta=0.01, weights="wtspo2yr",
                    Y_learners = .default_continuous_learners_big()[1:5],
                    Xdensity_learners = .default_density_learners_big()[c(1:3,8)],
                    Xbinary_learners = list(Lrnr_stepwise$new()), estimator="TMLE")
aipwfit <- vibr::varimp_refit(tmlefit,X=Xi,Y=Y,delta=0.01, estimator="AIPW")
gcompfit <- vibr::varimp_refit(tmlefit,X=Xi,Y=Y,delta=0.01, estimator="GCOMP")
ipwfit <- vibr::varimp_refit(tmlefit,X=Xi,Y=Y,delta=0.01, estimator="IPW")
save(file="/Users/akeil/temp/nhanes.Rdata", list=c( "tmlefit", "aipwfit", "gcompfit", "ipwfit","gcompxfit1", "ipwxfit1",  "gcompxfit2", "ipwxfit2"))

tmlexfit3 <- vibr::varimp(X=Xi,Y=Y, V=V,delta=0.01, weights="wtspo2yr",
                    Y_learners = .default_continuous_learners_big()[1:5],
                    Xdensity_learners = .default_density_learners_big()[c(1:3,8)],
                    Xbinary_learners = list(Lrnr_stepwise$new()), estimator="TMLEX",
                    xfitfolds = 3,foldrepeats = 50)
save(file="/Users/akeil/temp/nhanes.Rdata", list=c("tmlexfit3", "tmlefit", "aipwfit", "gcompfit", "ipwfit","gcompxfit1", "ipwxfit1",  "gcompxfit2", "ipwxfit2"))

tmlexfit1 <- vibr::varimp(X=Xi,Y=Y, V=V,delta=0.01, weights="wtspo2yr",
                    Y_learners = .default_continuous_learners_big()[1:5],
                    Xdensity_learners = .default_density_learners_big()[c(1:3,8)],
                    Xbinary_learners = list(Lrnr_stepwise$new()), estimator="TMLEX",
                    xfitfolds = 1,foldrepeats = 50)

#save(file="/Users/akeil/temp/nhanes2.Rdata", list=c("tmlexfit1"))
save(file="/Users/akeil/temp/nhanes.Rdata", list=c("tmlexfit1", "tmlexfit3", "tmlefit", "aipwfit", "gcompfit", "ipwfit","gcompxfit1", "ipwxfit1",  "gcompxfit2", "ipwxfit2"))
print(tmlefit)
print(aipwfit)
print(gcompfit)
print(ipwfit)
print(tmlexfit1)
print(tmlexfit3)

tmlexfit5 <- vibr::varimp(X=Xi,Y=Y, V=V,delta=0.01, weights="wtspo2yr",
                    Y_learners = .default_continuous_learners_big()[1:5],
                    Xdensity_learners = .default_density_learners_big()[c(1:3,8)],
                    Xbinary_learners = list(Lrnr_stepwise$new()), estimator="TMLEX",
                    xfitfolds = 5,foldrepeats = 50)

#save(file="/Users/akeil/temp/nhanes2.Rdata", list=c("tmlexfit1"))
save(file="/Users/akeil/temp/nhanes.Rdata", list=c("tmlexfit1", "tmlexfit3", "tmlexfit5", "tmlefit", "aipwfit", "gcompfit", "ipwfit","gcompxfit1", "ipwxfit1",  "gcompxfit2", "ipwxfit2"))


round(cor(X), 2)
X = dat[,c(mixturela, "ridageyr")]
plotshift_scatter(tmlefit, Acol=12, Bcol=2, delta=.01)
plotshift_dens(tmlefit, X, Acol=c(12), delta=.01)
plotshift_dens(tmlefit, X, Acol=c(12,1:11,13:35), delta=.01)


plotshift_scatter(tmlefit, Acol=13, Bcol=2)
plotshift_dens(tmlefit, X, Acol=1, delta=.01)
plotshift_wt(tmlefit, X, Acol=c(12), delta=.01)
plotshift_wt(tmlefit, X, Acol=c(12,1:11), delta=0.00001)
dx_dens(tmlefit, X, Acol=c(1), delta=.0001)



print(qgcft <- qgcomp::qgcomp(telomean~., expnms = mixturela, data=dat[,c("telomean", mixturela)]))
print(qgcft <- qgcomp::qgcomp(telomean~., expnms = pcbsla,  data=dat[,c("telomean", mixturela)]))
print(qgcft <- qgcomp::qgcomp(telomean~., expnms = furansla, data=dat[,c("telomean", mixturela)]))
print(qgcft <- qgcomp::qgcomp(telomean~., expnms = dioxinsla, data=dat[,c("telomean", mixturela)]))

print(qgcftadj <- qgcomp::qgcomp.noboot(telomean~.+I(ridageyr^2), expnms = mixturela, data=dat[,c("telomean", mixturela, "ridageyr")]))
print(qgcftadj <- qgcomp::qgcomp.noboot(telomean~.+I(ridageyr^2), expnms = pcbsla, data=dat[,c("telomean", mixturela, "ridageyr")]))
print(qgcftadj <- qgcomp::qgcomp.noboot(telomean~.+I(ridageyr^2), expnms = furansla, data=dat[,c("telomean", mixturela, "ridageyr")]))
print(qgcftadj <- qgcomp::qgcomp.noboot(telomean~.+I(ridageyr^2), expnms = dioxinsla, data=dat[,c("telomean", mixturela, "ridageyr")]))


wqft <- gWQS::gwqs(telomean~wqs+ridageex+I(ridageex^2), rs=TRUE, n_vars=6, mix_name=mixture, data=dat[,c("telomean", mixture, "ridageex")])
summary(wqft)
head(wqft$final_weights)

wqft2<- gWQS::gwqs(telomean~wqs, b1_pos = FALSE, rs=TRUE, n_vars=6, mix_name=mixture, data=dat[,c("telomean", mixture)])
summary(wqft2)



