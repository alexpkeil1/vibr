library("vibr")
library('rstan')
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

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
X = vibr:::.scale_continuous(X)
apply(X, 2, sd)

stanmod <- "
data{
 int N;
 int P;
 real Y[N];
 matrix[N,P] X;
 real delta;
}
transformed data{}
parameters{
 real b0;
 vector[P] beta;
 real<lower=0> sigma;
}
transformed parameters{}
model{
 Y ~ normal(b0 + X * beta, sigma);
}
generated quantities{
 // variable importance
 real vimp[P];
 // probability of being most important
 real pimp[P];
 real ovimp;
 real pcbimp;
 real dioxinimp;
 real furanimp;
{
 real meanY = mean(Y);
 real invn = 1.0/N;
 matrix[N,P] Xt;
 real maxv;
 for (v in 1:P){
   vimp[v] = 0.0;
   Xt = X;
   Xt[,v] = X[,v] + delta;
   vimp[v] += mean(b0 + Xt * beta);
   vimp[v] -= mean(b0 + X * beta);
   vimp[v] *= 100.;
 }
 ovimp = sum(vimp[1:34]);
 pcbimp = sum(vimp[1:18]);
 dioxinimp = sum(vimp[19:25]);
 furanimp = sum(vimp[26:34]);
 maxv = max([max(vimp), -min(vimp)]);
 for (v in 1:P){
   pimp[v] = (fabs(vimp[v]) == maxv);
 }
}
}
"

standat <- list(
  X=X,
  Y=Y,
  N=length(Y),
  P=ncol(X),
  delta = 0.01
)

smod <- rstan::stan_model(model_code=stanmod)
sfit <- rstan::sampling(smod, data = standat)
sfit

print(sfit, pars = "pimp") # 35, 11, 6
mixturela[c(35, 11, 6)]
print(sfit, pars = "vimp") # 35, 11, 27,
mixturela[c(35, 11, 27)]
print(sfit, pars = c("ovimp", "pcbimp", "dioxinimp", "furanimp")) # 9, 27, 8


lmx = cbind(intercept=rep(1, length(standat$Y)), standat$X)
lmy = standat$Y

vi <- vibr::varimp(standat$X, standat$Y, delta=0.01,
             Y_learners = .default_continuous_learners()[1],
             Xdensity_learners = .default_density_learners(),
             Xbinary_learners = NULL,
             estimator = "TMLE"
             )
100*cbind(
   100*vi$res$est,
   glm.fit(x=lmx, y=lmy)$coefficients[-1],
   summary(sfit, pars=c("beta"))$summary[,'mean']
)
stan_dens(sfit,  c("ovimp", "pcbimp", "dioxinimp", "furanimp"))
stan_dens(sfit,  c("beta[9]"))
stan_dens(sfit,  c("vimp[9]"))
stan_dens(sfit,  c("beta[27]"))
stan_dens(sfit,  c("vimp[27]"))
