# estimator of Haneuse and Rotnitzky (2012)
# requires parametric models

# created out of laziness after I used the function, incorrectly, without defining it
.ones <- function(fun, var){
  rep(1.0, fun(var))
}

.rql <- function(){
}

.lambdaql <- function(){
}

.gamma_step_continuous <- function(Y, lambda_a, r_a){
  ft <- lm(Y ~ -1 + offset(r_a) + lambda_a)
  gamma <- ft$coefficients[1]
  gamma
}

.muq <- function(lambda_q, r_q, gamma, Il = .ones(length, r_q)){
  mui <- r_q + gamma * lambda_q
  mu_q <- sum(Il * mui) / sum(Il)
  mu_q
}

.delta_dr_U <- function(Y, mu_q, lambda_a, r_a, lambda_q, r_q, gamma, Il = .ones(length, r_q)){
  delta <- mean(Y) - mu_q
  muqi <- r_q + gamma * lambda_q
  mui <- r_a + gamma * lambda_a
  c1 <- Il * (Y - muqi) - delta
  c2 <- lambda_a*(Y - mui)
  U <- c1 - c2
  U
}

.delta_dr_S <- function(Y, mu_q, lambda_a, r_a, lambda_q, r_q, gamma, Il = .ones(length, r_q)){
  mui <- r_a + gamma * lambda_a
  c2 <- lambda_a*(Y - mui)
  St = stop("should be the estimating function sum(S(tau)) = 0, where tau is used to estimate parameters of outcom regression.")
  S = c(St, c2)
  S
}


.delta_dr <- function(Y, mu_q, lambda_a, r_a, lambda_q, r_q, gamma, Il = .ones(length, r_q)){
  Ela = mean(Il)
  U <- .delta_dr_U(Y, mu_q, lambda_a, r_a, lambda_q, r_q, gamma, Il)
}
