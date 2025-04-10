# Implements Lasso for logistic regression, both classical/bayesian ways

## 1. SIMULATION

# for reproducibility
set.seed(666)

# sample size
n <- 100

# explanatory variables
x1 <- rnorm(n)
x2 <- rnorm(n)
x3 <- rnorm(n)

# regression parameters (X1 and X3 are significant, X2 is not)
beta <- c(0.5, 1, 0, 1) # beta0, beta1, beta2, beta3

# prob of success
logit_prob <- beta[1] + beta[2] * x1 + beta[3] * x2 + beta[4] * x3
prob <- plogis(logit_prob)

# response var
y <- rbinom(n, 1, prob)

## 2. AJUSTEMENT

# classical framework, no lasso
summary(glm(y~x1+x2+x3, family = "binomial"))

# classical framework, lasso
library(glmnet)
x.mat <- cbind(x1,x2,x3)
cv.fit <- cv.glmnet(x.mat, y, family = "binomial", type.measure = "class")
plot(cv.fit)
cv.fit$lambda.min
coef(cv.fit, s = "lambda.min")
#fit <- glmnet(x.mat, y)
#plot(fit)

# bayesian framework, no lasso
library(nimble)

# model
model <- nimbleCode({
  for (i in 1:n){
    logit(p[i]) <- beta[1] + beta[2] * X[i,1] + beta[3] * X[i,2] + beta[4] * X[i,3]
    y[i] ~ dbern(p[i])
  }
  for(j in 1:4){
    beta[j] ~ dnorm(0,0.1) # jags uses precision = 1 / var
  }
}
)

# data
data <- list(y = y , X = cbind(x1, x2, x3), n = n)

# initial values (2 chains)
init <- list(list(beta = rnorm(4)),
             list(beta = rnorm(4)))

# fit model
out <- nimbleMCMC(data = data,
            inits = init,
            monitors = c("beta"),
            code = model,
            nchains = 2,
            niter = 6000,
            nburnin = 1000,
            thin = 1)

# print results
MCMCsummary(out)

MCMCtrace(out, pdf = FALSE, ind = TRUE, params = "beta")
MCMCplot(out, params = "beta")



# bayesian framework, lasso

# model (https://stats.stackexchange.com/questions/28609/regularized-bayesian-logistic-regression-in-jags)
model <- nimbleCode({
  for (i in 1:n){
    logit(p[i]) <- beta[1] + beta[2] * X[i,1] + beta[3] * X[i,2] + beta[4] * X[i,3]
    y[i] ~ dbern(p[i])
  }
  beta[1] ~ dnorm(0, sd = 10) 
  # L1 regularization == a Laplace (double exponential) prior 
  for (j in 2:4) {
    beta[j] ~ ddexp(0, lambda)  
  }
  lambda ~ dunif(0.001,10)
})

# data
data <- list(y = y , X = cbind(x1, x2, x3), n = n)

# initial values (2 chains)
init <- list(list(beta = rnorm(4)),
             list(beta = rnorm(4)))

# fit model
out <- nimbleMCMC(data = data,
            inits = init,
            monitors = c("beta", "lambda"),
            code = model,
            nchains = 2,
            niter = 6000,
            nburnin = 1000,
            thin = 1)

# print results
MCMCsummary(out)

# check convergence
MCMCtrace(out, pdf = FALSE, ind = TRUE, params = "beta")
MCMCplot(out, params = "beta")


# bayesian framework, lasso (avec sélection de variables ?!)

# model (https://stats.stackexchange.com/questions/28609/regularized-bayesian-logistic-regression-in-jags)
model <- nimbleCode({
  for (i in 1:n){
    logit(p[i]) <- beta[1] + betagamma[1] * X[i,1] + betagamma[2] * X[i,2] + betagamma[3] * X[i,3]
    y[i] ~ dbern(p[i])
  }
  beta[1] ~ dnorm(0,0.1) # jags uses precision = 1 / var
  # L1 regularization == a Laplace (double exponential) prior 
  for (j in 2:4) {
    beta[j] ~ ddexp(0, lambda)  
  }
  lambda ~ dunif(0.001,10)
  for (j in 1:3) {
    betagamma[j] <- beta[j+1] * gamma[j]
    gamma[j] ~ dbern(0.5)
  }
})

# data
data <- list(y = y , X = cbind(x1, x2, x3), n = n)

# initial values (2 chains)
init <- list(list(beta = rnorm(4), gamma = rbernoulli(3, 0.5)),
             list(beta = rnorm(4), gamma = rbernoulli(3, 0.5)))

# fit model
out <- nimbleMCMC(data = data,
                  inits = init,
                  monitors = c("beta", "lambda", "betagamma", "gamma"),
                  code = model,
                  nchains = 2,
                  niter = 6000,
                  nburnin = 1000,
                  thin = 1)

# print results
MCMCsummary(out)

# check convergence
MCMCtrace(out, pdf = FALSE, ind = TRUE, params = "beta")
MCMCplot(out, params = "beta")

# for Bayesian Lasso, see
# https://people.eecs.berkeley.edu/~jordan/courses/260-spring09/other-readings/park-casella.pdf
# http://www2.stat-athens.aueb.gr/~jbn/papers2/23b_Lykou_Ntzoufras_2011_Wires_WinBUGS_final.pdf