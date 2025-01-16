library(spatstat)
library(cubature)
library(data.table)
library(doParallel)
# library("fastHawkes")

# =====================
# Tests to do
# =====================
# for all methods = consistently recovers hawkes process parameters
# - optimized/not
# - exponential kernal or not (numeric + density estimate)
# speed check on optimization

PARAMS <- list(mu = 5,alpha = 1.5,beta = 10,K = .75)
TIME<- 200
OMEGA <- c(0,10^2,0,10^2)
N_CORES = 7
N = N_CORES*1

print(paste0("Expected points is ",PARAMS$mu*TIME*(1/(1-PARAMS$K))))

make_cluster <- function(N_CORES){
  # setup the cluster:
  cl <- makeCluster(N_CORES)
  registerDoParallel(cl)
  estimatedparams<-list()
  # export libraries to cluster:
  clusterEvalQ(cl, {
    library(data.table)
    library(dplyr)
    library(spatstat)
    library(cubature)
  })
  clusterExport(cl, c("PARAMS",
                      "TIME",
                      "OMEGA",
                      "N_CORES",
                      "N"))
  return(cl)
}
cl <- make_cluster(N_CORES)
registerDoParallel(cl)

print("Doing estimates 1")
# Fullt optimized expoenntial kernel
t1<-proc.time()[3]
estimates_1 <- foreach(i = 1:N, .packages = c("data.table", "dplyr", "spatstat",'cubature')) %dopar% {
  set.seed(i)
  onerealiz <- sim_hawkes(params = PARAMS,
                          windowT = c(0, TIME),
                          windowS = OMEGA)
  onerealiz <- as.data.frame(onerealiz)
  onerealiz <- onerealiz[-c(1,5)]
  guess <- list(mu = runif(1),
               alpha = runif(1),
               beta = runif(1),
               K = runif(1))
  params_fitted <- fit_hawkes(params_init = guess,
                              realiz = onerealiz,
                              windowT = c(0, TIME),
                              windowS = as.owin(OMEGA),
                              maxit = 1000,
                              exponential_kernel = TRUE,
                              kernel = NULL,
                              optimized = TRUE,
                              trace = 0,
                              poisson_flag = FALSE,
                              zero_background_region = NULL
                              )
  params_fitted
  return(params_fitted)
}
t1 <- proc.time()[3] - t1

print("Doing estimates 2")
# Non optimized exponential kernel
t2<-proc.time()[3]
estimates_2 <- foreach(i = 1:N, .packages = c("data.table", "dplyr", "spatstat",'cubature')) %dopar% {
  set.seed(i)
  onerealiz <- sim_hawkes(params = PARAMS,
                          windowT = c(0, TIME),
                          windowS = OMEGA)
  onerealiz <- as.data.frame(onerealiz)
  onerealiz <- onerealiz[-c(1,5)]
  guess <- list(mu = runif(1),
                alpha = runif(1),
                beta = runif(1),
                K = runif(1))
  params_fitted <- fit_hawkes(params_init = guess,
                              realiz = onerealiz,
                              windowT = c(0, TIME),
                              windowS = as.owin(OMEGA),
                              maxit = 1000,
                              exponential_kernel = TRUE,
                              kernel = NULL,
                              optimized = FALSE,
                              trace = 1,
                              poisson_flag = FALSE,
                              zero_background_region = NULL
  )
  params_fitted
  return(params_fitted)
}
t2 <- proc.time()[3] - t2

print("Doing estimates 3")
t3<-proc.time()[3]
# Non optimized non-exponential kernel (density integral)
estimates_3 <- foreach(i = 1:N, .packages = c("data.table", "dplyr", "spatstat",'cubature')) %dopar% {
  set.seed(i)
  onerealiz <- sim_hawkes(params = PARAMS,
                          windowT = c(0, TIME),
                          windowS = OMEGA)
  onerealiz <- as.data.frame(onerealiz)
  onerealiz <- onerealiz[-c(1,5)]
  guess <- list(mu = runif(1),
                alpha = runif(1),
                beta = runif(1),
                K = runif(1))
  params_fitted <- fit_hawkes(params_init = guess,
                              realiz = onerealiz,
                              windowT = c(0, TIME),
                              windowS = as.owin(OMEGA),
                              maxit = 1000,
                              exponential_kernel = FALSE,
                              kernel = function(t,x,beta,alpha){return(exp(-beta*t -alpha*x))},
                              optimized = FALSE,
                              trace = 0,
                              poisson_flag = FALSE,
                              zero_background_region = NULL,
                              density_approx = TRUE
  )
  params_fitted
  return(params_fitted)
}
t3 <- proc.time()[3] - t3

print("Doing estimates 4")
t4<-proc.time()[3]
# Non optimized non-exponential kernel (numeric integral)
estimates_4 <- foreach(i = 1:N, .packages = c("data.table", "dplyr", "spatstat",'cubature')) %dopar% {
  set.seed(i)
  onerealiz <- sim_hawkes(params = PARAMS,
                          windowT = c(0, TIME),
                          windowS = OMEGA)
  onerealiz <- as.data.frame(onerealiz)
  onerealiz <- onerealiz[-c(1,5)]
  guess <- list(mu = runif(1),
                alpha = runif(1),
                beta = runif(1),
                K = runif(1))
  params_fitted <- fit_hawkes(params_init = guess,
                              realiz = onerealiz,
                              windowT = c(0, TIME),
                              windowS = as.owin(OMEGA),
                              maxit = 1000,
                              exponential_kernel = FALSE,
                              kernel = function(t,x,alpha,beta){return(exp(-beta*t -alpha*x))},
                              optimized = TRUE,
                              trace = 1,
                              poisson_flag = FALSE,
                              zero_background_region = NULL,
                              density_approx = FALSE,
                              numeric_integral = TRUE
  )
  params_fitted
  return(params_fitted)
}
t4 <- proc.time()[3] - t4
stopCluster(cl)

# ==================================
# CHECK ESTIMATED PARAMS and speed
# ==================================
for(i in 1:3){
  params_df <- do.call(rbind, lapply(eval(parse(text = paste0("estimates_",i))),function(x) {
    return(as.data.frame(t(x$par)))
  }))
  print(paste0("Method ",i," took ",round(eval(parse(text = paste0("t",i))),2)," seconds"))
  print(paste0(" For method ",i," mean estimated params as follows: "))
  print(apply(params_df, 2, mean))
}


# =============================
# USe method 1 or a big PP to show consistency
# =============================
TIME <- 200
print(paste0("Expected points is ",PARAMS$mu*TIME*(1/(1-PARAMS$K))))
print("Doing estimates 1")
# Fullt optimized expoenntial kernel
t5<-proc.time()[3]
estimates_5 <- foreach(i = 1:N, .packages = c("data.table", "dplyr", "spatstat",'cubature')) %dopar% {
  set.seed(i)
  onerealiz <- sim_hawkes(params = PARAMS,
                          windowT = c(0, TIME),
                          windowS = OMEGA)
  onerealiz <- as.data.frame(onerealiz)
  onerealiz <- onerealiz[-c(1,5)]
  guess <- list(mu = runif(1),
                alpha = runif(1),
                beta = runif(1),
                K = runif(1))
  params_fitted <- fit_hawkes(params_init = guess,
                              realiz = onerealiz,
                              windowT = c(0, TIME),
                              windowS = as.owin(OMEGA),
                              maxit = 1000,
                              exponential_kernel = TRUE,
                              kernel = NULL,
                              optimized = TRUE,
                              trace = 0,
                              poisson_flag = FALSE,
                              zero_background_region = NULL
  )
  params_fitted
  return(params_fitted)
}
t5 <- proc.time()[3] - t5
# ~5 minuts for 4000 point realization

i = 5
params_df <- do.call(rbind, lapply(eval(parse(text = paste0("estimates_",i))),function(x) {
  return(as.data.frame(t(x$par)))
}))
print(paste0("Method ",i," took ",round(eval(parse(text = paste0("t",i))),2)," seconds"))
print(paste0(" For method ",i," mean estimated params as follows: "))
print(apply(params_df, 2, mean))


