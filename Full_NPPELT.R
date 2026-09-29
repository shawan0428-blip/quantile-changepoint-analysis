rm(list = ls())
dev.off()
library(quantreg)
library(EnvCpt)



#getwd()
#setwd('/Users/shawan/Desktop/ChangePoint')
#data <- read.table("Tuscaloosa_annual.txt")
#data <- data[-1]



#############################################
# 1. Cost Function for a One Segment 
#############################################
find_tk <- function(data) {
  n   <- length(data)
  K0  <- ceiling(4 * log(n))
  c0  <- -log(2*n - 1)
  gamma <- c0 / K0
  
  # 1) compute all p_k
  ks  <- seq_len(K0)
  p_k <- 1 / (1 + (2*n - 1) * exp(gamma * (2*ks - 1)))
  
  # 2) keep only those above 25/n and less than (1-25/n)
  low=25/n
  sel <- p_k > low & p_k < (1-low)
  p_k   <- p_k[sel]
  K_new <- length(p_k)
  
  # 3) get the matching quantiles
  t_k <- quantile(data, probs = p_k)
  
  
  #list(K = K_new, p = p_k, t = t_k)
  return(t_k)
}

############################################
mll_nonparametric_ed <- function(segdata, t_k) {
  n <- length(segdata)
  if (n < 2) return(0)  # not enough data for lag-1 regression
  K <- length(t_k) # K is new K
  # Create the lag-1 relationship.
  y_current <- segdata[-1]
  y_lag1    <- segdata[-n]
  
  # Build a data frame for regression.
  dataset <- data.frame(y = y_current, y_lag1 = y_lag1)
  
  taus <- ecdf(segdata)(t_k)
  
  # Fit QAR models at each tau level.
  qar_models <- lapply(taus, function(tau_val) {
    rq(y ~ y_lag1, tau = tau_val, data = dataset)
  })
  
  # Combine residuals from all QAR models into a matrix.
  residuals_matrix <- do.call(cbind, lapply(qar_models, function(m) m$residuals))
  
  #Build Q matrix (empirical CDF function):  
  Q <- sapply(seq_along(t_k), function(i) {
    cumsum(residuals_matrix[, i] < 0) +
      0.5 * cumsum(residuals_matrix[, i] == 0)
  })
  
  # Ensure Q is a matrix with the correct number of rows.
  Q <- matrix(Q, nrow = nrow(residuals_matrix))
  
  # Effective number of observations (due to lag)
  
  if (n-1 <= 0) return(0)
  
  cost <- 0
  
  # For each quantile level, compute the segment’s cost
  for (i in 1:K) {
    # For a segment starting at the beginning (n-1 due to lag)
    # last row of Q (cumsum)
    sumstat <- Q[n-1, i]
    if (is.na(sumstat)) next
    # Empirical CDF for this segment.
    Fkl <- sumstat / (n-1)
    if (is.na(Fkl) || Fkl <= 0 || Fkl >= 1) next
    cost <- cost + (n-1)* (Fkl * log(Fkl) + (1 - Fkl) * log(1 - Fkl))
  }
  
  # Scale the cost
  -2 * log(2 * n - 1) * cost / K
  #return(list(cost = cost, Q = Q))
}


#############################################
# 2. PELT Algorithm
#############################################
pelt_nonparametric <- function(data, pen = 0, minseglen) {
  t_k <- find_tk(data)
  n <- length(data)
  K   <- length(t_k)
  
  lastchangecpts <- rep(0, n + 1)
  lastchangelike <- -pen
  checklist <- NULL
  
  for (i in minseglen:(2 * minseglen - 1)) {
    lastchangelike[i + 1] <- mll_nonparametric_ed(data[1:i], t_k)
    lastchangecpts[i + 1] <- 0
  }
  
  # Initialize candidate change-point checklist.
  checklist <- c(0, minseglen)
  
  for(tstar in (2*minseglen):n){
    tmplike=unlist(lapply(checklist,FUN=function(tmpt){
      return(lastchangelike[tmpt+1]+mll_nonparametric_ed(data[(tmpt + 1):tstar], t_k)+pen)}))
    lastchangelike[tstar+1]=min(tmplike,na.rm=TRUE)
    lastchangecpts[tstar+1]=checklist[which.min(tmplike)[1]]
    checklist=checklist[tmplike<=(lastchangelike[tstar+1]+pen)]
    checklist=c(checklist,tstar-minseglen+1)
    #if(verbose){if(tstar%%10==0){print(paste("Finished",tstar))}}
  }
  fcpt=NULL
  last=n
  while(last!=0){
    fcpt=c(fcpt,lastchangecpts[last+1])
    last=lastchangecpts[last+1]
  }
  return(sort(fcpt)[-1])
}

#############################################
# 3. a1b1 
#############################################
library(changepoint)
GDPQAR=function(n=500,a1=.1,b1=.8,a2=.11,b2=0.2332299){
  
  #GET STARTED
  
  upast=runif(10000)#Generates 10,000 uniform random numbers between 0 and 1 
  
  utpast=qnorm(upast)#Converts them to a standard normal distribution using the quantile function
  
  alphatpast=pmin(a1+b1*upast,1) #  \alpha_t = a1 + b1 U_t 
  
  y0=(utpast[1]+sum(utpast[2:10000]*cumprod(alphatpast[-10000])))
  # y_0 = u_1 + \sum_{t=2}^{10000} \left( u_t \times \prod_{j=1}^{t-1} \alpha_j \right)
  
  
  u=runif(n) #Generates a new sequence of random number
  
  y=u
  
  ut=qnorm(u)
  
  cpt=floor(n/2)#cpt is the change point, where the structural break happens (halfway through the series).
  
  alphat=1:n #autoregressive coefficient 
  
  alphat[1:cpt]=pmin(a1+b1*u[1:cpt],1) #Assigns first-half persistence
  
  y[1]=alphat[1]*y0+ut[1] #Initializes the first value
  
  for(i in 2:cpt)
    
    y[i]=alphat[i]*y[(i-1)]+ut[i] #Each value depends on the previous value, weighted by  \alpha_t 
  
  #START SHIFT in STATIONARY DISTRIBUTION
  
  upast=runif(10000)
  
  utpast=qnorm(upast)
  
  alphatpast=pmin(a2+b2*upast,1)#The Structural Shift (New Stationary Process), parameters  a2, b2 .
  
  ystart=(utpast[1]+sum(utpast[2:10000]*cumprod(alphatpast[-10000])))
  
  alphat[(cpt+1):n]=pmin(a2+b2*u[(cpt+1):n],1)#Creates a new starting value after the structural break.
  
  y[cpt]=alphat[cpt]*ystart+ut[cpt]
  
  for(i in (cpt+1):n)#Creates a new starting value after the structural break.
    
    y[i]=alphat[i]*y[(i-1)]+ut[i]
  
  return(y)                             
  
}
#############################################
# 4. Run PELT & Plot Results
#############################################


set.seed(123)
y <- GDPQAR(n=500,a1=0.1,b1=0.8,a2=0.5,b2=0.03614288)
#y<-data

n=length(y)

tmp <- find_tk(y)
K   <- length(tmp)
penalty_val <-(2*K+1) * log(n)

##Quick sanity check
print(tmp)
print(K)


cpts_idx <- pelt_nonparametric(y, pen = penalty_val, minseglen = 20)
print(cpts_idx)


# Plot the series
plot(y, type="l", lwd=2)

# Add vertical line for structural break
abline(v=cpts_idx, col="red", lty=2, lwd=2)




#############################################
# envcpt
#############################################

fit_env <- envcpt(y)
fit_env
# 2) Pick “best” model by BIC 
best_by_bic <- names(which.min(BIC(fit_env)))
best_fit    <- fit_env[[best_by_bic]]
best_by_bic
best_fit 
cp_envcpt <- if (inherits(best_fit, "cpt")) cpts(best_fit) else integer(0)
cp_envcpt
# compare to your custom detector:
cpts_idx


