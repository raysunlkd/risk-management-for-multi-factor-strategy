setwd('F:\\Courses\\HKUST_Semester_Spring\\MAFS5210\\project2\\data_for_R')

# ERASE DATA AND CALL GARBAGE COLLECTOR
rm(list=ls())
gc()

# LOAD PACKAGES
list.of.packages <- c("data.table", "fasttime",'plyr',"PerformanceAnalytics",
                      'imputeTS',"parallel","doParallel","doMC",'lubridate',
                      'anytime','xts','TTR','missRanger','mvtnorm','MASS')
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages, require, character.only = TRUE)

# LOAD ENVIRONMENT
if(Sys.info()['sysname'] == "Windows" ) {
  library(doParallel)
  registerDoParallel(cores=detectCores()-1)
} else {
  library(doMC)
  registerDoMC(cores=detectCores()-1)
}
setDTthreads(detectCores()*2-2)


# LOAD DATA
Data= fread('csi_prices_and_industry.csv',header = TRUE,encoding = "UTF-8")
Data=  Data[order(index),]
Data[,RETURN.FUTURE:=c(diff((log(1+close))),rep(0,1)),by=symbol]
#Data[,size:= log(1+size)]

exposure.vars = c('Sector','MOM','vol_volatility_1m','volatility_1m',
                  'pe','circ_mv','dv_ratio',
                  'shibor_1m', 'shibor_6m', 'libor_1m', 'libor_6m','shibor_lpr_1y')
ret.var = 'RETURN.FUTURE'

# TS standardization
zScore <- function(i) {
  meanExp = mean(i,na.rm=T)
  sigmaExp = sd(i,na.rm=T)
  sigmaEWMA <- stdExpo <- exposures <- i
  ts <- (i - meanExp)^2
  var_past_2 <- sigmaExp ^ 2
  sigmaEWMA <- sapply(ts, function(x) var_past_2 <<- 0.10 * x + 0.90 * var_past_2)
  sigmaEWMA[which(sigmaEWMA==0)]<- 1
  as.vector((i -  meanExp) / sqrt(sigmaEWMA))
}
for(i in exposure.vars[-which(exposure.vars %in% c('Sector'))]) Data[,(i):=zScore(get(i)) ,by=symbol]

### DATA PREPROCESSING
# SUBSET DATA OR IT WILL CRASH FOR MISSRANGER
# Data <- missRanger(Data, pmm.k = 3, num.trees = 100)
for (j in 1:ncol(Data)) set(Data, which(is.infinite(Data[[j]])|is.nan(Data[[j]])), j, NA)
Data = Data[complete.cases(Data)]
Data = Data[-which(index %in% sort(unique(Data$index))[length(unique(Data$index))])]
which.numeric <- sapply(Data[,exposure.vars,with=FALSE], is.numeric)
exposures.num <- exposure.vars[which.numeric]
exposures.char <- exposure.vars[!which.numeric]
formula.expochar = as.formula(paste(ret.var, "~", exposures.char, "-1"))
factor.names <- c("Market",paste(levels(Data[,exposures.char,with=F]),sep=" "), exposures.num)
beta.expochar <- model.matrix(formula.expochar, data=Data)
rownames(beta.expochar) <- Data$symbol

# Beta for the whole model (generally without intercept)
fm.formula <- as.formula(paste(ret.var, "~", paste(exposure.vars, collapse="+")))
beta <- model.matrix(fm.formula, data=Data)
rownames(beta) <- Data$symbol

#Define beta.star as Beta of the whole model with Intercept/Market represtend by col of ones
beta.star <- cbind("Market" = rep(1, nrow(beta.expochar)), beta.expochar)

beta.style<- matrix(beta[,exposures.num], ncol = length(exposures.num))
colnames(beta.style) = exposures.num

#DEFINE RESTRICTION MATRIX 
K <- dim(beta.star)[2]
R_matrix = rbind(diag(K-1), c(0,rep(-1,K-2)))
uni = beta.star[unique(row.names(beta.star)),]
B.mod = (as.matrix(uni)) %*% R_matrix

# GET FORMULA TO REGRESS
fmSI.formula = paste(ret.var, "~", "B.mod+", paste(exposures.num, collapse="+"),"-1" )

# CONCATENATE ALL REGS
reg = function(data) {
  #data = Data[which(index==unique(Data$index)[2])]
  data = data[which(RETURN.FUTURE!=0)]
  data = unique(data)
  B.mod.tmp = B.mod[which(row.names(B.mod) %in% data$symbol),]
  fmSI.formula.tmp = as.formula(gsub('B.mod','B.mod.tmp',fmSI.formula))
  reg.time = lm(fmSI.formula.tmp,data=data)
  res = list(residuals(reg.time))
  names(res[[1]]) = data$symbol
  list(residuals=res,
       beta=list(coefficients(reg.time)),
       r.squared=summary(reg.time)$r.squared)
}
reg.list = Data[,reg(.SD),by=index]

# GET FACTOR RETURNS
g = sapply(reg.list$beta,c)
for(i in 1:ncol(g)) g[which(is.na(g[,i])),i] = 0
factor.returns  = R_matrix %*% g[1:(K-1), ]
factor.returns = t(rbind(factor.returns, g[K:nrow(g), ]))
uni = unique(Data$Sector)
colnames(factor.returns)[1:(length(uni)+1)] = c('Market',uni)
row.names(factor.returns) = unique(reg.list$index)
save(factor.returns, file = "Rtn_fac.rda")

# CUMULATIVE PERFORMANCE RETURNS
cumPerf = apply(factor.returns,2,cumsum)
cumPerf = data.table(cumPerf)
cumPerf$time = fastPOSIXct(sort(unique(Data$index)))
cumPerf = xts(cumPerf[,-c('time'),with=F],order.by=fastPOSIXct(cumPerf$time))
print(tail(cumPerf,1))

# GET RESIDUALS
res = do.call(rbind,apply(reg.list,1,function(x) cbind(x$index,names(x$residuals),x$residuals)))
res = data.table(res)
residuals = dcast(res,V1~ V2)
resid.var <- apply(coredata(residuals), 2, var, na.rm=T)
names(resid.var) = names(residuals)
resid.cov <- diag(resid.var)
colnames(resid.cov) = names(residuals)

# GET COVARIANCE MATRIX
beta.combine = cbind(beta.star, beta.style)
factor.cov <- cov(coredata(factor.returns))
s = Data[which(index==index[length(index)]),]
symbol = unique(s$symbol)
beta = beta.combine[(nrow(beta.combine) - length(symbol) + 1):nrow(beta.combine), ]
factor.cov <- cov(coredata(factor.returns))
return.cov <-  beta %*% factor.cov %*% t(beta) 
colnames(resid.cov)[which(!colnames(resid.cov)%in%colnames(return.cov))]
return.cov = return.cov[-which(!colnames(return.cov)%in%colnames(resid.cov)),-which(!colnames(return.cov)%in%colnames(resid.cov))]
resid.cov = resid.cov[-which(!colnames(resid.cov)%in%colnames(return.cov)),-which(!colnames(resid.cov)%in%colnames(return.cov))]
return.cov <-return.cov  + resid.cov



# FFM RISK DECOMPOSITION: fmEsDecomp.R
beta.star = beta
resid.var2 = resid.var[-which(!row.names(beta) %in% names(resid.var))]
beta.star <- as.matrix(cbind(beta, sqrt(resid.var)))
colnames(beta.star)[dim(beta.star)[2]] <- "Residuals"
K <- ncol(beta)
factor.star.cov <- diag(K+1)
factor.star.cov[1:K, 1:K] <- factor.cov
colnames(factor.star.cov) <- c(colnames(factor.cov),"Residuals")
rownames(factor.star.cov) <- c(colnames(factor.cov),"Residuals")


MU <- c(colMeans(factor.returns, na.rm=TRUE), 0)
SIGB <-  beta.star %*% factor.star.cov
N <- length(row.names(beta))
K = ncol(factor.returns)
ES.fm <- rep(NA, N)
mES <- matrix(NA, N, K+1)
cES <- matrix(NA, N, K+1)
pcES <- matrix(NA, N, K+1)
rownames(mES)=rownames(cES)=rownames(pcES)=row.names(beta)
colnames(mES)=colnames(cES)=colnames(pcES)=c(colnames(factor.returns),"Residuals")

p = 0.05
for (i in sort(row.names(beta))) {
  #i =  sort(row.names(beta))[1]
  print(i)
  R.xts <- Data[which(symbol==i),]$RETURN.FUTURE
  beta.i <- beta.star[i,,drop=F]
  ES.fm[i] <- -(beta.star[i,] %*% MU + sqrt(beta.i %*% factor.star.cov %*% t(beta.i))*dnorm(qnorm(p))/(p)) 
  # COMPUTE MARGINAL EXPECTED SHORTFALL
  mES[i,] <- -(t(MU) + SIGB[i,]/sd(R.xts, na.rm=TRUE) * dnorm(qnorm(p))/(p))
  cf <- as.numeric( ES.fm[i] / sum(mES[i,]*beta.star[i,], na.rm=TRUE) )
  
  # compute marginal, component and percentage contributions to ES
  # each of these have dimensions: N x (K+1)
  mES[i,] <- cf * mES[i,]
  cES[i,] <- mES[i,] * beta.star[i,]
  pcES[i,] <- 100* cES[i,] / ES.fm[i]
}

