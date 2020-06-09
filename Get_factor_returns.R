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
Rtn_fac=data.table(factor.returns)
Rtn_fac[,Date:=rownames(factor.returns)]
save(Rtn_fac, file = "Rtn_fac.rda")