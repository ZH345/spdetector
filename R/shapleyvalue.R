library(readr)
library(lme4)
library(dplyr)
library(nlme)
library(merTools)
library(data.table)

rm(list = ls())
VarG <- list(VarG1 = c("dem", "sdem", "slope", "tempe", "npp", "pcrop", "preci"),
           VarG2 = c("river"),
           VarG3 = c("applyfor"),
           VarG4 = c("landarea", "road"))

data1995 <- read_csv("D:/0 job/0 q-vpc/data/yrdata1995.csv", col_names = T)
data1995$applyfor <- log(data1995$applyfor + 1)
#data2015 <- read_csv("yrdata2015.csv", col_names = T)
#data2015$applyfor <- log(data2015$applyfor + 1)

res1995 <- shapley_value(VarG, data = data1995)

#res2015 <- shapley_value(VarG, data = data2015)





### customized Shapley's value function
library(gtools)
shapley_value <- function(
                      VarGroups, ## a named list defining variable groups
                      data       ## data to fit models
                      ) {
  
  ## extract variable groups and formula
Varfull <- VarGroups 
J <- length(Varfull)
  
  ##generate possible combinations as input profiles
  index<-eval(parse(text = paste0("expand.grid(",paste(rep("c(0,1)",J), collapse=", "),")")))
  
  ##transform the combinations to the variable names included 
  namesindex<-apply(index,1,function(x) unlist(Varfull[which(x==1)]))
  
  ##the model list
  modellist<-lapply(namesindex, function(x) 
    paste("model <- lmer(ntlvalue ~ ",
          paste(x, collapse = " + "), ##paste variables entered
          " + (1|ProAdCode) + (1|ProAdCode:CityAdCode) + (1|ProAdCode:CityAdCode:AdminCode), data=data)" 
    )
  )
  
  ##function to compute R2
  computeR2<-function(model){
    modelfit<-eval(parse(text = model))
    ### Calculated the variance of y explained by variables
    fitted_y <- predict(modelfit, re.form = NA)
    var_y_f <- var(fitted_y)
    ### Extract variances of random effects 
    temp.var <- as.data.frame(VarCorr(modelfit))$vcov
    # Note the data has a four-level structure so there would be 4 variances
    R2 <- (var_y_f ) / (var_y_f + sum(temp.var))
    
    return(R2)
  }
  
  ##compute 2^p R2
  R2<-sapply(modellist,computeR2)
  
  S<-c()
  for (k in 1:J){ ## equation (4)
    
    xindex<-index[apply(index,1,function(x) x[k]==0),] ##excluding k
    xkindex<-xindex 
    xkindex[,k]<-1  ##including k
    
    gxindex<-as.numeric(rownames(xindex))
    gxkindex<-(1:2^J)[-gxindex]
    
    #alternative
    #gxindex<-which(apply(index,1,function(x) any(apply(xindex,1,function(y) all(x==y)))))
    #gxkindex<-which(apply(index,1,function(x) any(apply(xkindex,1,function(y) all(x==y)))))
    
    R2increase<-R2[gxkindex]-R2[gxindex]
    
    xbar<-apply(xindex,1,sum)
    commonpart<-factorial(xbar)*factorial(J-xbar-1)/factorial(J)
    
    S[k]<-sum(R2increase*commonpart)
    
  }
  
  return(S)
  
}
