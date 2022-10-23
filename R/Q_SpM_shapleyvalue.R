library(readr)
library(lme4)
library(dplyr)
library(nlme)
library(merTools)
library(data.table)

#library(tmap)
library(geodetector)
library(spatialreg)
library(spdep)

data(DiseaseData_shp)
data(SoilType_shp)
data(Watershed_shp)
data(Elevation_shp)



data(CollectData)
nrow(geodetector::CollectData)
DiseaseData_shp$SP_ID = as.numeric(DiseaseData_shp$SP_ID)
CollectData2<-maps2dataframe(DiseaseData_shp,c(SoilType_shp, Watershed_shp, Elevation_shp),namescolomn= c('SP_ID',
                                                                            'soiltype', 'watershed', 'elevation'))
nrow(CollectData2)

CollectData2$SP_ID = as.numeric(CollectData2$SP_ID)
DiseaseData_sf = as(DiseaseData_shp,'sf')
DiseaseData_sf = DiseaseData_sf[DiseaseData_sf$SP_ID %in% CollectData2$SP_ID,]

DiseaseData_sf = DiseaseData_sf[order(DiseaseData_sf$SP_ID),]
CollectData2 = CollectData2[order(CollectData2$SP_ID),]
DiseaseData_sf$SP_ID - CollectData2$SP_ID
CollectData2 = cbind(incidence = DiseaseData_sf$incidence,CollectData2)
nrow(DiseaseData_sf)
#DiseaseData_sf = merge(DiseaseData_sf,CollectData2,by.x = 'SP_ID',by.y ='SP_ID')


#tm4<-tm_shape(DiseaseData_sf)+tm_polygons(col = 'incidence',breaks=c(5.66,6.13,6.38,6.58),palette="RdYlGn")
#tm4

## factor_detector
factor_detector('incidence','elevation',CollectData)
factor_detector('incidence','elevation',CollectData2)
x.mat = model.matrix(lm(incidence~ factor(elevation),data = CollectData2))
sar.data = as.data.frame(cbind(y = CollectData2$incidence, x.mat))
nb_DiseaseData_sf <- spdep::poly2nb(DiseaseData_sf)
lisw_DiseaseData_sf <- spdep::nb2listw(nb_DiseaseData_sf, style = "W", zero.policy = TRUE)
mat.w = listw2mat(lisw_DiseaseData_sf)
sar.model <- lagsarlm(y ~x.mat[,-1], data = sar.data,
                      listw = lisw_DiseaseData_sf, zero.policy = T)

mi = moran.test(as.vector(sar.data$y), lisw_DiseaseData_sf)
mi.num= as.numeric(unlist(mi$estimate[1]))

var(sar.model$rho * mat.w %*% as.vector(sar.data$y))/var(sar.data$y)



summary(sar.model)
inv.rhoW = invIrM(neighbours = nb_DiseaseData_sf,rho = as.numeric(sar.model$rho))
incidence_hat = predict(sar.model,listw = lisw_DiseaseData_sf, pred.type = 'TS')
incidence_hat = predict(sar.model,listw = lisw_DiseaseData_sf, pred.type = 'TC')
factor_detector('incidence','elevation',CollectData)
R.squared = 1-(var(CollectData2$incidence - incidence_hat))/var(CollectData2$incidence)
R.squared

## interaction_detector

### customized Shapley's value function
library(gtools)
shapley_value.sar <- function(
    VarGroups, ## a named list defining variable groups
    data,       ## data to fit models
    listW
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
    paste0("model <- lagsarlm(y ~ ",
           paste(x, collapse = " + "),",data= data, listw = listW, zero.policy = T)" ##paste variables entered
    )
  )
  modellist[[1]] = "model <- lagsarlm(y ~ 1, data= data, listw = listW, zero.policy = T)"
  #modellist[[1]] = "model <- lm(y ~ 1 ,data= data)"
  ##function to compute R2
  computeR2<-function(model){
    modelfit<-eval(parse(text = model))
    ### Calculated the variance of y explained by variables
    #fitted_y <- predict(modelfit, re.form = NA)
    fitted_y <- predict(modelfit, pred.type = 'TC',listw = listW, re.form = NA)
    #var_y_f <- var(fitted_y)
    ### Extract variances of random effects
    #temp.var <- as.data.frame(VarCorr(modelfit))$vcov
    # Note the data has a four-level structure so there would be 4 variances
    R2 <- 1 - var(data$y - fitted_y)/var(data$y)

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
    print(k)
  }

  return(S)

}

shapley_value.ols <- function(
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
    paste0("model <- lm(y ~  ",
           paste(x, collapse = " + "),",data= data)" ##paste variables entered
    )
  )
  #modellist[[1]] = "model <- lagsarlm(y ~ 1, data= data, listw = listW, zero.policy = T)"
  modellist[[1]] = "model <- lm(y ~ 1 ,data= data)"
  ##function to compute R2
  computeR2<-function(model){
    modelfit<-eval(parse(text = model))
    ### Calculated the variance of y explained by variables
    #fitted_y <- predict(modelfit, re.form = NA)
    fitted_y <- predict(modelfit,re.form = NA)
    #var_y_f <- var(fitted_y)
    ### Extract variances of random effects
    #temp.var <- as.data.frame(VarCorr(modelfit))$vcov
    # Note the data has a four-level structure so there would be 4 variances
    R2 <- 1 - var(data$y - fitted_y)/var(data$y)

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
    print(k)
  }

  return(S)

}



###Double
interaction_detector('incidence',c('elevation','soiltype'),CollectData)

CollectData3 = CollectData2
CollectData3 = CollectData3[,-2]
names(CollectData3) = c('y',  'x1', 'x2', 'x3')
x.mat = model.matrix(lm(y~ factor(x1)*factor(x2),data = CollectData3))

sar.data = as.data.frame(cbind(y = CollectData3$y, x.mat[,-1]))
nb_DiseaseData_sf <- spdep::poly2nb(DiseaseData_sf)
lisw_DiseaseData_sf <- spdep::nb2listw(nb_DiseaseData_sf, style = "W", zero.policy = TRUE)

sar.model <- lagsarlm(sar.data$y ~x.mat[,-1],
                      listw = lisw_DiseaseData_sf, zero.policy = T)

summary(sar.model)
inv.rhoW = invIrM(neighbours = nb_DiseaseData_sf,rho = as.numeric(sar.model$rho))
incidence_hat = predict(sar.model,listw = lisw_DiseaseData_sf, pred.type = 'TS')
incidence_hat = predict(sar.model,listw = lisw_DiseaseData_sf, pred.type = 'TC')

R.squared = 1-(var(CollectData2$incidence - incidence_hat))/var(CollectData2$incidence)
R.squared


library(dplyr)

library(stringr)
names.data = names(sar.data)

names.data = str_replace_all(names.data,'[()]','.')
names(sar.data) = names.data

names.data = names.data[-1]
VarG.list = NULL


BitMatrix <- function(n){
  # 作用：返还Bit矩阵
  # Args:n：数据长度
  set <- 0:(2^n-1)
  rst <- matrix(0,ncol = n,nrow = 2^n)
  for (i in 1:n){
    rst[, i] = ifelse((set-rowSums(rst*rep(c(2^((n-1):0)), each=2^n)))/(2^(n-i))>=1, 1, 0)
  }
  rst
}

k <- c('x1','x2')

names(CollectData2)[3:5]
var_name = c("soiltype", "elevation")
selectMatrix <- BitMatrix(length(k))
Varg.name = apply(selectMatrix, 1, function(x){ var_name[which(x==1)] })
Varg.name = str_replace_all(Varg.name,'"','')


#Construct VarG based on selectMatrix
for (ri in 2:nrow(selectMatrix)){
  VarG.temp = names.data
  for (ci in 1:length(selectMatrix[ri,])){
    var_name = k[ci]
    if (selectMatrix[ri,ci] == 1){
      VarG.temp = VarG.temp[str_detect(string = VarG.temp ,pattern = var_name)]
    }else{
      VarG.temp = VarG.temp[str_detect(string = VarG.temp ,pattern = var_name,negate = T)]
    }
  }
  VarG.list = append(VarG.list,list(VarG.temp))
}


VarGroups = VarG.list
data = sar.data
listW = lisw_DiseaseData_sf
ols.shapley = shapley_value.ols(VarGroups = VarG.list, data = sar.data)
rbind(Varg.name[-1],ols.shapley)
interaction_detector('incidence',c('elevation','soiltype'),CollectData)

sar.shapley = shapley_value.sar(VarGroups = VarG.list, data = sar.data,listW  = lisw_DiseaseData_sf)
Varg.name[-1]
sar.shapley
rbind(Varg.name[-1],sar.shapley)









###Multi (3)
interaction_detector('incidence',c('elevation','soiltype','watershed'),CollectData)

CollectData3 = CollectData2
names(CollectData2)
CollectData3 = CollectData3[,-2]
names(CollectData3) = c('y',  'x1', 'x2', 'x3')
x.mat = model.matrix(lm(y~ factor(x1)*factor(x2)*factor(x3),data = CollectData3))

sar.data = as.data.frame(cbind(y = CollectData3$y, x.mat[,-1]))
nb_DiseaseData_sf <- spdep::poly2nb(DiseaseData_sf)
if(F){
  DiseaseData_sf_plot = merge(DiseaseData_sf,CollectData2,by.x = 'SP_ID', by.y = 'SP_ID')
  DiseaseData_sf_plot$plot.id = paste0(as.character(DiseaseData_sf_plot$soiltype),as.character(DiseaseData_sf_plot$watershed),as.character(DiseaseData_sf_plot$elevation))
  #DiseaseData_sf_plot$plot.id = as.numeric(DiseaseData_sf_plot$plot.id)
  ggplot()+
    geom_sf(data = DiseaseData_sf_plot,aes(fill = plot.id))
  ggplot()+
    geom_sf(data = DiseaseData_sf_plot,aes(fill =incidence.x ))+
    scale_fill_continuous(type = "viridis")

  tm4<-tm_shape(DiseaseData_sf_plot)+tm_polygons(col = 'plot.id',palette="RdYlGn",)
  tm4

}
lisw_DiseaseData_sf <- spdep::nb2listw(nb_DiseaseData_sf, style = "W", zero.policy = TRUE)

sar.model <- lagsarlm(sar.data$y ~x.mat[,-1],
                      listw = lisw_DiseaseData_sf, zero.policy = T)
summary(lm(sar.data$y ~x.mat[,-1]))$r.squared

summary(sar.model)
inv.rhoW = invIrM(neighbours = nb_DiseaseData_sf,rho = as.numeric(sar.model$rho))
incidence_hat = predict(sar.model,listw = lisw_DiseaseData_sf, pred.type = 'TC')
R.squared = 1-(var(CollectData2$incidence - incidence_hat))/var(CollectData2$incidence)
R.squared


library(dplyr)

library(stringr)
names.data = names(sar.data)

names.data = str_replace_all(names.data,'[()]','.')
names(sar.data) = names.data

names.data = names.data[-1]
VarG.list = NULL


BitMatrix <- function(n){
  # 作用：返还Bit矩阵
  # Args:n：数据长度
  set <- 0:(2^n-1)
  rst <- matrix(0,ncol = n,nrow = 2^n)
  for (i in 1:n){
    rst[, i] = ifelse((set-rowSums(rst*rep(c(2^((n-1):0)), each=2^n)))/(2^(n-i))>=1, 1, 0)
  }
  rst
}

k <- c('x1','x2','x3')
var_name = names(CollectData2)[3:5]
selectMatrix <- BitMatrix(length(k))
Varg.name = apply(selectMatrix, 1, function(x){ var_name[which(x==1)] })
Varg.name = str_replace_all(Varg.name,'"','')

sapply(1:length(Varg.name), function(i){paste0(Varg.name[[i]])})

#Construct VarG based on selectMatrix
for (ri in 2:nrow(selectMatrix)){
  VarG.temp = names.data
  for (ci in 1:length(selectMatrix[ri,])){
    var_name = k[ci]
    if (selectMatrix[ri,ci] == 1){
      VarG.temp = VarG.temp[str_detect(string = VarG.temp ,pattern = var_name)]
    }else{
      VarG.temp = VarG.temp[str_detect(string = VarG.temp ,pattern = var_name,negate = T)]
    }
  }
  VarG.list = append(VarG.list,list(VarG.temp))
}





VarGroups = VarG.list
data = sar.data

ols.shapley = shapley_value.ols(VarGroups = VarG.list, data = sar.data)
ols.shapley
Varg.name[-1]
interaction_detector('incidence',c('elevation','soiltype'),CollectData)
rbind(Varg.name[-1],ols.shapley)


sar.shapley = shapley_value.sar(VarGroups = VarG.list, data = sar.data,listW  = lisw_DiseaseData_sf)
Varg.name[-1]
sar.shapley
rbind(Varg.name[-1],sar.shapley)

