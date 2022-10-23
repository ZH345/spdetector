

###Multi


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
names(CollectData2)


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


## interaction_detector
VarGroups
### customized Shapley's value function
require(gtools)
f = incidence ~   elevation + soiltype
nb_DiseaseData_sf <- spdep::poly2nb(DiseaseData_sf)
listW <- spdep::nb2listw(nb_DiseaseData_sf, style = "W", zero.policy = TRUE)
data = CollectData2



listW = listW
y = 'incidence'
x = c('elevation','soiltype')
inter.term = F

shapley_value <- function(
    y, ## y
    x, ## x name list
    data,       ## data to fit models
    type = 'sar',
    listW = NULL,
    inter.term = T
) {


  require(gtools)
  #interaction_detector('incidence',c('elevation','soiltype','watershed'),CollectData)
  model.data = data[,c(y,x)]

  data.names = names(model.data)
  VarN = ncol(model.data)-1
  if (inter.term == T){
    ff = paste0(data.names[1],'~')
    for (xi in 1:(VarN)){
      if (xi < (VarN)){
        ff = paste0(ff,'factor(',data.names[xi+1],')','*')
      }else{
        ff =paste0(ff,'factor(',data.names[xi+1],')')
      }
    }
  }else{
    ff = paste0(data.names[1],'~')
    for (xi in 1:(VarN)){
      if (xi < (VarN)){
        ff = paste0(ff,'factor(',data.names[xi+1],')','+')
      }else{
        ff =paste0(ff,'factor(',data.names[xi+1],')')
      }
    }
  }

  ff = as.formula(ff)

  x.mat = model.matrix(lm(ff,data = model.data))

  sar.data = data.frame(y = model.data[,1],  x.mat[,-1])


  #summary(sar.model)
  #inv.rhoW = invIrM(neighbours = nb_DiseaseData_sf,rho = as.numeric(sar.model$rho))
  #incidence_hat = predict(sar.model,listw = lisw_DiseaseData_sf, pred.type = 'TC')
  #R.squared = 1-(var(CollectData2$incidence - incidence_hat))/var(CollectData2$incidence)
  #R.squared

  if (inter.term == T){

    require(dplyr)
    require(stringr)
    require(spatialreg)
    names.sardata = names(sar.data)

    names.sardata = str_replace_all(names.sardata,'[()]','.')
    names(sar.data) = names.sardata

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

    Var.Name <- data.names[-1]
    selectMatrix <- BitMatrix(length(Var.Name))
    Varg.name = apply(selectMatrix, 1, function(x){  Var.Name[which(x==1)] })
    Varg.name = unlist(Varg.name)[-1]

    VarG.list = NULL
    #Construct VarG based on selectMatrix
    for (ri in 2:nrow(selectMatrix)){
      VarG.temp = names.sardata
      for (ci in 1:length(selectMatrix[ri,])){
        var_name =  Var.Name[ci]
        if (selectMatrix[ri,ci] == 1){
          VarG.temp = VarG.temp[str_detect(string = VarG.temp ,pattern = var_name)]
        }else{
          VarG.temp = VarG.temp[str_detect(string = VarG.temp ,pattern = var_name,negate = T)]
        }
      }
      VarG.list = append(VarG.list,list(VarG.temp))
    }
  }else{
    require(dplyr)
    require(stringr)
    require(spatialreg)
    names.sardata = names(sar.data)

    names.sardata = str_replace_all(names.sardata,'[()]','.')
    names(sar.data) = names.sardata
    BitMatrix <- function(n){
      # 作用：返还Bit矩阵
      # Args:n：数据长度
      diag(1,n)
    }

    Var.Name <- data.names[-1]
    selectMatrix <- BitMatrix(length(Var.Name))
    Varg.name = apply(selectMatrix, 1, function(x){  Var.Name[which(x==1)] })
    Varg.name = unlist(Varg.name)


    VarG.list = NULL
    #Construct VarG based on selectMatrix
    for (ri in 1:length(x)){
      VarG.temp = names.sardata
      for (ci in 1:length(selectMatrix[ri,])){
        var_name =  Var.Name[ci]
        if (selectMatrix[ri,ci] == 1){
          VarG.temp = VarG.temp[str_detect(string = VarG.temp ,pattern = var_name)]
        }else{
          VarG.temp = VarG.temp[str_detect(string = VarG.temp ,pattern = var_name,negate = T)]
        }
      }
      VarG.list = append(VarG.list,list(VarG.temp))
    }

  }



  ## extract variable groups and formula
  Varfull <- VarG.list
  J <- length(Varfull)

  ##generate possible combinations as input profiles
  index<-eval(parse(text = paste0("expand.grid(",paste(rep("c(0,1)",J), collapse=", "),")")))

  ##transform the combinations to the variable names included
  namesindex<-apply(index,1,function(x) unlist(Varfull[which(x==1)]))


  if (type == 'sar'){
    ##the model list
    modellist<-lapply(namesindex, function(x)
      paste0("model <- lagsarlm(","y"," ~ ",
             paste(x, collapse = " + "),",data= sar.data, listw = listW, zero.policy = T)" ##paste variables entered
      )
    )
    modellist[[1]] = paste0("model <- lagsarlm(", "y", " ~ 1, data= sar.data, listw = listW, zero.policy = T)")
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
      R2 <- 1 - var(model.data[,1] - fitted_y)/var(model.data[,1])

      return(R2)
    }
  }else if (type == 'ols'){
    ##the model list
    modellist<-lapply(namesindex, function(x)
      paste0("model <- lm(","y"," ~ ",
             paste(x, collapse = " + "),",data= sar.data)" ##paste variables entered
      )
    )

    modellist[[1]] = paste0("model <- lm(", "y", " ~ 1, data= sar.data)")
    ##function to compute R2
    computeR2<-function(model){
      modelfit<-eval(parse(text = model))
      ### Calculated the variance of y explained by variables
      #fitted_y <- predict(modelfit, re.form = NA)
      fitted_y <- predict(modelfit)
      #var_y_f <- var(fitted_y)
      ### Extract variances of random effects
      #temp.var <- as.data.frame(VarCorr(modelfit))$vcov
      # Note the data has a four-level structure so there would be 4 variances
      R2 <- 1 - var(model.data[,1] - fitted_y)/var(model.data[,1])

      return(R2)
    }

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
  S.percentage = S/sum(S)*100
  return(t(rbind(Varg.name, S,S.percentage)))

}


pwy = rnorm(1000,0,sqrt(3))
e = rnorm(1000,0,sqrt(5))

cor.c = as.numeric(cor.test(pwy,e)$estimate)

var(pwy+e) - (var(e)+var(pwy)+2*cor.c*sqrt(var(e))*sqrt(var(pwy)))
var(pwy*e)
