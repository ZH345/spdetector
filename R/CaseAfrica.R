library(rgdal)
library(raster)
library(sf)
library(sp)
library(spdep)
library(dplyr)
library(terra)
library(classInt)
setwd('D:/BaiduSyncdisk/0 job/0 q-vpc/')
Africa_GreenWall = readOGR('./data/Africa_GreenWall.shp')
Africa_GreenWall = as(Africa_GreenWall,'sf')
ddi = readGDAL('./data/case_study/AFRICA_ddi2020_GW_Rs01.tif')


if (F){
  name.list = list.files('./data/case_study/')
  for (fn in name.list){
    data = readGDAL(paste0('./data/case_study/',fn))
    data = raster(data) %>%
      crop(.,y = Africa_GreenWall) %>%
      mask(.,mask =  Africa_GreenWall)
    raster::writeRaster(data,filename =  paste0('./dataout/',fn))
  }
}






raster2Category <- function(d.path,fn){
  data = readGDAL(paste0(d.path,fn))
  data = terra::rast(data)
  data = terra::resample(data,ddi.terra)
  #plot(data)

  # classIfy method for SpatRaster
  rcl = values(data$band1)[!is.na(values(data$band1))] %>% as.numeric() %>%
    classIntervals(var = ., n = 6, style = 'jenks')

  rcl = as.matrix(cbind(rcl$brks[1:(length(rcl$brks)-1)],
                        rcl$brks[2:(length(rcl$brks))],
                        1:(length(rcl$brks)-1)))
  rcl[1,1] = rcl[1,1]-1
  rcl[6,2] = rcl[6,2]+1
  data = classify(data, rcl, include.lowest=FALSE, right=TRUE,
                  others=NULL, brackets=TRUE)


  #data = raster(data)
  #plot(data)

  dna = strsplit(fn,'.t')[[1]][1]
  names(data) = dna
  data
  #rasterToPolygons(data)

}



ddi.terra = rast(ddi)
names(ddi.terra) = 'ddi'
if (F){
  d.path = './dataout/depenV/'
  files = list.files(d.path)
  for (fn in files){
    if (endsWith(fn,'.tif')){
      res = raster2Category(d.path,fn)
      ddi.terra = append(ddi.terra, res)
    }
  }
  names(ddi.terra) = c('ddi','dem','lc','pet','pop','pr','sm')
  terra::writeRaster(x = ddi.terra,filename = './dataout/africa_ddi.terra.tif')
}else{
  ddi.terra = readGDAL('./dataout/africa_ddi.terra.tif')
  ddi.terra = rast(ddi.terra)
}

names(ddi.terra) = c('ddi','dem','lc','pet','pop','pr','sm')
var.df = as.data.frame(values(ddi.terra))
names(var.df)
var.df = var.df[!is.na(var.df$ddi),]
nrow(var.df)

#ddi = crop(x = ddi,y = Africa_GreenWall)
#ddi = mask(x = ddi,mask =  Africa_GreenWall)
plot(ddi)
ddi = raster(ddi.terra$ddi)
ddi.poly = rasterToPolygons(ddi)
#writeOGR(obj = ddi.poly,dsn = './dataout/',layer = 'africa_ddipoly',driver = 'ESRI Shapefile')

ddi.poly = as(ddi.poly,'sf')
nrow(ddi.poly)


nrow(ddi.poly)
ddi.poly = ddi.poly[!is.na(ddi.poly$ddi),]
data.logit = is.na(var.df)[,1] | is.na(var.df)[,2] | is.na(var.df)[,3] | is.na(var.df)[,4]| is.na(var.df)[,5] | is.na(var.df)[,6]|is.na(var.df)[,7]
data.logit = !data.logit
var.df = var.df[data.logit,]
ddi.poly = ddi.poly[data.logit,]
ddi.nb = spdep::poly2nb(ddi.poly)
listw_ddi <- spdep::nb2listw(ddi.nb, style = "W", zero.policy = TRUE)
#print.listw(listw_ddi,zero.policy = T)
nrow(ddi.poly)
nrow(var.df)


x = c("dem", "lc", "pet", "pop", "pr" , "sm")

y = "ddi"
x = c("dem", "lc", "sm",'pop')
data = var.df
listW = listw_ddi
type = 'sar'
inter.term = T

shapley_value(y = "ddi",x = c("dem", "lc", "sm",'pop')
              ,data = var.df, listW = listw_ddi, type = 'sar'
              ,inter.term = T
              )





shapley_value <- function(
    y, ## y
    x, ## x name list
    data,       ## data to fit models
    type = 'sar',
    listW = NULL,
    inter.term = T
) {
  require(spatialreg)
  require(gtools)
  require(dplyr)
  require(stringr)

  #interaction_detector('incidence',c('elevation','soiltype','watershed'),CollectData)
  model.data = data[,c(y,x)]
  x.mat = NULL
  data.names = names(model.data)
  VarN = ncol(model.data)-1
  if (inter.term == T){
    ff = paste0(data.names[1],'~')
    for (xi in 1:(VarN)){
      if (xi < (VarN)){
        ff = paste0(ff,'factor(',data.names[xi+1],')','*')
        f =  paste0( paste0(data.names[1],'~'),'factor(',data.names[xi+1],')')
        x.mat = cbind(x.mat, model.matrix(lm(f,data = model.data))[,-1])
      }else{
        ff =paste0(ff,'factor(',data.names[xi+1],')')
        f =  paste0( paste0(data.names[1],'~'),'factor(',data.names[xi+1],')')
        x.mat = cbind(x.mat, model.matrix(lm(f,data = model.data))[,-1])
      }
    }
  }else{
    ff = paste0(data.names[1],'~')
    for (xi in 1:(VarN)){
      if (xi < (VarN)){
        ff = paste0(ff,'factor(',data.names[xi+1],')','+')
        f =  paste0( paste0(data.names[1],'~'),'factor(',data.names[xi+1],')')
        x.mat = cbind(x.mat, model.matrix(lm(f,data = model.data))[,-1])

      }else{
        ff =paste0(ff,'factor(',data.names[xi+1],')')
        f =  paste0( paste0(data.names[1],'~'),'factor(',data.names[xi+1],')')
        x.mat = cbind(x.mat, model.matrix(lm(f,data = model.data))[,-1])
      }
    }

  }

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

  getSubset2 = function(Var.Name){
    selectMatrix = BitMatrix(length(Var.Name))
    selectMatrix = selectMatrix[rowSums(selectMatrix)==2,]
    apply(selectMatrix,1, function(x){Var.Name[which(x==1)]})
  }


  if (inter.term == T){
    x.mat = NULL
    var.inter.2 = getSubset2(x)
    for (ci in 1:ncol(var.inter.2)){
      f =  paste0(paste0(data.names[1],'~'),
                  paste0('factor(',var.inter.2[1,ci],")*factor(",var.inter.2[2,ci]),")")
      x.mat = cbind(x.mat, model.matrix(lm(f,data = model.data))[,-1])

    }
  }else{
    for (xi in 1:(VarN)){
      if (xi < (VarN)){
        f =  paste0( paste0(data.names[1],'~'),'factor(',data.names[xi+1],')')
        x.mat = cbind(x.mat, model.matrix(lm(f,data = model.data))[,-1])
      }else{
        f =  paste0( paste0(data.names[1],'~'),'factor(',data.names[xi+1],')')
        x.mat = cbind(x.mat, model.matrix(lm(f,data = model.data))[,-1])
      }

    }
  }


  #x.mat = model.matrix(lm(ff,data = model.data))

  sar.data = data.frame(y = model.data[,1],  x.mat[,-1])



  #summary(sar.model)
  #inv.rhoW = invIrM(neighbours = nb_DiseaseData_sf,rho = as.numeric(sar.model$rho))
  #incidence_hat = predict(sar.model,listw = lisw_DiseaseData_sf, pred.type = 'TC')
  #R.squared = 1-(var(CollectData2$incidence - incidence_hat))/var(CollectData2$incidence)
  #R.squared

  if (inter.term == T){

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

    getSubset = function(Var.Name){
      selectMatrix = BitMatrix(length(Var.Name))
      Varg.name = apply(selectMatrix, 1, function(x){  Var.Name[which(x==1)] })
      Varg.name[-1]
    }
    names.sardata = names(sar.data)

    names.sardata = str_replace_all(names.sardata,'[()]','.')
    names(sar.data) = names.sardata
    Var.Name <- x
    Varg.name = getSubset(Var.Name = Var.Name)
    selectMatrix <- BitMatrix(length(Var.Name))
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
      print(model)
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
      print(model)
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














