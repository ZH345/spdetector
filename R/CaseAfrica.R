#git config --global http.sslVerify "false"

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
ddi = readGDAL('./dataout/AFRICA_ddi2020_GW_Rs05.tif')

# mask V.
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


raster2Category <- function(d.path,fn,togrid,typen){
  data = readGDAL(paste0(d.path,fn))
  data = terra::rast(data)
  data = terra::resample(data,togrid)

  #plot(data)

  # classIfy method for SpatRaster
  rcl = values(data$band1)[!is.na(values(data$band1))] %>% as.numeric() %>%
    classIntervals(var = ., n = typen, style = 'jenks')

  rcl = as.matrix(cbind(rcl$brks[1:(length(rcl$brks)-1)],
                        rcl$brks[2:(length(rcl$brks))],
                        1:(length(rcl$brks)-1)))

  rcl[1,1] = rcl[1,1]-1
  rcl[typen,2] = rcl[typen,2]+1
  print(rcl)
  data = classify(data, rcl, include.lowest=FALSE, right=TRUE,
                  others=NULL, brackets=TRUE)

  #plot(data)
  dna = strsplit(fn,'.t')[[1]][1]
  names(data) = dna
  data
}


#resample to D.V and lassify
ddi.terra = rast(ddi)
names(ddi.terra) = 'ddi'
if (F){
  d.path = './dataout/depenV/'
  files = list.files(d.path)
  for (fn in files){
    if (endsWith(fn,'.tif') & fn != 'AFRICA_LC.tif'){
      res = raster2Category(d.path,fn,togrid = ddi.terra, typen = 4)
      ddi.terra = append(ddi.terra, res)
    }else if (fn == 'AFRICA_LC.tif'){
      lc = readGDAL(paste0(d.path,fn))
      lc = terra::rast(lc)
      lc = terra::resample(lc,ddi.terra,method = 'near')
      ddi.terra = append(ddi.terra, lc)
    }
  }
  names(ddi.terra) = c('ddi','dem','lc','pet','pop','pr','sm')
  terra::writeRaster(x = ddi.terra,filename = './dataout/africa_ddi.terra.tif',overwrite=TRUE)
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
x = c("pr", "sm", "pop")
x = c("pr")
data = var.df
listW = listw_ddi
type = 'Durbin'
inter.term = 3

spdep::moran.test(var.df$ddi,listw_ddi,zero.policy=TRUE)
moran.plot(var.df$ddi,listw_ddi,zero.policy=TRUE)
lm.t = lm(ddi~ factor(pr) * factor(sm) *factor(pop) ,data = var.df)
spdep::lm.LMtests(lm.t, listw_ddi, test="all",zero.policy=TRUE)




model <- lagsarlm(ddi~ factor(pr) , data= var.df, listw = listW, zero.policy = T, Durbin=TRUE)
fity = predict(model, pred.type = 'TC',listw = listW, re.form = NA)
1 - var(var.df$ddi-fity)/var(var.df$ddi)

model <- lagsarlm(ddi~ factor(sm) , data= var.df, listw = listW, zero.policy = T, Durbin=TRUE)
fity = predict(model, pred.type = 'TC',listw = listW, re.form = NA)
1 - var(var.df$ddi-fity)/var(var.df$ddi)


model <- lagsarlm(ddi~ factor(pop) , data= var.df, listw = listW, zero.policy = T, Durbin=TRUE)
fity = predict(model, pred.type = 'TC',listw = listW, re.form = NA)
1 - var(var.df$ddi-fity)/var(var.df$ddi)

model <- lagsarlm(ddi~ factor(sm)*factor(pop), data= var.df, listw = listW, zero.policy = T, Durbin=TRUE)
fity = predict(model, pred.type = 'TC',listw = listW, re.form = NA)
1 - var(var.df$ddi-fity)/var(var.df$ddi)

model <- lagsarlm(ddi~ factor(pr)*factor(pop), data= var.df, listw = listW, zero.policy = T, Durbin=TRUE)
fity = predict(model, pred.type = 'TC',listw = listW, re.form = NA)
1 - var(var.df$ddi-fity)/var(var.df$ddi)

model <- lagsarlm(ddi~ factor(pr)*factor(sm), data= var.df, listw = listW, zero.policy = T, Durbin=TRUE)
fity = predict(model, pred.type = 'TC',listw = listW, re.form = NA)
1 - var(var.df$ddi-fity)/var(var.df$ddi)


## plot data

if(F){
  plot.var = 'ddi'
  ddi.terra.plot = raster(ddi.terra[plot.var]) %>% as(., "SpatialPixelsDataFrame") %>% as.data.frame()
  names(ddi.terra.plot) = c('act', "x", "y" )
  g1 = ggplot()+
    geom_raster(data = ddi.terra.plot, aes(x = x, y = y, fill = act))+
    scale_fill_gradient2(low = "steelblue",mid = 'yellow',high = 'red',midpoint = 200)+
    annotation_north_arrow(
      mapping = NULL,
      data = NULL,
      height = unit(1.5, "cm"),
      width = unit(1.5, "cm"),
      pad_x = unit(18, "cm"),
      pad_y = unit(22, "cm"),
      rotation = NULL,
      style = north_arrow_orienteering
    )+
    guides(color = guide_legend(override.aes = list(size = 4)))+

    #geom_spatial_point(aes(x, y),data = boxedge ,size=0)+
    annotation_scale(width_hint = 0.3,
                     plot_unit = 'km',
                     style = "ticks",
                     line_width = unit(3, "cm"),
                     height = unit(0.4, "cm"),
                     pad_x = unit(0.1, "cm"),
                     pad_y = unit(0.3, "cm")) +
    coord_sf(crs = 4326)+
    theme_bw()+
    theme(legend.position = c(0.95,0.7))+
    theme(legend.background = element_blank())+
    theme(panel.grid = element_blank(),
          #axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.title = element_text(size = 10,colour = 'white'),
          legend.text = element_text(size = 10))
  #maptheme
  #g1
  ggsave(paste0('D:/BaiduSyncdisk/0 job/0 q-vpc/image/Africa_', plot.var,'.pdf'),device = 'pdf',plot = g1,width = 8,height = 3,dpi = 1200)

  plot.var = 'pop'
  ddi.terra.plot = raster(ddi.terra[plot.var]) %>% as(., "SpatialPixelsDataFrame") %>% as.data.frame()
  names(ddi.terra.plot) = c('act', "x", "y" )
  g2 = ggplot()+
    geom_raster(data = ddi.terra.plot, aes(x = x, y = y, fill = as.character(act)))+
    scale_fill_discrete(labels = c("1: < 114.6", "2: <447.7", "3: < 1650.8", "4: < 3707.4"))+
    annotation_north_arrow(
      mapping = NULL,
      data = NULL,
      height = unit(1.5, "cm"),
      width = unit(1.5, "cm"),
      pad_x = unit(18, "cm"),
      pad_y = unit(22, "cm"),
      rotation = NULL,
      style = north_arrow_orienteering
    )+
    guides(color = guide_legend(override.aes = list(size = 4)))+

    #geom_spatial_point(aes(x, y),data = boxedge ,size=0)+
    annotation_scale(width_hint = 0.3,
                     plot_unit = 'km',
                     style = "ticks",
                     line_width = unit(3, "cm"),
                     height = unit(0.4, "cm"),
                     pad_x = unit(0.1, "cm"),
                     pad_y = unit(0.3, "cm")) +
    coord_sf(crs = 4326)+
    theme_bw()+
    theme(legend.position = c(0.92,0.79))+
    theme(legend.background = element_blank())+
    theme(panel.grid = element_blank(),
          #axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.title = element_text(size = 10,colour = 'black'),
          legend.text = element_text(size = 10))
  #maptheme
  #g2
  ggsave(paste0('D:/BaiduSyncdisk/0 job/0 q-vpc/image/Africa_', plot.var,'.pdf'),device = 'pdf',plot = g2,width = 8,height = 3,dpi = 1200)

  plot.var = 'pr'
  ddi.terra.plot = raster(ddi.terra[plot.var]) %>% as(., "SpatialPixelsDataFrame") %>% as.data.frame()
  names(ddi.terra.plot) = c('act', "x", "y" )
  g2 = ggplot()+
    geom_raster(data = ddi.terra.plot, aes(x = x, y = y, fill = as.character(act)))+
    scale_fill_discrete(labels = c("1: < 8.3", "2: <23.4", "3: < 49.4", "4: < 147.5"))+
    annotation_north_arrow(
      mapping = NULL,
      data = NULL,
      height = unit(1.5, "cm"),
      width = unit(1.5, "cm"),
      pad_x = unit(18, "cm"),
      pad_y = unit(22, "cm"),
      rotation = NULL,
      style = north_arrow_orienteering
    )+
    guides(color = guide_legend(override.aes = list(size = 4)))+

    #geom_spatial_point(aes(x, y),data = boxedge ,size=0)+
    annotation_scale(width_hint = 0.3,
                     plot_unit = 'km',
                     style = "ticks",
                     line_width = unit(3, "cm"),
                     height = unit(0.4, "cm"),
                     pad_x = unit(0.1, "cm"),
                     pad_y = unit(0.3, "cm")) +
    coord_sf(crs = 4326)+
    theme_bw()+
    theme(legend.position = c(0.92,0.79))+
    theme(legend.background = element_blank())+
    theme(panel.grid = element_blank(),
          #axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.title = element_text(size = 10,colour = 'black'),
          legend.text = element_text(size = 10))
  #maptheme
  #g2
  ggsave(paste0('D:/BaiduSyncdisk/0 job/0 q-vpc/image/Africa_', plot.var,'.pdf'),device = 'pdf',plot = g2,width = 8,height = 3,dpi = 1200)

  plot.var = 'sm'
  ddi.terra.plot = raster(ddi.terra[plot.var]) %>% as(., "SpatialPixelsDataFrame") %>% as.data.frame()
  names(ddi.terra.plot) = c('act', "x", "y" )
  g2 = ggplot()+
    geom_raster(data = ddi.terra.plot, aes(x = x, y = y, fill = as.character(act)))+
    scale_fill_discrete(labels = c("1: < 60.6", "2: < 157.6", "3: < 269.3", "4: < 445.0"))+
    annotation_north_arrow(
      mapping = NULL,
      data = NULL,
      height = unit(1.5, "cm"),
      width = unit(1.5, "cm"),
      pad_x = unit(18, "cm"),
      pad_y = unit(22, "cm"),
      rotation = NULL,
      style = north_arrow_orienteering
    )+
    guides(color = guide_legend(override.aes = list(size = 4)))+

    #geom_spatial_point(aes(x, y),data = boxedge ,size=0)+
    annotation_scale(width_hint = 0.3,
                     plot_unit = 'km',
                     style = "ticks",
                     line_width = unit(3, "cm"),
                     height = unit(0.4, "cm"),
                     pad_x = unit(0.1, "cm"),
                     pad_y = unit(0.3, "cm")) +
    coord_sf(crs = 4326)+
    theme_bw()+
    theme(legend.position = c(0.92,0.79))+
    theme(legend.background = element_blank())+
    theme(panel.grid = element_blank(),
          #axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.title = element_text(size = 10,colour = 'black'),
          legend.text = element_text(size = 10))
  #maptheme
  #g2
  ggsave(paste0('D:/BaiduSyncdisk/0 job/0 q-vpc/image/Africa_', plot.var,'.pdf'),device = 'pdf',plot = g2,width = 8,height = 3,dpi = 1200)

}




summary(model)
library(geodetector)

geodetector::interaction_detector(y,x,tabledata =var.df)


shapley_value(y = y,x = x
              ,data = var.df, listW = listw_ddi, type = 'sar'
              ,inter.term = inter.term
              )




shapley_value <- function(
    y, ## y
    x, ## x name list
    data,       ## data to fit models
    type = 'sar',
    listW = NULL,
    inter.term = 2
) {
  require(spatialreg)
  require(gtools)
  require(dplyr)
  require(stringr)
  require(data.table)
  require(stringr)
  model.data = data[,c(y,x)]
  x.mat = NULL
  data.names = names(model.data)
  VarN = ncol(model.data)-1
  if (VarN == 1){
    inter.term = 0
  }

  x.mat = NULL
  for (xi in 1:(VarN)){
    if (xi < (VarN)){
      f =  paste0( paste0(data.names[1],'~'),'factor(',data.names[xi+1],')')
      x.mat = cbind(x.mat, model.matrix(lm(f,data = model.data))[,-1])
    }else{
      f =  paste0( paste0(data.names[1],'~'),'factor(',data.names[xi+1],')')
      x.mat = cbind(x.mat, model.matrix(lm(f,data = model.data))[,-1])
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

  getSubset = function(Var.Name){
    selectMatrix = BitMatrix(length(Var.Name))
    Varg.name = apply(selectMatrix, 1, function(x){
      Varg.t= Var.Name[which(x==1)][1]
      if (length(Var.Name[which(x==1)]) == 1){
        return(Varg.t)
      }
      for (Varg in 2:length(Var.Name[which(x==1)])){
        Varg.t = paste0(Varg.t,',',Var.Name[which(x==1)][Varg])

      }
      Varg.t
      })
   Varg.name[-1]
  }
  x.mat = as.data.frame(x.mat)

  if (inter.term > 1){
    names.x.mat = names(x.mat)
    selectMatrix <- BitMatrix(length(x))
    selectMatrix <- selectMatrix[rowSums(selectMatrix)>1,]
    selectMatrix <- selectMatrix[rowSums(selectMatrix)<=inter.term,]
    for (ri in 1:nrow(selectMatrix)){
      Var.temp = x[which(selectMatrix[ri,] == 1)]
      x.mat.sublist = c()
      for (vari in 1:length(Var.temp)){
        VarG.temp1 = names.x.mat[str_detect(string =  names.x.mat, pattern = Var.temp[vari])]

        x.mat.1 = x.mat[,VarG.temp1]

        x.mat.sublist[vari] = list(x.mat.1)
      }
      if(F){#test
        a = matrix(1:12,3,4)
        b = matrix(2:13,3,4)
        x.mat.t = apply(a,2,function(x){
          apply(b, 2, function(y){
            x*y
          })
        })
        matrix(data = x.mat.t,nrow = 3,ncol = 16)
        a
        b
      }
      x.mat.1 = x.mat.sublist[[1]]

      for (vari in 2:length(Var.temp)){
        x.mat.2 = x.mat.sublist[[vari]]
        ncol = ncol(x.mat.1)*ncol(x.mat.2)
        names.x.mat.1 = names(x.mat.1)
        names.x.mat.2 = names(x.mat.2)
        names.x.mat.1 = sapply(X = names.x.mat.1,FUN = function(x){
          sapply(X = names.x.mat.2, FUN = function(y){
            paste0(x,':',y)
          })
        })

        names.x.mat.1 = as.character(names.x.mat.1)
        x.mat.1 = apply(x.mat.1,2,function(x){
                    apply(x.mat.2, 2, function(y){
                      x*y
                    })
                })

        x.mat.1 = matrix(data = x.mat.1,nrow = nrow(x.mat),ncol = ncol)
        x.mat.1 = data.frame(x.mat.1)
        names(x.mat.1) = names.x.mat.1
        x.mat.1 = x.mat.1[,colSums(x.mat.1)>0]
        print(vari)
      }

      x.mat.1 = data.frame(x.mat.1)
      names(x.mat.1)
      x.mat.1 = x.mat.1[,colSums(x.mat.1)>0]
      x.mat = cbind(x.mat,x.mat.1)

    }
  }

  ncol(x.mat)
  names(x.mat)
  #x.mat = model.matrix(lm(ff,data = model.data))

  sar.data = data.frame(y = model.data[,1],  x.mat[,-1])



  #summary(sar.model)
  #inv.rhoW = invIrM(neighbours = nb_DiseaseData_sf,rho = as.numeric(sar.model$rho))
  #incidence_hat = predict(sar.model,listw = lisw_DiseaseData_sf, pred.type = 'TC')
  #R.squared = 1-(var(CollectData2$incidence - incidence_hat))/var(CollectData2$incidence)
  #R.squared

  if (inter.term >1){


    names.sardata = names(sar.data)

    names.sardata = str_replace_all(names.sardata,'[()]','.')
    names(sar.data) = names.sardata
    Var.Name <- x
    Varg.name = getSubset(Var.Name = Var.Name)
    selectMatrix <- BitMatrix(length(Var.Name))
    selectMatrix <- selectMatrix[rowSums(selectMatrix)<=inter.term,]

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
    modellist[[1]] = paste0("model <- lagsarlm(y ~ 1, data= sar.data, listw = listW, zero.policy = T)")
    #modellist[[1]] = "model <- lm(y ~ 1 ,data= data)"
    ##function to compute R2
    computeR2<-function(model,sar.data){
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
    computeR2<-function(model,sar.data){
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
  }else if (type == 'Durbin'){

    ##the model list
    modellist<-lapply(namesindex, function(x)
      paste0("model <- lagsarlm(","y"," ~ ",
             paste(x, collapse = " + "),",data= sar.data, listw = listW, zero.policy = T, Durbin=TRUE)" ##paste variables entered
      )
    )
    modellist[[1]] = paste0("model <- lagsarlm(y ~ 1, data= sar.data, listw = listW, zero.policy = T, Durbin=TRUE)")
    #modellist[[1]] = "model <- lm(y ~ 1 ,data= data)"
    ##function to compute R2
    computeR2<-function(model,sar.data){
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
  }
  x2l = length(modellist)

  ##compute 2^p R2

  library(foreach)
  library(doParallel)
  # 创建一个集群并注册
  cl <- makeCluster(10)
  registerDoParallel(cl)

  # 启动并行计算
  tt = NULL
  time1 <- Sys.time()
  R2 = NULL
  for (moi in seq(10,length(modellist)+10,10)){
    x2 <- foreach(i = seq((moi-9),moi,1), .combine = cbind,
            .packages = c('geodetector','spatialreg','spdep'),
            .errorhandling = "stop") %dopar% {
              if (i <= length(modellist)){
                computeR2(modellist[[i]],sar.data)
              }
            }
    R2 = append(R2, x2)
    print((moi)/x2l)
  }

  #R2<-sapply(modellist,computeR2)


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










