#install.packages("geodetector")
library(geodetector)
library(dplyr)
library(ggplot2)
library(tmap)
library(sf)
library(RColorBrewer)
library(maptools)
library(spatialreg)

library(RandomFields)
library('lattice')

simu2pm = function(simu,n=1){
  pm_grid = simu
  pm_grid$act = pm_grid@data[,n]
  pm_grid$act[is.na(pm_grid$act)] = 0
  
  pm_grid = sp::rbind.SpatialPixelsDataFrame(pm_grid)
  
  pm_point <- SpatialPointsDataFrame(coords = pm_grid@coords,
                                     data = dplyr::select(pm_grid@data, act),
                                     proj4string = pm_grid@proj4string)
  proj4string(pm_point) =
    CRS(" +proj=longlat +datum=WGS84 +no_defs")
  proj4string(pm_grid) =
    CRS(" +proj=longlat +datum=WGS84 +no_defs")
  pm_point$lon = pm_point@coords[,1]
  pm_point$lat = pm_point@coords[,2]
  return(c(pm_grid,pm_point))
}
rmse <- function(y, f) {
  sqrt(mean((y - f)^2))
}

construct_W <- function(data,lon_name,lat_name){
  
  #install.packages("geosphere")
  library(geosphere)
  W = distm(data[,c(lon,lat)],data[,c(lon,lat)],fun=distVincentyEllipsoid)

  d = distm(data[1,c(lon,lat)],data[5,c(lon,lat)],fun=distVincentyEllipsoid)
  W = exp(-(W^2)/(d^2))
  W[W<exp(-1)] = 0
  as.matrix(W)
  rw = rowSums(W)
  rw[rw==0] = 0.0000000001
  rw = 1/rw
  W*rw
  write.csv(as.matrix(W),'D:/0 job/0 q-vpc/data/W_dist.csv')
}
read_W <- function(path_W){
  W = read.table(path_W,header = T,sep = ',')
  W = W[,-1]
  W = as.matrix(W)

  #W <- as(W,"dgCMatrix")
}
#mat.list = construct_W(data = pm_point@coords,lon_name = 'coords.x1',lat__name = 'coords.x2')
path_W = 'D:/0 job/0 q-vpc/data/W_dist.csv'
mat.list = read_W(path_W)
mat.listW=mat2listw(mat.list)
rpi = 0
RFoptions(seed = NULL, height = 4)
breakP = 2
simuList_q = matrix(0,10,100)
simuList_r2 = matrix(0,10,100)
simuList_r2lm = matrix(0,10,100)
simuList_realr2 = matrix(0,10,100)
simuList_morans = matrix(0,10,100)

rp = 1
######## y_hat = γ_1 I+gama2*X2+gama3*x3
for (lamda in seq(0.1,0.9,0.1)) {
  for (simi in 1:10) {
    lon <- lat <- seq(1, 3, 0.05)
    
    model <- RMstable(alpha = rp, scale = 1)
    #model = RPgauss(model)
    
    simu_x <- RFsimulate(model, lon, lat)
    x_list = simu2pm(simu_x,1)
    x_point_1 = x_list[[2]]
    IW_solve = solve(diag(1,sim_size)-lamda*mat.list)
    #levelplot(x_point_1$act~x_point_1$lon+x_point_1$lat,
     #         xlab = 'lon',ylab = 'lat',pretty = T)
    x = x_point_1$act
    sim_size = nrow(x_point_1)

    mu = rnorm(sim_size,0,sqrt(0.5))

    y = 1 +  0.7 * x + IW_solve%*%mu
    
    
    
    #总方差
    sst <-  var(0.7 * x) + var(IW_solve%*%mu)
    
    #真实解释度
    real_R2 = var( 0.7 * x) / sst
    simuList_realr2[as.integer(lamda*10),simi] = real_R2
    #hist(y)
    
    
    sim_df = data.frame(y = y)
    
    library(classInt)
    x_break = classIntervals(var = x, n = 5, style = 'jenks', warnLargeN = F)
    x_break$var[x_break$var>=(x_break$brks[1]-1) & x_break$var<x_break$brks[2]] = 100
    x_break$var[x_break$var>=x_break$brks[2] & x_break$var<x_break$brks[3]] = 101
    x_break$var[x_break$var>=x_break$brks[3] & x_break$var<x_break$brks[4]] = 102
    x_break$var[x_break$var>=x_break$brks[4] & x_break$var<x_break$brks[5]] = 103
    x_break$var[x_break$var>=x_break$brks[5] & x_break$var<=(x_break$brks[6]+1)] = 104

    sim_df$x =  x_break$var
    q = factor_detector('y','x',sim_df)
    simuList_q[as.integer(lamda*10),simi] = unlist(q)[1]
    
   

    

    x.mat = model.matrix(lm(y ~ factor(x),data = sim_df))
    semData = data.frame(x0 = x.mat[,1],x1 = x.mat[,2],x2 = x.mat[,3],x3 = x.mat[,4],x4 = x.mat[,5])
    semData$y = sim_df$y
    names(semData)
    
    
    
    
    if (F) {
      se_model <- errorsarlm(y ~x1+x2+x3+x4, data = semData,
                             listw = mat.listW, zero.policy = T)
      y_hat =x.mat %*% se_model$coefficients
    }else{
      se_model <- errorsarlm(y ~ x,
                             listw =  mat.listW, zero.policy = T)
      y_hat = se_model$coefficients[1] + se_model$coefficients[2]*x
    }

    R_squared = 1 - sum((y-y_hat)^2)/sum((y-mean(y))^2)
    simuList_r2[as.integer(lamda*10),simi] = R_squared
    
    
  }
  print(lamda)
}

plot(rowMeans(simuList_q[1:9,]),rowMeans(simuList_r2[1:9,]))
plot(seq(0.1,0.9,0.1),rowMeans(simuList_q[1:9,])-rowMeans(simuList_r2[1:9,]),xlab = 'lambda',ylab = 'q - R^2',main = 'Difference between R^2 and q')

plot(seq(0.1,0.9,0.1),rowMeans(simuList_realr2[1:9,])-rowMeans(simuList_r2[1:9,]),xlab = 'lambda',ylab = 'Real_R^2 - R^2_estimated',main = 'Difference between Real_R^2 and R^2_estimated')
plot(seq(0.1,0.9,0.1),rowMeans(simuList_realr2[1:9,])-rowMeans(simuList_q[1:9,]),xlab = 'lambda',ylab = 'Real_R^2 - q',main = 'Difference between Real R^2 and q')

mean(abs(rowMeans(simuList_realr2[1:9,])-rowMeans(simuList_r2[1:9,])))
mean(abs(rowMeans(simuList_realr2[1:9,])-rowMeans(simuList_q[1:9,])))


mean((rowMeans(simuList_realr2[1:9,])-rowMeans(simuList_r2[1:9,])))
mean((rowMeans(simuList_realr2[1:9,])-rowMeans(simuList_q[1:9,])))

plot(seq(0.1,0.9,0.1),rowMeans(simuList_morans[1:9,]),xlab = 'lamda',ylab = "morans'I", main = "Relation between lamda and Morans'I")

simuList_q[1:9,1]
simuList_r2[1:9,1]
simuList_r2lm[1:9,1]



