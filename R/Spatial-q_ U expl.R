#install.packages("geodetector")
library(geodetector)
library(dplyr)
library(ggplot2)
library(tmap)
library(sf)
library(RColorBrewer)
library(maptools)
library(spatialreg)



data(DiseaseData_shp)
data(SoilType_shp)
data(Watershed_shp)
data(Elevation_shp)

nrow(DiseaseData_shp)
DiseaseData_sf = as(DiseaseData_shp,'sf')
#writeOGR(obj = DiseaseData_shp,dsn = 'D:/0 job/0 q-vpc/data/',layer = 'DieDiseaseData',driver = 'ESRI Shapefile')
tm4<-tm_shape(DiseaseData_sf)+tm_polygons(col = 'incidence',breaks=c(5.66,6.13,6.38,6.58),palette="RdYlGn")
tm4
DiseaseData_sf = DiseaseData_sf[,-2]
library(Matrix)
library(spdep)
library(classInt)
nb.list <- poly2nb(DiseaseData_sf)
mat.list <- nb2mat(nb.list,style="W")
listW = nb2listw(nb.list ,zero.policy = T)
W <- as(mat.list,"dgCMatrix")
plot(st_geometry(DiseaseData_sf),border = 'grey')
plot( st_centroid(DiseaseData_sf$geometry), add=TRUE )
plot(nb.list, st_centroid(DiseaseData_sf$geometry), add= TRUE )


######## y_hat
simuList_q = matrix(0,10,100)
simuList_r2 = matrix(0,10,100)
simuList_r2lm = matrix(0,10,100)
simuList_realr2 = matrix(0,10,100)
simuList_morans = matrix(0,10,100)

expl = 60
sim_size = nrow(DiseaseData_sf)


for (lamda in seq(0.1,0.9,0.1)) {
  for (simi in 1:100) {
    sim_size = nrow(DiseaseData_sf)
    x <- sample(1:3, size = sim_size, replace = TRUE)
    #解释度对应的方差
    
    x.mat = model.matrix(lm(mu ~ factor(x)))
    sigma_e = var(0.5 * x.mat[,2]+ 0.7 * x.mat[,3])/(expl/100)-var(0.5 * x.mat[,2]+ 0.7 * x.mat[,3])
    
    mu = rnorm(sim_size,0,sqrt(0.1/(lamda*10)))
    
    IW_solve = solve(diag(1,sim_size)-lamda*mat.list)
    
    
    y <- 5 + 0.5 * x.mat[,2]+ 0.7 * x.mat[,3] + IW_solve%*%mu
    
    y_hat_real =5 + 0.5 * x.mat[,2]+ 0.7 * x.mat[,3]
    
    real_R2 = 1 - sum((y-y_hat_real)^2)/sum((y-mean(y))^2)
    simuList_realr2[as.integer(lamda*10),simi] = real_R2
    
    sim_df = data.frame(y = y)
    sim_df$x = x
    q = factor_detector('y','x',sim_df)
    simuList_q[as.integer(lamda*10),simi] = unlist(q)[1]
    
    listW = nb2listw(nb.list ,zero.policy = T)
    
    
    semData = data.frame(x0 = x.mat[,1],x1 = x.mat[,2],x2 = x.mat[,3])
    semData$y = y
    names(semData)
    if (T){
      se_model <- errorsarlm(y ~ x, data = semData,
                             listw = listW, zero.policy = T)
      y_hat = se_model$coefficients[1] + se_model$coefficients[2]*x
      R_squared_con = 1 - sum((y-y_hat)^2)/sum((y-mean(y))^2)
      simuList_r2_continu[as.integer(lamda*10),simi] = R_squared_con
    }
    
    se_model <- errorsarlm(y ~x1+x2, data = semData,
                           listw = listW, zero.policy = T)
    
    y_hat =x.mat %*% se_model$coefficients
    
    R_squared = 1 - sum((y-y_hat)^2)/sum((y-mean(y))^2)
    simuList_r2[as.integer(lamda*10),simi] = R_squared
  }
  print(lamda)
}

plot(rowMeans(simuList_q[1:9,]),rowMeans(simuList_r2[1:9,]))
plot(seq(0.1,0.9,0.1),rowMeans(simuList_realr2[1:9,]),xlab = 'lambda',ylab = 'Real_R^2',main = 'Real_R^2')
plot(seq(0.1,0.9,0.1),rowMeans(simuList_realr2[1:9,])-rowMeans(simuList_q[1:9,]),xlab = 'lambda',ylab = 'Real_R^2 - Q',main = 'Difference between Real_R^2 and Q')
plot(seq(0.1,0.9,0.1),rowMeans(simuList_realr2[1:9,])-rowMeans(simuList_r2[1:9,]),xlab = 'lambda',ylab = 'Real_R^2 - SEM_R^2',main = 'Difference between Real_R^2 and SEM_R^2')
plot(seq(0.1,0.9,0.1),rowMeans(simuList_r2[1:9,])-rowMeans(simuList_q[1:9,]),xlab = 'lambda',ylab = 'SEM_R^2 - Q',main = 'Difference between SEM_R^2 and Q')


realR2_array = array(unlist(simuList_realr2[1:9,]))
q_array = array(unlist(simuList_q[1:9,]))
realR2_array_remove = realR2_array[realR2_array>0]
q_array_remove = q_array[realR2_array>0]

plot(realR2_array_remove,q_array_remove,xlab = 'Real R^2',ylab = 'Q',xlim = c(0,1),ylim = c(0,1),main = 'Correlation between Real_R^2 and Q')
abline(h = 0, col = "gray")
abline(v = 0, col = "gray")
abline(a = 0 , b = 1, col = "red")
mean(abs(realR2_array_remove-q_array_remove))



realR2_array = array(unlist(simuList_realr2[1:9,]))
r2_array = array(unlist(simuList_r2[1:9,]))
realR2_array_remove = realR2_array[realR2_array>0]
r2_array_remove = r2_array[realR2_array>0]

plot(realR2_array_remove,r2_array_remove,xlab = 'Real R^2',ylab = 'SEM R^2',xlim = c(0,1),ylim = c(0,1),main = 'Correlation between Real_R^2 and SEM_R^2')
abline(h = 0, col = "gray")
abline(v = 0, col = "gray")
abline(a = 0 , b = 1, col = "red")
mean(abs(realR2_array_remove-r2_array_remove))

# Non Factor
plot(seq(0.1,0.9,0.1),rowMeans(simuList_realr2[1:9,])-rowMeans(simuList_r2_continu[1:9,]),xlab = 'lambda',ylab = 'Real_R^2 - SEM_R^2_NonFactor',main = 'Difference between Real_R^2 and SEM_R^2_NonFactor')
plot(seq(0.1,0.9,0.1),rowMeans(simuList_r2_continu[1:9,])-rowMeans(simuList_q[1:9,]),xlab = 'lambda',ylab = 'SEM_R^2_NonFactor - Q',main = 'Difference between SEM_R^2_NonFactor and Q')

realR2_array = array(unlist(simuList_realr2[1:9,]))
r2_array_continu = array(unlist(simuList_r2_continu[1:9,]))
realR2_array_remove = realR2_array[realR2_array>0]
r2_array_continu_remove = r2_array_continu[realR2_array>0]

plot(realR2_array_remove,r2_array_continu_remove,xlab = 'Real R^2',ylab = 'SEM_R^2_NonFactor',xlim = c(0,0.6),ylim = c(0,0.6),main = 'Correlation between Real_R^2 and SEM_R^2_NonFactor')
abline(h = 0, col = "gray")
abline(v = 0, col = "gray")
abline(a = 0 , b = 1, col = "red")
mean(abs(realR2_array_remove-r2_array_remove))




simuList_q[1:9,1]
simuList_r2[1:9,1]
simuList_r2lm[1:9,1]
























