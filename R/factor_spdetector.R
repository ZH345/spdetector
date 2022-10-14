library(spatialreg)
library(geodetector)


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



## interaction_detector

### customized Shapley's value function
require(gtools)
f2 = incidence ~   elevation * soiltype
f1 = incidence ~ elevation

nb_DiseaseData_sf <- spdep::poly2nb(DiseaseData_sf)
listW <- spdep::nb2listw(nb_DiseaseData_sf, style = "W", zero.policy = TRUE)
data = CollectData2





f= f1
listw = listW
factor_spdector < - function(f, data,type='sar',listw = NULL){
  model.data = model.frame(f,data)
  data.names = names(model.data)

  VarN = ncol(model.data)-1

  ff = paste0(data.names[1],'~')
  for (xi in 1:(VarN)){
    if (xi < (VarN)){
      ff = paste0(ff,'factor(',data.names[xi+1],')','*')
    }else{
      ff =paste0(ff,'factor(',data.names[xi+1],')')
    }

  }

  ff = as.formula(ff)
  if (type == 'sar'){
    sar.fit = lagsarlm(ff,data= data, listw = listW, zero.policy = T)
    y_hat = predict(sar.fit,listw = listw, pred.type = 'TC')
    R.squared = 1-(var(model.data[,1] - y_hat))/var(model.data[,1])
    return(R.squared)

  }else{
    if (type == 'ols'){
      lm.fit = lm(ff,data= data)
      y_hat = predict(lm.fit)
      R.squared = 1-(var(model.data[,1] - y_hat))/var(model.data[,1])
      return(R.squared)
    }
  }

}













