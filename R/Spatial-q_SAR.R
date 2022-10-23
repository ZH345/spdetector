#install.packages("geodetector")

  library(dplyr)
  library(ggplot2)
  library(tmap)
  library(sf)
  library(RColorBrewer)
  library(maptools)


library(geodetector)
library(spatialreg)
library(spdep)

setwd('D:/0 job/0 research/spdetector/')
data(DiseaseData_shp)
data(SoilType_shp)
data(Watershed_shp)
data(Elevation_shp)

nrow(DiseaseData_shp)
DiseaseData_sf = as(DiseaseData_shp,'sf')
#writeOGR(obj = DiseaseData_shp,dsn = 'D:/0 job/0 q-vpc/data/',layer = 'DieDiseaseData',driver = 'ESRI Shapefile')
#tm4<-tm_shape(DiseaseData_sf)+tm_polygons(col = 'incidence',breaks=c(5.66,6.13,6.38,6.58),palette="RdYlGn")
#tm4
DiseaseData_sf = DiseaseData_sf[,-2]

nb.list <- poly2nb(DiseaseData_sf)
mat.list <- nb2mat(nb.list,style="W")
#test for different W
if (F){
  dim(mat.list)
  mat.list = matrix(data = sample(x =0:1,size =  (dim(mat.list)[1])^2,replace = T,prob = c(0.99,0.01)),nrow = dim(mat.list)[1],ncol = dim(mat.list)[1])
  rsm = rowSums(mat.list)
  mat.list[,7][rsm==0] = 1
  rsm = rowSums(mat.list)
  sum(rsm)
  mat.list = mat.list/rsm
  rowSums(mat.list)
}

listW = mat2listw(mat.list)
plot(st_geometry(DiseaseData_sf),border = 'grey')
plot( st_centroid(DiseaseData_sf$geometry), add=TRUE )
plot(listW, st_centroid(DiseaseData_sf$geometry), add= TRUE )



data(SoilType_shp)
data(Watershed_shp)
data(Elevation_shp)
SoilType_sf = as(SoilType_shp,'sf')
Watershed_sf = as(Watershed_shp,'sf')
Elevation_sf = as(Elevation_shp,'sf')

plot(st_geometry(SoilType_sf$geometry),border = 'black')
plot(st_geometry(Watershed_sf$geometry), add=TRUE ,border = 'black')
plot(st_geometry(Elevation_sf$geometry), add= TRUE ,border = 'black')


######## y_hat tttttttttttttttttttttttttttttttttttttttttt
simuList_q = matrix(0,100,1000)
simuList_r2 = matrix(0,100,1000)
simuList_MI = matrix(0,100,1000)
simuList_rho = matrix(0,100,1000)


SAR.Sim <- function(rho){
  lapply(1:1000,function(simi){
    sim_size = dim(mat.list)[1]
    x <- sample(1:3, size = sim_size, replace = TRUE)
    #解释度对应的方差
    #sigma_e = 0.25 * var(x)/(expl/100)-0.25 * var(x)
    mu = rnorm(sim_size,0,sqrt(0.1))
    x.mat = model.matrix(lm(mu ~ factor(x)))

    IW_solve = solve(diag(1,sim_size)-rho*mat.list)

    y <- IW_solve%*%(5 + 0.5 * x.mat[,2]+ 0.7 * x.mat[,3] + mu)
    #y <- IW_solve%*%(5 + 0.3 * x.mat[,2]+ 0.5 * x.mat[,3] + mu)

    sim_df = data.frame(y = y)
    sim_df$x = x
    q = factor_detector('y','x',sim_df)
    simuList_q[as.integer(rho*100),simi] = unlist(q)[1]
    stop('STOP!')
    if(T){
      #moran.plot(as.vector(y), listW,
      #          zero.policy = T,
      #         labels = F,
      #        pch = 20, cex = 0.1)
      mi = moran.test(as.vector(y), listW)
      simuList_MI[as.integer(rho*100),simi] = as.numeric(unlist(mi$estimate[1]))
    }

    sarData = data.frame(x0 = x.mat[,1],x1 = x.mat[,2],x2 = x.mat[,3])
    sarData$y = y

    sar_model <- lagsarlm(y ~x1+x2, data = sarData,
                          listw = listW, zero.policy = T)

    simuList_rho[as.integer(rho*100),simi] = sar_model$rho
    y_hat = predict(sar_model,listw = listW, pred.type = 'TC')
    R_squared = 1 - var((y-y_hat))/var(y)
    simuList_r2[as.integer(rho*100),simi] = R_squared
    print(simi)
  })
  print(rho)
}




simuList_q = matrix(0,1,1000)
simuList_r2 = matrix(0,1,1000)
simuList_MI = matrix(0,1,1000)
simuList_rho = matrix(0,1,1000)

simuList_q_all = matrix(0,100,1000)
simuList_r2_all = matrix(0,100,1000)
simuList_MI_all = matrix(0,100,1000)
simuList_rho_all = matrix(0,100,1000)

rholist =  seq(0.01,0.99,0.01)
library(foreach)
library(doParallel)
# 创建一个集群并注册
cl <- makeCluster(15)
registerDoParallel(cl)

# 启动并行计算
tt = NULL
time1 <- Sys.time()
for (rho in seq(1,1,0.1)){
  simuList_q = matrix(0,1,1000)
  simuList_r2 = matrix(0,1,1000)
  simuList_MI = matrix(0,1,1000)
  simuList_rho = matrix(0,1,1000)
  x2 <- foreach(i = seq((rho-0.09),rho,0.01), .combine = rbind,.packages = c('geodetector','spatialreg','spdep'),.errorhandling = "stop") %dopar% {
    IW_solve = invIrM(neighbours = nb.list,style = 'W',rho = i,method = 'solve',)
    for (simi in 1:1000){
      if (i != 1){
        sim_size = dim(mat.list)[1]
        x <- sample(1:3, size = sim_size, replace = TRUE)
        #解释度对应的方差
        #sigma_e = 0.25 * var(x)/(expl/100)-0.25 * var(x)

        x.mat = model.matrix(lm(mu ~ factor(x)))
        mu = rnorm(sim_size,0,sqrt(0.1))
        #var(mu)
        IW_solve = invIrM(neighbours = nb.list,style = 'W',rho = i,method = 'solve',)
        y <- IW_solve%*%(5 + 0.5 * x.mat[,2]+ 0.7 * x.mat[,3] + mu)
        #y <- IW_solve%*%(5 + 0.3 * x.mat[,2]+ 0.5 * x.mat[,3] + mu)
        SST = var(IW_solve%*%(5 + 0.5 * x.mat[,2]+ 0.7 * x.mat[,3]))+var(IW_solve%*%mu)

        design.R2 = 1-var(IW_solve%*%mu)/SST


        sim_df = data.frame(y = y)
        sim_df$x = x
        q = factor_detector('y','x',sim_df)
        q = as.numeric(unlist(q)[1])
        simuList_q[1,simi] = q

        if(T){
          #moran.plot(as.vector(y), listW,
          #          zero.policy = T,
          #         labels = F,
          #        pch = 20, cex = 0.1)
          mi = moran.test(as.vector(y), listW)
          mi.num= as.numeric(unlist(mi$estimate[1]))
          simuList_MI[1,simi] = mi.num
        }

        sarData = data.frame(x0 = x.mat[,1],x1 = x.mat[,2],x2 = x.mat[,3])
        sarData$y = y

        sar_model <- lagsarlm(y ~x1+x2, data = sarData,
                              listw = listW, zero.policy = T)
        lm.model = lm(y ~x1+x2, data = sarData)

        summary(lm.model)
        simuList_rho[1,simi] = sar_model$rho
        rho.num = as.numeric(sar_model$rho)
        y_hat = predict(sar_model,listw = listW, pred.type = 'TC')
        R_squared = as.numeric(1 - var((y-y_hat))/var(y))
        simuList_r2[1,simi] = R_squared

        if (F){
          lm.model$coefficients
          sar_model$coefficients
          var(mu)
          rho.num
          q-summary(lm.model)$r.squared
          (q-design.R2)/design.R2
          (q-R_squared)/R_squared
        }

      }
    }
    return_data = cbind(simuList_q,simuList_r2,simuList_rho,simuList_MI)
    return(return_data)
  }

  simuList_q_all[((rho-0.09)*100):((rho)*100),seq(1,1000)] = x2[1:10,seq(1,1000)]
  simuList_r2_all[((rho-0.09)*100):((rho)*100),seq(1,1000)] = x2[1:10,seq(1001,2000)]
  simuList_rho_all[((rho-0.09)*100):((rho)*100),seq(1,1000)] = x2[1:10,seq(2001,3000)]
  simuList_MI_all[((rho-0.09)*100):((rho)*100),seq(1,1000)] = x2[1:10,seq(3001,4000)]
  print(rho)
}




time2 <- Sys.time()
print(time2-time1)
# Time difference of 19.15744 secs

# 在计算结束后别忘记关闭集群
stopImplicitCluster()
stopCluster(cl)

if(F){

  simuList_q_csv = read.table('./data/simulation/0 simuList_q_sar_1000_simg0.1.csv',header = T,sep = ',')
  simuList_q_csv = simuList_q_csv[,-1]
  simuList_r2_csv = read.table('./data/simulation/0 simuList_r2_sar_1000_simg0.1.csv',header = T,sep = ',')
  simuList_r2_csv = simuList_r2_csv[,-1]
  simuList_rho_csv = read.table('./data/simulation/0 simuList_rho_sar_1000_simg0.1.csv',header = T,sep = ',')
  simuList_rho_csv = simuList_rho_csv[,-1]
  simuList_MI_csv = read.table('./data/simulation/0 simuList_MI_sar_1000_simg0.1.csv',header = T,sep = ',')
  simuList_MI_csv = simuList_MI_csv[,-1]


  if (F){ #test map
    #MI & rho
    plot(rowMeans(simuList_rho_csv[2:98,]),rowMeans(simuList_MI_csv[2:98,]),
         xlab = 'lambda',ylab = 'SEM_R^2 - Q',
         main = 'Difference between SEM_R^2 and Q',
         #xlim = c(0,1),ylim = c(0,-100)
    )
    #real rho & diff
    plot(seq(0.01,0.99,0.01),(rowMeans(simuList_q_csv[2:98,])-rowMeans(simuList_r2_csv[2:98,]))/rowMeans(simuList_r2_csv[2:98,])*100,
         xlab = 'lambda',ylab = 'SEM_R^2 - Q',
         main = 'Difference between SEM_R^2 and Q',
         #xlim = c(0,1),ylim = c(0,-100)
    )
    #rho_hat & diff
    plot(rowMeans(simuList_rho_csv[2:98,]),(rowMeans(simuList_q_csv[2:98,])-rowMeans(simuList_r2_csv[2:98,]))/rowMeans(simuList_r2_csv[2:98,])*100,
         xlab = 'lambda',ylab = 'SEM_R^2 - Q',
         main = 'Difference between SEM_R^2 and Q',
         #xlim = c(0,1),ylim = c(0,-100)
    )

    #MI & diff
    plot(rowMeans(simuList_MI_csv[2:98,]),(rowMeans(simuList_q_csv[2:98,])-rowMeans(simuList_r2_csv[2:98,])),
         xlab = 'lambda',ylab = 'SEM_R^2 - Q',
         main = 'Difference between SEM_R^2 and Q',
         #xlim = c(0,1),ylim = c(0,-100)
    )

    trans.Lambda = -5*rowMeans(simuList_MI_csv[2:98,])^1.8
    trans.Lambda.norm = (trans.Lambda-min(trans.Lambda))/(max(trans.Lambda)-min(trans.Lambda))
    diff.r2q = (rowMeans(simuList_q_csv[2:98,])-rowMeans(simuList_r2_csv[2:98,]))
    diff.r2q.norm = (diff.r2q-min(diff.r2q))/(max(diff.r2q)-min(diff.r2q))
    plot(trans.Lambda,diff.r2q
         ,xlab = 'lambda',ylab = 'SEM_R^2 - Q',
         main = 'Difference between SEM_R^2 and Q',
         #xlim = c(0,1),ylim = c(0,-80)
    )
    abline(h = 0, col = "gray")
    abline(v = 0, col = "gray")
    abline(a = 0 , b =1 , col = "red")


    eigen.W.solve = eigen(IW_solve)
    var(eigen.W.solve$values)
    t(eigen.W.solve$values)%*%eigen.W.solve$values


    x = c()
    for (rho in seq(0.01,0.90,0.01)){
      IW_solve = solve(diag(1,sim_size)-rho*mat.list)
      eigen.W.solve = eigen(IW_solve)
      var(eigen.W.solve$values)
      wst = unlist(t(eigen.W.solve$values)%*%eigen.W.solve$values)
      x = append(x,wst)
    }


    y = rowMeans(simuList_q_csv[2:98,])-rowMeans(simuList_r2_csv[2:98,])
    x = rowMeans(simuList_MI_csv[2:98,])


    m <- nls(y ~ a*x^(b),
             start = list(a= -1.9,b=2))


    plot(x,y,
         xlab = 'rho',ylab = 'Q -R^2 (%)',
         main = 'Difference between Q and SAR_R^2',
         #xlim = c(0,1),ylim = c(0,-100)
    )
    #lines(x,y = trans.Lambda,col='green')
    lines(x,y = predict(m),col='red')
    cor.test(predict(m),y)$estimate^2

  }else{

    #MI & diff


    library(ggplot2)
    ggdata = data.frame(x=rowMeans(simuList_MI_csv[2:98,]))
    ggdata$y = (rowMeans(simuList_q_csv[2:98,])-rowMeans(simuList_r2_csv[2:98,]))
    m <- nls(y ~ a*x^(b),data = ggdata,
             start = list(a= -1.9,b=2))

    g1 = ggplot()+
      geom_point(aes(x = x,y = y),data = ggdata,size=0.5)+
      geom_smooth(aes(x = x,y = y), data = ggdata,color = 'red',method='lm',formula = y ~ I(m$m$getAllPars()[1]*(x^m$m$getAllPars()[2])),
                  se = T,level = 0.99,size = 0.3
                  )

    #ggplot2::ggsave('D:/0 job/0 q-vpc/results_simulation/MI_qr2_diff.pdf',device = 'pdf',plot = g1,width = 8,units = 'cm',height = 7.5,dpi = 1200)


    #rho & diff


    ggdata = data.frame(x=rowMeans(simuList_rho_csv[2:98,]))
    ggdata$y = (rowMeans(simuList_q_csv[2:98,])-rowMeans(simuList_r2_csv[2:98,]))/rowMeans(simuList_r2_csv[2:98,])
    m <- nls(y ~ a*x^(b),data = ggdata,
             start = list(a= -1.9,b=2))

    m$m$getAllPars()[1]
    m$m$getAllPars()[2]
    cor.test(predict(m),ggdata$y)$estimate^2
    g2 = ggplot()+
      geom_point(aes(x = x,y = y),data = ggdata,size=0.5)+
      geom_smooth(aes(x = x,y = y), data = ggdata,color = 'red',method='lm', formula= y~ I(m$m$getAllPars()[1]*(x^m$m$getAllPars()[2])),
                  se = T,level = 0.99,size = 0.3
      )
    ggplot2::ggsave('D:/0 job/0 q-vpc/results_simulation/rho_qr2_diff.pdf',device = 'pdf',plot = g2,width = 8,units = 'cm',height = 7.5,dpi = 1200)

    #rho & MI

    summary.diff <- apply(simuList_MI_csv,1,quantile,probs=c(0.5,0.025,0.975),na.rm = T)[,2:98]

    plot.data <- data.frame(rho = rho,diff50 = summary.diff[1,],
                            diff25 = summary.diff[2,],
                            diff975 = summary.diff[3,])
    rho =  seq(0.02,0.98,0.01)
    ## Fit a nonlinear least square model
    diff.model <- nls(diff50 ~ b * rho ^ c, data = plot.data,
                      start=list(b = -0.5, c = 2))

    summary(diff.model)
    xx <- predict(diff.model)
    cor(xx,plot.data$diff50)^2

    plot.data$diff50.fit <- xx

    g1 = ggplot(data = plot.data, aes(x = rho)) +
      geom_linerange(aes(ymin=diff25, ymax=diff975),show.legend = T)+
      geom_line(aes(y=diff50.fit,color="Fitted curve"),size=0.5) +
      geom_point(aes(y=diff50,shape="Diff with lower limit (2.5%) and upper limit (97.5%)"),color = 'black')+
      #geom_line(aes(y=diff25.fit,color="Lower limit (2.5%)"),linetype="dashed",size=0.1) +
      #geom_line(aes(y=diff975.fit,color="Upper limit (97.5%)"),linetype="dashed",size=0.1) +
      #ylim(-0.8,0.01)+
      #xlim(0.05,1)+
      scale_x_continuous(breaks = scales::breaks_width(0.2, offset = 0.05))+
      #geom_hline(yintercept = 0,color = 'red4',size = 0.3,lty = "dashed")+
      #geom_hline(yintercept = -1,color = 'red4',size = 0.3,lty = "dashed")+
      theme_bw()+
      theme(legend.position = c(0.4,0.3))+
      theme(legend.background = element_rect(fill = "transparent",colour = NA))

    g1



    ggplot2::ggsave('D:/0 job/0 q-vpc/results_simulation/rho_MI_diff.pdf',device = 'pdf',plot = g3,width = 8,units = 'cm',height = 7.5,dpi = 1200)


    ### quantile
    diff.quantile = data.frame(rho = (2:98)/100,diff.2.5 = rep(NA,97), diff.25 = rep(NA,97), diff.50 = rep(NA,97),diff.75 = rep(NA,97),diff.97.5 = rep(NA,97))
    diff.mat = ((simuList_q_csv)-(simuList_r2_csv))/(simuList_r2_csv)
    for (rho in 2:98){
      tq =  quantile(diff.mat[rho,],probs = c(0.025,0.20,0.50,0.80,0.975))

      diff.quantile[(rho-1),2:6] = tq

    }

    plot.data <- data.frame(rho = diff.quantile[,1],diff50 = diff.quantile[,4],
                            diff25 = diff.quantile[,3],
                            diff975 = diff.quantile[,5])

    ## Fit a nonlinear least square model
    diff.model <- nls(diff50 ~ b * rho^c, data = plot.data,
                      start=list(b = 0.5, c = -0.2))

    summary(diff.model)
    xx <- predict(diff.model)
    cor(xx,plot.data$diff50)^2

    plot.data$diff50.fit <- xx
    ggplot(data = plot.data, aes(x = rho)) +
      geom_linerange(aes(ymin=diff25, ymax=diff975))+
      geom_point(aes(y=diff50))+
      geom_line(aes(y=diff50.fit), col = "red") +
      theme_classic()



    plot(diff.quantile$rho,diff.quantile$diff.97.5)
    write.csv(diff.quantile,file = 'D:/0 job/0 q-vpc/data/simulation/q_r2.quantile.csv')

    #denoise
    for (rho in seq(0.01,0.99,0.01)){
      simuList_rho_csv[as.integer(rho*100),][simuList_rho_csv[as.integer(rho*100),]< (rho-0.01) | simuList_rho_csv[as.integer(rho*100),]> (rho+0.01)]
      diff.quantile[(rho-1),2:6] = tq
    }

    #vector
    diff.mat = ((simuList_q_csv)-(simuList_r2_csv))/(simuList_r2_csv)
    diff.mat = ((simuList_q_csv))/(simuList_r2_csv)
    diff.vector= unlist(diff.mat)
    diff.vector[diff.vector>2] = NA
    diff.vector[diff.vector< -2] = NA
    simuList_MI_vector = unlist(simuList_MI_csv)
    plot(simuList_MI_vector,diff.vector)
    data.vector = data.frame(x = simuList_MI_vector,y = diff.vector)
    data.vector = data.vector[!is.na(data.vector$y),]

    data.vector = data.vector[data.vector$x>0,]
    #data.vector = data.vector[data.vector$x>0,]

    fit <- nls(y ~ a*x^(b),data = data.vector,
             start = list(a= -1.9,b=2))
    cor.test(predict(fit),data.vector$y)$estimate^2

    library(propagate) #加载propagate包用于误差传播
    library(investr) #调用predFit函数
    pred_nls<-function(m){
      r=data.frame(x=m) # 此处注意数据框的列名"x"应与模型中的自变量名称"x"相同，否则报错
      result=predictNLS(fit,r,nsim=10000, interval="conf")
      #其中，nsim为蒙特卡洛模拟次数，该值默认为10000，某些情况下可设置为100000，该值越大，运算速度越慢；interval参数格式只接受字符串，"conf"表示置信区间，"pred"表示预测区间
      return(result$summary)
    }

    x_pred = seq(0,1,0.01)
    y_prednls<-pred_nls(x_pred) #将自变量取值代入函数，运行可看到误差传播过程
    y_predfit<-as.data.frame(predFit(fit,newdata=data.frame(x=x_pred),interval="confidence"))

    g2 = ggplot() +
      geom_point(aes(x=x,y=y),data.vector, size = 0.001)+
      #geom_density(aes(x=x,y=y))+
      #geom_line(size = 0.5) +
      #scale_fill_continuous(type = "viridis") +
      #geom_abline(intercept=0.75,slope=0,color = 'black')+
      geom_hline(yintercept = 0,color = 'red4',size = 0.3,lty = "dashed")+
      geom_hline(yintercept = -1,color = 'red4',size = 0.3,lty = "dashed")+
      geom_line(aes(x=x_pred,y=y_prednls$Sim.Median,color="Fitted curve"),size=0.1) +
      geom_line(aes(x=x_pred,y=y_prednls$'Sim.2.5%',color="lower limit (2.5%)"),linetype="dashed",size=0.1) +
      geom_line(aes(x=x_pred,y=y_prednls$'Sim.97.5%',color="lower limit (2.5%)"),linetype="dashed",size=0.1) +

      #geom_line(aes(x=x_pred,y=y_predfit$fit,color="Fitted curve"),size=0.1) +
      #geom_line(aes(x=x_pred,y=y_predfit$lwr,color=,linetype="dashed",size=0.1) +
      #geom_line(aes(x=x_pred,y=y_predfit$upr,color="lower limit (2.5%)"),linetype="dashed",size=0.1) +

      #geom_smooth(method = 'lm',na.rm = T,formula=y~I(x^2.0224 ),se = T,level = 0.95,color = 'red')+
      #geom_smooth(method = 'glm',na.rm = T,formula=y~x,se = T,color = 'red')+
      #scale_color_gradient2(low = 'light green',mid = 'light blue',high = 'red')+
      #scale_color_gradient2(values = c("red","blue","aquamarine1","aquamarine4",'blue1',"blue4"),name = 'Model(Level)',labels = c('FRK(Mix-Level)','Kriging(Mix-Level)','FRK','Kriging','FRK','Kriging'))+
      #ylim(-0.7,0.5)+
      #xlim(0,0.9)+
      ylab(~Q-R^2)+
      xlab("Morans'I")+
      #scale_color_manual(name ="Values", values = value, labels = c("* 0.1","* 0.3")) +
      #facet_grid(~ value)+
      #scale_colour_gradient2(low = 'red', mid = 'cyan',high='#F08080',midpoint = 0.7,name = 'Distance Power',breaks = c(0.01,0.5,1.0,1.5,2.0),labels = c('0.01','0.5','1.0','1.5','2.0'))+
      theme_bw()+
      theme(panel.grid =element_blank())+
      theme(axis.title.x = element_text(size = 7,color="black"),axis.title.y = element_text(size = 7,color="black"))+
      theme(axis.text.x = element_text(size = 7,color="black"),axis.text.y = element_text(size = 7,color="black"))+
      theme(legend.key.size = unit(7, "pt"),legend.title = element_text(size = 7,color="black"),axis.text.y = element_text(size = 7,color="black"))+
      theme(legend.text = element_text(size = 6,color="black"),axis.text.y = element_text(size = 6,color="black"))+
      theme(legend.position = c(0.6,0.1))
    g2



    ggplot2::ggsave('D:/0 job/0 q-vpc/results_simulation/MI_diff_mat_point.pdf',device = 'pdf',plot = g2,width = 16,units = 'cm',height = 15,dpi = 1200)



    library(tidyr)
    library(data.table)
    qr2.diff.mat = (simuList_q_csv-simuList_r2_csv)/simuList_r2_csv
    data1 = data.frame(x=rowMeans(simuList_MI_csv[2:98,]))
    data1 = cbind(data1,qr2.diff.mat[2:98,])
    setDT(data1)
    data1 = melt(data1,id='x')

    data1 = data1[data1$value<0.5 & data1$value > -1,]
    g2 = ggplot(data1, aes(x=x,y=value)) +
      geom_point(size = 0.05)+
      #geom_line(size = 0.5) +
      #scale_fill_continuous(type = "viridis") +
      #geom_abline(intercept=0.75,slope=0,color = 'black')+
      geom_hline(yintercept = 0,color = 'red4',size = 0.3,lty = "dashed")+
      geom_hline(yintercept = -1,color = 'red4',size = 0.3,lty = "dashed")+
      #geom_smooth(method = 'loess',formula = y~x,size = 0.2,se= T,na.rm = T,level = 0.90)+
      scale_color_gradient2(low = 'light green',mid = 'light blue',high = 'red')+
      #scale_color_gradient2(values = c("red","blue","aquamarine1","aquamarine4",'blue1',"blue4"),name = 'Model(Level)',labels = c('FRK(Mix-Level)','Kriging(Mix-Level)','FRK','Kriging','FRK','Kriging'))+
      #ylim(-0.7,0.5)+
      #xlim(0,0.9)+
      ylab(~R^2)+
      xlab('Alpha')+
      #scale_color_manual(name ="Values", values = value, labels = c("* 0.1","* 0.3")) +
      #facet_grid(~ value)+
      #scale_colour_gradient2(low = 'red', mid = 'cyan',high='#F08080',midpoint = 0.7,name = 'Distance Power',breaks = c(0.01,0.5,1.0,1.5,2.0),labels = c('0.01','0.5','1.0','1.5','2.0'))+
      theme_bw()+
      theme(panel.grid =element_blank())+
      theme(axis.title.x = element_text(size = 7,color="black"),axis.title.y = element_text(size = 7,color="black"))+
      theme(axis.text.x = element_text(size = 7,color="black"),axis.text.y = element_text(size = 7,color="black"))+
      theme(legend.key.size = unit(7, "pt"),legend.title = element_text(size = 7,color="black"),axis.text.y = element_text(size = 7,color="black"))+
      theme(legend.text = element_text(size = 6,color="black"),axis.text.y = element_text(size = 6,color="black"))+
      theme(legend.position = c(0.8,0.1))
    g2
    ggplot2::ggsave('D:/0 job/0 q-vpc/results_simulation/qr2_diff_mat.jpg',plot = g2,width = 8,units = 'cm',height = 7.5,dpi = 1200)


  }






}else{

  #write.csv(x = mat.list,file = 'D:/0 job/0 q-vpc/data/simulation/mat.list.random.csv')
  write.csv(x = simuList_q_all,file = 'D:/0 job/0 q-vpc/data/simulation/0 simuList_q_sar_1000_simg0.1.csv')
  write.csv(x = simuList_r2_all,file = 'D:/0 job/0 q-vpc/data/simulation/0 simuList_r2_sar_1000_simg0.1.csv')
  write.csv(x = simuList_rho_all,file = 'D:/0 job/0 q-vpc/data/simulation/0 simuList_rho_sar_1000_simg0.1.csv')
  write.csv(x = simuList_MI_all,file = 'D:/0 job/0 q-vpc/data/simulation/0 simuList_MI_sar_1000_simg0.1.csv')

}


plot(seq(0.1,0.9,0.1),(rowMeans(simuList_q[1:9,1:10])-rowMeans(simuList_r2[1:9,1:10]))/rowMeans(simuList_r2[1:9,1:10])*100,
     xlab = 'lambda',ylab = 'SEM_R^2 - Q',
     main = 'Difference between SEM_R^2 and Q',
     #xlim = c(0,1),ylim = c(0,-100)
     )

trans.Lambda = -exp(4.5*seq(0.01,0.99,0.01))
plot((trans.Lambda-min(trans.Lambda))/(max(trans.Lambda)-min(trans.Lambda)),(rowMeans(simuList_q[2:98,])-rowMeans(simuList_r2[2:98,]))/rowMeans(simuList_r2[2:98,])*100
     ,xlab = 'lambda',ylab = 'SEM_R^2 - Q',
     main = 'Difference between SEM_R^2 and Q',
     #xlim = c(0,1),ylim = c(0,-80)
     )
  abline(h = 0, col = "gray")
  abline(v = 0, col = "gray")
  abline(a = 0 , b =1 , col = "red")

tt



















































y_hat2 = (x.mat %*% sar_model$coefficients)
R_squared2 = 1 - var((y-y_hat2))/var(y)

y_line = y - sar_model$rho*mat.list%*%y
y_hat_line =(x.mat %*% sar_model$coefficients)
R_squared_line = 1 - var((y_line-y_hat_line))/var(y_line)

######## y_hat
simuList_q = matrix(0,10,100)
simuList_r2 = matrix(0,10,100)
simuList_r2lm = matrix(0,10,100)
simuList_realr2 = matrix(0,10,100)
simuList_morans = matrix(0,10,100)

expl = 60
sim_size = nrow(DiseaseData_sf)


for (rho in seq(0.1,0.9,0.1)) {
  for (simi in 1:50) {
    sim_size = nrow(DiseaseData_sf)
    x <- sample(1:3, size = sim_size, replace = TRUE)
    #解释度对应的方差
    #sigma_e = 0.25 * var(x)/(expl/100)-0.25 * var(x)


    mu = rnorm(sim_size,0,0.3)
    x.mat = model.matrix(lm(mu ~ factor(x)))

    IW_solve = solve(diag(1,sim_size)-rho*mat.list)
    y <- IW_solve%*%(1 + 0.5 * x.mat[,2]+ 0.7 * x.mat[,3] + mu)


    #总方差
    sst <-  var(IW_solve%*%(0.5 * x.mat[,2]+ 0.7 * x.mat[,3])) + var(IW_solve%*%mu)

    #真实解释度
    real_R2 = var(IW_solve%*%(0.5 * x.mat[,2]+ 0.7 * x.mat[,3])) / sst
    simuList_realr2[as.integer(rho*10),simi] = real_R2
    #hist(y)

    sim_df = data.frame(y = y)
    sim_df$x = x
    q = factor_detector('y','x',sim_df)
    simuList_q[as.integer(rho*10),simi] = unlist(q)[1]

    listW = nb2listw(nb.list ,zero.policy = T)


    semData = data.frame(x0 = x.mat[,1],x1 = x.mat[,2],x2 = x.mat[,3])
    semData$y = y
    names(semData)
    library(HSAR)
    se_model <- errorsarlm(y ~x1+x2, data = semData,
                           listw = listW, zero.policy = T)

    y_hat =(x.mat %*% se_model$coefficients)

    R_squared = 1 - sum((y-y_hat)^2)/sum((y-mean(y))^2)
    simuList_r2[as.integer(rho*10),simi] = R_squared
  }
  print(rho)
}

plot(rowMeans(simuList_q[1:9,]),rowMeans(simuList_r2[1:9,]))
plot(seq(0.1,0.9,0.1),rowMeans(simuList_realr2[1:9,]),xlab = 'lambda',ylab = 'Real_R^2',main = 'Real_R^2')
plot(seq(0.1,0.9,0.1),rowMeans(simuList_realr2[1:9,])-rowMeans(simuList_q[1:9,]),xlab = 'lambda',ylab = 'Real_R^2 - Q',main = 'Difference between Real_R^2 and Q')
plot(seq(0.01,0.99,0.01),rowMeans(simuList_realr2[1:9,])-rowMeans(simuList_r2[1:9,]),xlab = 'lambda',ylab = 'Real_R^2 - SEM_R^2',main = 'Difference between Real_R^2 and SEM_R^2')
plot(seq(0.1,0.9,0.1),abs(rowMeans(simuList_q[1:9,])-rowMeans(simuList_r2[1:9,]))/abs(rowMeans(simuList_r2[1:9,])),xlab = 'lambda',ylab = 'SAR_R^2 - Q (%)',main = 'Difference between SAR_R^2 and Q (%)')


realR2_array = array(unlist(simuList_realr2[1:9,]))
q_array = array(unlist(simuList_q[1:9,]))
realR2_array_remove = realR2_array[realR2_array>0]
q_array_remove = q_array[realR2_array>0]

plot(realR2_array_remove,q_array_remove,xlab = 'Real R^2',ylab = 'Q',xlim = c(0,0.6),ylim = c(0,0.6),main = 'Correlation between Real_R^2 and Q')
abline(h = 0, col = "gray")
abline(v = 0, col = "gray")
abline(a = 0 , b = 1, col = "red")
mean(abs(realR2_array_remove-q_array_remove))



realR2_array = array(unlist(simuList_realr2[1:9,]))
r2_array = array(unlist(simuList_r2[1:9,]))
realR2_array_remove = realR2_array[realR2_array>0]
r2_array_remove = r2_array[realR2_array>0]

plot(realR2_array_remove,r2_array_remove,xlab = 'Real R^2',ylab = 'SEM R^2',xlim = c(0,0.6),ylim = c(0,0.6),main = 'Correlation between Real_R^2 and SEM_R^2')
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




