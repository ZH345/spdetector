#install.packages("geodetector")
library(geodetector)
library(dplyr)
data(CollectData)
names(CollectData)

# single variable
length(CollectData$incidence)
var(CollectData$incidence[CollectData$elevation==1])
var(CollectData$incidence[CollectData$elevation==2])
var(CollectData$incidence[CollectData$elevation==3])
var(CollectData$incidence[CollectData$elevation==4])
var(CollectData$incidence[CollectData$elevation==5])
var(CollectData$incidence[CollectData$elevation==6])
var(CollectData$incidence[CollectData$elevation==7])
hist(CollectData$incidence)

model.1 = lm(incidence~factor(elevation),data = CollectData)
MoMa1 = model.matrix(model.1)
elevation = MoMa1[,1]+MoMa1[,2]*1+MoMa1[,3]*2+MoMa1[,4]*3+MoMa1[,5]*4+MoMa1[,6]*5+MoMa1[,7]*6
summary(model.1)
summary(model.1)
factor_detector("incidence","elevation",CollectData) 

model.3 = lm(CollectData$incidence~(CollectData$elevation))
MoMa3 = model.matrix(model.3)

summary(model.3)
plot(CollectData$incidence~CollectData$elevation)
incidence_groupmean = group_by(CollectData,elevation) %>% summarise_each(mean)
model.3 = lm(incidence_groupmean$incidence~(incidence_groupmean$elevation))
summary(model.3)


model.2 = lm(CollectData$incidence~CollectData$elevation-1)
MoMa2 = model.matrix(model.2)
summary(model.2)


#multi-variable
length(CollectData$incidence)
model.1 = lm(incidence~factor(elevation)*factor(watershed),data = CollectData)
MoMa2 = model.matrix(model.1)
#CollectData$elevation-elevation
summary(model.1)


model.3 = lm(incidence~elevation+watershed+soiltype,data = CollectData)
MoMa3 = model.matrix(model.3)
summary(model.3)

model.2 = lm(incidence~elevation+watershed-1,data = CollectData)
MoMa2 = model.matrix(model.2)
summary(model.2)
interaction_detector("incidence",c("elevation","watershed"),CollectData) 


#interaction_detector
model.4 = lm(incidence~factor(elevation)*factor(watershed),data = CollectData)
MoMa1 = model.matrix(model.4)
summary(model.4)
mean(CollectData$incidence[CollectData$elevation==1 & CollectData$watershed==1])


rowSums(MoMa1)

model.3 = lm(incidence~elevation+watershed-1,data = CollectData)
MoMa3 = model.matrix(model.3)
summary(model.3)

interaction_detector("incidence",c("elevation","watershed"),CollectData) 



#simulation1 - Homovariance
id = c(real_R2= NULL,r2_lm = NULL,r2_lm_factor = NULL,q = NULL)
simu_data = data.frame(id=1:99,real_R2= 0,r2_lm = 0,r2_lm_factor = 0,q = 0)

simuList_real_R2 = matrix(0,100,100)
simuList_r2_lm = matrix(0,100,100)
simuList_r2_lm_factor = matrix(0,100,100)
simuList_q = matrix(0,100,100)
for (expl in 99:1){

  for (simi in 1:100){
    
    sim_size = 100
    x <- sample(1:5, size = sim_size, replace = TRUE)
    #解释度对应的方差
    sigma_e = 0.25 * var(x)/(expl/100)-0.25 * var(x)
    e <- rnorm(sim_size,0,sqrt(sigma_e))
    var(e[x==1])
    var(e[x==2])
    var(e[x==3])
    var(e[x==4])
    var(e[x==5])
    hist(e)
    #总方差
    sst <- 0.25 * var(x) + var(e)
    
    #真实解释度
    real_R2 = 0.25 * var(x) / sst
    simuList_real_R2[expl,simi] = real_R2
    y <- 1 + 0.5 * x + e
    var(y[x==1])
    var(y[x==2])
    var(y[x==3])
    var(y[x==4])
    var(y[x==5])
    hist(y)
    
    r2_lm = summary(lm(y ~ x))$r.squared
    simuList_r2_lm[expl,simi] = r2_lm
    mean(summary(lm(y ~ x))$residuals)
    mean(y[x==2])
    1.16426+0.44695*2
    
    r2_lm_factor = summary(lm(y ~ factor(x)))$r.squared
    simuList_r2_lm_factor[expl,simi]= r2_lm_factor
    mean(summary(lm(y ~ factor(x)))$residuals)
    
    sim_df = as.data.frame(y)
    sim_df$x = x
    q = factor_detector('y','x',sim_df)
    simuList_q[expl,simi] = unlist(q)[1]
  }
  simu_data[expl,] = c(expl,mean(simuList_real_R2[expl,]),mean(simuList_r2_lm[expl,]),mean(simuList_r2_lm_factor[expl,]),mean(simuList_q[expl,]))
  print(expl)
}

plot(simu_data$real_R2- simu_data$q)
plot(simu_data$real_R2- simu_data$r2_lm)



#single 
x2 = (as.array(MoMa_sim_1[,2]))
solve(t(x2)%*%x2)%*%(t(x2)%*%y)

#multi
MoMa_sim_1

solve(t(MoMa_sim_1)%*%MoMa_sim_1)%*%(t(MoMa_sim_1)%*%y)


#simulation1 - interprection
id = c(real_R2= NULL,r2_lm = NULL,r2_lm_factor = NULL,q = NULL)
simu_data = data.frame(id=1:99,real_R2= 0,r2_lm = 0,r2_lm_factor = 0,q = 0)

simuList_real_R2 = matrix(0,100,100)
simuList_r2_lm = matrix(0,100,100)
simuList_r2_lm_factor = matrix(0,100,100)
simuList_q = matrix(0,100,100)
for (expl in 99:1){
  
  for (simi in 1:100){
    
    sim_size = 100
    x <- sample(1:5, size = sim_size, replace = TRUE)
    #解释度对应的方差
    sigma_e = 0.25 * var(x)/(expl/100)-0.25 * var(x)
    e <- rnorm(sim_size,0,sqrt(sigma_e))
    var(e[x==1])
    var(e[x==2])
    var(e[x==3])
    var(e[x==4])
    var(e[x==5])
    hist(e)
    #总方差
    sst <- 0.25 * var(x) + var(e)
    
    #真实解释度
    real_R2 = 0.25 * var(x) / sst
    simuList_real_R2[expl,simi] = real_R2
    y <- 1 + 0.5 * x + e
    var(y[x==1])
    var(y[x==2])
    var(y[x==3])
    var(y[x==4])
    var(y[x==5])
    hist(y)
    
    r2_lm = summary(lm(y ~ x))$r.squared
    simuList_r2_lm[expl,simi] = r2_lm
    
    mean(y[x==2])
    1.16426+0.44695*2
    mean(y[x==1])
    1.16426+0.44695*1
    mean(y[x==3])
    1.16426+0.44695*3
    
    r2_lm_factor = summary(lm(y ~ factor(x)))$r.squared
    simuList_r2_lm_factor[expl,simi]= r2_lm_factor
    
    sim_df = as.data.frame(y)
    sim_df$x = x
    q = factor_detector('y','x',sim_df)
    simuList_q[expl,simi] = unlist(q)[1]
  }
  simu_data[expl,] = c(expl,mean(simuList_real_R2[expl,]),mean(simuList_r2_lm[expl,]),mean(simuList_r2_lm_factor[expl,]),mean(simuList_q[expl,]))
  print(expl)
}

plot(simu_data$real_R2- simu_data$q)
plot(simu_data$real_R2- simu_data$r2_lm)

#simulation1 - Heteroscedasticity

id = c(real_R2= NULL,r2_lm = NULL,r2_lm_factor = NULL,q = NULL)
simu_data = data.frame(id=1:99,real_R2= 0,r2_lm = 0,r2_lm_factor = 0,q = 0)
simuList_real_R2 = matrix(0,100,100)
simuList_r2_lm = matrix(0,100,100)
simuList_r2_lm_factor = matrix(0,100,100)
simuList_q = matrix(0,100,100)
for (expl in 99:1){

  for (simi in 1:100){
    
    sim_size = 1000
    x <- sample(1:5, size = sim_size, replace = TRUE)
    #解释度对应的方差 具有异方差
    sigma_e = 0.25 * var(x)/(expl/100)-0.25 * var(x)
    e = 1:sim_size
    e[x==1] <- rnorm(length(e[x==1]),0,sqrt(sigma_e/2))
    e[x==2] <- rnorm(length(e[x==2]),0,sqrt(sigma_e/3))
    e[x==3] <- rnorm(length(e[x==3]),0,sqrt(sigma_e*2))
    e[x==4] <- rnorm(length(e[x==4]),0,sqrt(sigma_e/4))
    e[x==5] <- rnorm(length(e[x==5]),0,sqrt(sigma_e*3))
    var(e[x==1])
    var(e[x==2])
    var(e[x==3])
    var(e[x==4])
    var(e[x==5])
    #总方差
    sst <- 0.25 * var(x) + var(e)
    
    #真实解释度
    real_R2 = 0.25 * var(x) / sst
    simuList_real_R2[expl,simi] = real_R2
    y <- 1 + 0.5 * x + e
    
    r2_lm = summary(lm(y ~ x))$r.squared
    simuList_r2_lm[expl,simi] = r2_lm
    r2_lm_factor = summary(lm(y ~ factor(x)))$r.squared
    simuList_r2_lm_factor[expl,simi]= r2_lm_factor
    
    sim_df = as.data.frame(y)
    sim_df$x = x
    q = factor_detector('y','x',sim_df)
    simuList_q[expl,simi] = unlist(q)[1]
  }
  simu_data[expl,] = c(expl,mean(simuList_real_R2[expl,]),mean(simuList_r2_lm[expl,]),mean(simuList_r2_lm_factor[expl,]),mean(simuList_q[expl,]))
  print(expl)
}
plot(simu_data$real_R2- simu_data$q)
plot(simu_data$real_R2- simu_data$r2_lm)

#single 
x2 = (as.array(MoMa_sim_1[,2]))
solve(t(x2)%*%x2)%*%(t(x2)%*%y)

#multi
MoMa_sim_1

solve(t(MoMa_sim_1)%*%MoMa_sim_1)%*%(t(MoMa_sim_1)%*%y)




#simulation1 - Nonlinear

id = c(real_R2= NULL,r2_lm = NULL,r2_lm_factor = NULL,q = NULL)
simu_data = data.frame(id=1:99,real_R2= 0,r2_lm = 0,r2_lm_factor = 0,q = 0)
simuList_real_R2 = matrix(0,100,100)
simuList_r2_lm = matrix(0,100,100)
simuList_r2_lm_factor = matrix(0,100,100)
simuList_q = matrix(0,100,100)
for (expl in 99:1){
  
  for (simi in 1:100){
    
    sim_size = 1000
    x <- sample(1:5, size = sim_size, replace = TRUE)
    #解释度对应的方差 具有异方差
    sigma_e = 0.25 * var(x)/(expl/100)-0.25 * var(x)
    e = 1:sim_size
    e[x==1] <- rnorm(length(e[x==1]),0,sqrt(sigma_e/2))
    e[x==2] <- rnorm(length(e[x==2]),0,sqrt(sigma_e/3))
    e[x==3] <- rnorm(length(e[x==3]),0,sqrt(sigma_e*2))
    e[x==4] <- rnorm(length(e[x==4]),0,sqrt(sigma_e/4))
    e[x==5] <- rnorm(length(e[x==5]),0,sqrt(sigma_e*3))
    var(e[x==1])
    var(e[x==2])
    var(e[x==3])
    var(e[x==4])
    var(e[x==5])
    #总方差
    sst <- 0.25 * var(x) + var(e)
    
    #真实解释度
    real_R2 = 0.25 * var(x) / sst
    simuList_real_R2[expl,simi] = real_R2
    y <- exp(1 + 0.5 * x + e)
    ln_y = log(y)
    r2_lm = summary(lm(y ~ x))$r.squared
    simuList_r2_lm[expl,simi] = r2_lm
    r2_lm_factor = summary(lm(y ~ factor(x)))$r.squared
    simuList_r2_lm_factor[expl,simi]= r2_lm_factor
    
    sim_df = as.data.frame(y)
    sim_df$x = x
    q = factor_detector('y','x',sim_df)
    simuList_q[expl,simi] = unlist(q)[1]
  }
  simu_data[expl,] = c(expl,mean(simuList_real_R2[expl,]),mean(simuList_r2_lm[expl,]),mean(simuList_r2_lm_factor[expl,]),mean(simuList_q[expl,]))
  print(expl)
}
plot(simu_data$real_R2- simu_data$q)
plot(simu_data$real_R2- simu_data$r2_lm)




#simulation1 - Weak explanatory power at individual scale and strong explanatory power at high scale

id = c(real_R2= NULL,r2_lm = NULL,r2_lm_factor = NULL,q = NULL)
simu_data = data.frame(id=1:99,real_R2= 0,r2_lm = 0,r2_lm_factor = 0,q = 0)
simuList_real_R2 = matrix(0,100,100)
simuList_r2_lm = matrix(0,100,100)
simuList_r2_lm_factor = matrix(0,100,100)
simuList_q = matrix(0,100,100)

for (expl in 99:1){
  
  for (simi in 1:100){
    
    sim_size = 1000
    x <- sample(1:5, size = sim_size, replace = TRUE)
    #解释度对应的方差 具有异方差
    sigma_e = 0.25 * var(x)/(expl/100)-0.25 * var(x)
    e = 1:sim_size
    e <- rnorm(length(e),0,sqrt(sigma_e)*3)

    
    e[x==1] <- rnorm(length(e[x==1]),0,sqrt(sigma_e/2))
    e[x==2] <- rnorm(length(e[x==2]),0,sqrt(sigma_e/3))
    e[x==3] <- rnorm(length(e[x==3]),0,sqrt(sigma_e/1))
    e[x==4] <- rnorm(length(e[x==4]),0,sqrt(sigma_e/4))
    e[x==5] <- rnorm(length(e[x==5]),0,sqrt(sigma_e/5))
    var(e[x==1])
    var(e[x==2])
    var(e[x==3])
    var(e[x==4])
    var(e[x==5])
    #总方差
    sst <- 0.25 * var(x) + var(e)
    
    #真实解释度
    real_R2 = 0.25 * var(x) / sst
    simuList_real_R2[expl,simi] = real_R2
    y <- 1 + 0.5 * x + e
    
    r2_lm = summary(lm(y ~ x))$r.squared
    simuList_r2_lm[expl,simi] = r2_lm
    r2_lm_factor = summary(lm(y ~ factor(x)))$r.squared
    simuList_r2_lm_factor[expl,simi]= r2_lm_factor
    
    sim_df = as.data.frame(y)
    sim_df$x = x
    q = factor_detector('y','x',sim_df)
    simuList_q[expl,simi] = unlist(q)[1]
  }
  simu_data[expl,] = c(expl,mean(simuList_real_R2[expl,]),mean(simuList_r2_lm[expl,]),mean(simuList_r2_lm_factor[expl,]),mean(simuList_q[expl,]))
  print(expl)
}
plot(simu_data$real_R2- simu_data$q)
plot(simu_data$real_R2- simu_data$r2_lm)





#simulation1 - X and y do not match

id = c(real_R2= NULL,r2_lm = NULL,r2_lm_factor = NULL,q = NULL)
simu_data = data.frame(id=1:99,real_R2= 0,r2_lm = 0,r2_lm_factor = 0,q = 0)
simuList_real_R2 = matrix(0,100,100)
simuList_r2_lm = matrix(0,100,100)
simuList_r2_lm_factor = matrix(0,100,100)
simuList_q = matrix(0,100,100)
for (expl in 99:1){
  
  for (simi in 1:100){
    
    sim_size = 1000
    x <- sample(1:5, size = sim_size, replace = TRUE)
    #解释度对应的方差 具有异方差
    sigma_e = 0.25 * var(x)/(expl/100)-0.25 * var(x)
    e = 1:sim_size
    e <- rnorm(length(e),0,sqrt(sigma_e))
    
    
    var(e[x==1])
    var(e[x==2])
    var(e[x==3])
    var(e[x==4])
    var(e[x==5])
    
    
    #总方差
    sst <- 0.25 * var(x) + var(e)
    
    #真实解释度
    real_R2 = 0.25 * var(x) / sst
    simuList_real_R2[expl,simi] = real_R2
    y <- 1 + 0.5 * x + e
    
    x[x==3] = 0
    x[x==1] = 3
    x[x==0] = 1
    
    x[x==5] = 0
    x[x==2] = 5
    x[x==0] = 2
    
    r2_lm = summary(lm(y ~ x))$r.squared
    simuList_r2_lm[expl,simi] = r2_lm
    r2_lm_factor = summary(lm(y ~ factor(x)))$r.squared
    simuList_r2_lm_factor[expl,simi]= r2_lm_factor
    
    sim_df = as.data.frame(y)
    sim_df$x = x
    q = factor_detector('y','x',sim_df)
    simuList_q[expl,simi] = unlist(q)[1]
  }
  simu_data[expl,] = c(expl,mean(simuList_real_R2[expl,]),mean(simuList_r2_lm[expl,]),mean(simuList_r2_lm_factor[expl,]),mean(simuList_q[expl,]))
  print(expl)
}
plot(simu_data$real_R2- simu_data$q)
plot(simu_data$real_R2- simu_data$r2_lm)





#simulation1 - random coefficients

id = c(real_R2= NULL,r2_lm = NULL,r2_lm_factor = NULL,q = NULL)
simu_data = data.frame(id=1:99,real_R2= 0,r2_lm = 0,r2_lm_factor = 0,q = 0)
simuList_real_R2 = matrix(0,100,100)
simuList_r2_lm = matrix(0,100,100)
simuList_r2_lm_factor = matrix(0,100,100)
simuList_q = matrix(0,100,100)
for (expl in 99:1){
  
  for (simi in 1:100){
    
    sim_size = 1000
    x <- sample(1:5, size = sim_size, replace = TRUE)
    b = rnorm(sim_size,0.5,0.5)
    #解释度对应的方差 具有异方差
    sigma_e = var(b*x)/(expl/100)-var(b*x)
    e = 1:sim_size
    e <- rnorm(length(e),0,sqrt(sigma_e))
    
    
    var(e[x==1])
    var(e[x==2])
    var(e[x==3])
    var(e[x==4])
    var(e[x==5])
    
    
    #总方差
    sst <- var(b*x) + var(e)
    
    #真实解释度
    real_R2 = var(b*x) / sst
    simuList_real_R2[expl,simi] = real_R2
    y <- 1 + b * x + e
    

    
    r2_lm = summary(lm(y ~ x))$r.squared
    simuList_r2_lm[expl,simi] = r2_lm
    r2_lm_factor = summary(lm(y ~ factor(x)))$r.squared
    simuList_r2_lm_factor[expl,simi]= r2_lm_factor
    
    sim_df = as.data.frame(y)
    sim_df$x = x
    q = factor_detector('y','x',sim_df)
    simuList_q[expl,simi] = unlist(q)[1]
  }
  simu_data[expl,] = c(expl,mean(simuList_real_R2[expl,]),mean(simuList_r2_lm[expl,]),mean(simuList_r2_lm_factor[expl,]),mean(simuList_q[expl,]))
  print(expl)
}
plot(simu_data$real_R2- simu_data$q)
plot(simu_data$real_R2- simu_data$r2_lm)



#simulation1 - Discrete continuous variable


simu_data = data.frame(id=1:99,real_R2= 0,r2_lm = 0,r2_lm_break = 0, r2_lm_factor = 0,q = 0)

simuList_real_R2 = matrix(0,100,100)
simuList_r2_lm = matrix(0,100,100)
simuList_r2_lm_factor = matrix(0,100,100)
simuList_q = matrix(0,100,100)
simuList_r2_lm_break = matrix(0,100,100)
for (expl in 99:1){
  
  for (simi in 1:100){
    
    sim_size = 100
    x <- runif(n = sim_size,min =  1, max = 100)
    hist(x)
    #解释度对应的方差
    sigma_e = 0.25 * var(x)/(expl/100)-0.25 * var(x)
    e <- rnorm(sim_size,0,sqrt(sigma_e))
    var(e[x==1])
    var(e[x==2])
    var(e[x==3])
    var(e[x==4])
    var(e[x==5])
    hist(e)
    #总方差
    sst <- 0.25 * var(x) + var(e)
    
    #真实解释度
    real_R2 = 0.25 * var(x) / sst
    simuList_real_R2[expl,simi] = real_R2
    y <- 1 + 0.5 * x + e
    var(y[x==1])
    var(y[x==2])
    var(y[x==3])
    var(y[x==4])
    var(y[x==5])
    hist(y)
    library(classInt)
    x_break = classIntervals(var = x, n = 5, style = 'jenks', warnLargeN = F)
    x_break$var[x_break$var>=(x_break$brks[1]-1) & x_break$var<x_break$brks[2]] = 1
    x_break$var[x_break$var>=x_break$brks[2] & x_break$var<x_break$brks[3]] = 2
    x_break$var[x_break$var>=x_break$brks[3] & x_break$var<x_break$brks[4]] = 3
    x_break$var[x_break$var>=x_break$brks[4] & x_break$var<x_break$brks[5]] = 4
    x_break$var[x_break$var>=x_break$brks[5] & x_break$var<=(x_break$brks[6]+1)] = 5
    
    r2_lm = summary(lm(y ~ x))$r.squared
    simuList_r2_lm[expl,simi] = r2_lm
    
    r2_lm_break = summary(lm(y ~ x_break$var))$r.squared
    simuList_r2_lm_break[expl,simi] = r2_lm_break
    
    r2_lm_factor = summary(lm(y ~ factor(x_break$var)))$r.squared
    simuList_r2_lm_factor[expl,simi]= r2_lm_factor
    
    sim_df = as.data.frame(y)
    sim_df$x = x_break$var
    q = factor_detector('y','x',sim_df)
    simuList_q[expl,simi] = unlist(q)[1]
  }
  simu_data[expl,] = c(expl,mean(simuList_real_R2[expl,]),mean(simuList_r2_lm[expl,]),mean(simuList_r2_lm_break[expl,]),mean(simuList_r2_lm_factor[expl,]),mean(simuList_q[expl,]))
  print(expl)
}

plot(simu_data$real_R2- simu_data$q)
plot(simu_data$real_R2- simu_data$r2_lm)













