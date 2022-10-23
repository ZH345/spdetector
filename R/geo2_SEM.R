
### Geodector paper
## Proof of the concept
# the simplest case with one variable --- discretized continuous variable

# cut a continuous variable into three categories with probability/proportion equal to the vector p
p <- c(0.3,0.3,0.4)

x <- rmultinom(3000,1,p)

rowSums(x) / 3000

n <- 10000
x <- sample(1:3,n,replace = TRUE)
d <- data.frame(y =rep(0,n), x = x)
x_mat<- model.matrix( ~ factor(x)- 1, data = d)
head(x_mat)
b <- c(1,1.5, 1)
#b_c <-c(0.5, b)
x_c <- x_mat

# make the first dummy variable (category) as base category
x_c[,1] <- rep(1,n)

d$y <- x_mat[] %*% b + rnorm(n, mean = 0, sd = sqrt(0.1))


## deleting the intercept term yields meaningless (at least noncomparable model fit)
#m0 <- lm(y ~ factor(x) - 1, data = d)
#summary(m0)

m1 <- lm(y ~ x_mat[,-1], data = d)
summary(m1)

## check the r2 from matrix calculations
b <- coefficients(m1)
new_b <- b
new_b[-1] <- new_b[-1] + new_b[1]

## calculate the variance explained the model
Var_exp <- as.numeric(t(new_b) %*% var(x_mat) %*% new_b)
Var_exp / var(d$y)

p <- table(d$x) / n
p * (1-p)

########## Equal to the OLS r2 estimate
d$f1 <- factor(d$x)

library(geodetector)

factor_detector("y","f1",tabledata = d)



#### Section 2 --- scenario where a spatial autocorrelation process was presented
####               in the true data generation process

### the SAR process case
library(geodetector)
library(spdep)
library(spatialreg)

data("DiseaseData_shp")
str(DiseaseData_shp, max.level = 2)
head(DiseaseData_shp@data)


nc <- st_read(system.file("shapes/sids.shp", package="spData")[1], quiet=TRUE)
st_crs(nc) <- "+proj=longlat +datum=NAD27"
row.names(nc) <- as.character(nc$FIPSNO)
plot(st_geometry(nc))

# neighbourhood structure and weights matrix
nb_nc <- spdep::poly2nb(nc)
lisw_nc <- spdep::nb2listw(nb_nc, style = "W", zero.policy = TRUE)


### a simulation experiment
n <- nrow(nc)
b <- c(1,1.5, 1) # note that this yields an medium degree of explanation power
sigma <- 0.1 # again, this relates to the calculation of R2
rho <- seq(0.05, 0.95, by = 0.05)
Nsim <- 1000

# save results
diff.res <- matrix(nrow = Nsim, ncol = length(rho))
true.R2  <- matrix(nrow = Nsim, ncol = length(rho))
true.R2.pvalue <- matrix(nrow = Nsim, ncol = length(rho))
R2.sar = matrix(nrow = Nsim, ncol = length(rho))
Q.mat  <- matrix(nrow = Nsim, ncol = length(rho))
Q.pvalue <- matrix(nrow = Nsim, ncol = length(rho))
MI.mat  <- matrix(nrow = Nsim, ncol = length(rho))
adjust.q = matrix(nrow = Nsim, ncol = length(rho))
adjust.q.sar = matrix(nrow = Nsim, ncol = length(rho))

for (i in 1:length(rho)) {
  rho.temp <- rho[i]
  inv.temp <- invIrW(lisw_nc,rho.temp)

  # for each rho value, we do 1000 simulations
  for (j in 1:Nsim) {
    # generate x
    x.temp <- sample(1:3,n,replace = TRUE)
    d.temp <- data.frame(y =rep(0,n), x = x.temp)
    x_mat.temp <- model.matrix( ~ factor(x.temp), data = d.temp)
    mu = rnorm(n, 0, sqrt(sigma))
    d.temp$y <- x_mat.temp %*% b + inv.temp %*% ( mu)

    # model estimation
    #sar <- spautolm(y ~ x_mat.temp[,-1],data = d.temp, listw = lisw_nc, family = "SAR")
    #sem.fit <- errorsarlm(y ~ x_mat.temp[,-1],data = d.temp, listw = lisw_nc)
    #fit.y <- predict(sem.fit) # include trend = Xβ, and signal = λW y − λW Xβ
    #fit.y = x_mat.temp %*% sem.fit$coefficients

    #R2.true = 1-(var(d.temp$y - fit.y))/var(d.temp$y)
    #adjust.true = var(sar$rho * lag.listw(lisw_nc,d.temp$y))/var(d.temp$y)


    # check what is actually predicted
    #x_mat.temp.1 <- x_mat.temp
    #x_mat.temp.1[,1] <- rep(1,times=n)
    #inv.temp.0 <- invIrW(lisw_nc,sar$rho)
    #ff2 <- inv.temp.0 %*% (x_mat.temp.1 %*% sar$coefficients)
    #sum(abs(ff2 - as.numeric(ff1)))
    ## OK, use the above

    ### The second thought: true R2 should be based on two parameters of rho and betas
    d.temp$y_tilde <- d.temp$y - rho.temp * lag.listw(lisw_nc,d.temp$y) + rho.temp * lag.listw(lisw_nc,(x_mat.temp %*% b))

    #adjust.true = var(rho.temp * lag.listw(lisw_nc,d.temp$y))/var(d.temp$y)

    # run an OLS model
    m.temp <- lm(y_tilde ~ x_mat.temp[,-1], data = d.temp)
    R2.true <- summary(m.temp)$r.squared

    mi = moran.test(as.vector(d.temp$y), lisw_nc)
    mi.num= as.numeric(unlist(mi$estimate[1]))


    #SST = var(d.temp$y)
    #SSE = var(inv.temp %*% (mu))
    #R2.true = 1- SSE/SST

    #x_mat.temp.1 <- x_mat.temp
    #x_mat.temp.1[,1] <- rep(1,times=n)
    #fit.y <- inv.temp %*% (x_mat.temp.1 %*% b)
    #R2.true <- cor(fit.y,d.temp$y)^2
    # Q-statistic
    d.temp$f1 <- factor(d.temp$x)
    Q.res <- factor_detector("y","f1",tabledata = d.temp)
    Q.temp <- as.numeric(Q.res[[1]][1])
    Q.p = as.numeric(Q.res[[1]][2])
    # save the difference
    adjust.q[j,i] = adjust.true
    R2.sar[j,i] = as.numeric(R.sar)
    MI.mat[j,i] = mi.num
    Q.pvalue[j,i] = Q.p
    diff <- (Q.temp - R2.true) / R2.true * 100
    diff.res[j,i] <- diff
    true.R2[j,i]  <- R2.true
    Q.mat[j,i]    <- Q.temp
  }
  print(i)
  # the end of the j loop
}

#QP = matrix(0,nrow = Nsim, ncol = length(rho))
#RP = matrix(0,nrow = Nsim, ncol = length(rho))
#QP[Q.pvalue>0.1] = NA
#RP[true.R2.pvalue>0.01] = NA
# diff.res_act = diff.res + QP + RP

#write.csv(diff.res,'D:/0 job/0 q-vpc/results_simulation/0 data.diff.res.csv')
#write.csv(MI.mat,'D:/0 job/0 q-vpc/results_simulation/0 data.MI.mat.csv')
#write.csv(true.R2,'D:/0 job/0 q-vpc/results_simulation/0 data.true.R2.csv')
#write.csv(Q.mat,'D:/0 job/0 q-vpc/results_simulation/0 data.Q.mat.csv')

#diff.res = read.csv('D:/0 job/0 q-vpc/results_simulation/0 data.diff.res.csv',header = T,sep = ',',row.names = 1)
#MI.mat = read.csv('D:/0 job/0 q-vpc/results_simulation/0 data.MI.mat.csv',header = T,sep = ',',row.names = 1)
#true.R2 = read.csv('D:/0 job/0 q-vpc/results_simulation/0 data.true.R2.csv',header = T,sep = ',',row.names = 1)
#true.R2 = read.csv('D:/0 job/0 q-vpc/results_simulation/0 data.Q.mat.csv',header = T,sep = ',',row.names = 1)


diff.res_act = diff.res


summary.diff <- apply(diff.res_act,2,quantile,probs=c(0.5,0.025,0.975),na.rm = T)
summary.r2 <- apply(true.R2,2,quantile,probs=c(0.5,0.025,0.975),na.rm = T)
#summary.diff <- apply(diff.res,2,quantile,probs=seq(0,1,0.05),na.rm = T)

plot(rho, summary.diff[1,])

library(ggplot2)

plot.data <- data.frame(rho = rho,diff50 = summary.diff[1,]/100,
                        diff25 = summary.diff[2,]/100,
                        diff975 = summary.diff[3,]/100)

## Fit a nonlinear least square model
diff.model <- nls(diff50 ~ b * rho^c, data = plot.data,
                  start=list(b = -0.5, c = 2))

summary(diff.model)
xx <- predict(diff.model)
cor(xx,plot.data$diff50)^2

plot.data$diff50.fit <- xx


## mapping
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


names(nc)


ggplot2::ggsave('D:/BaiduSyncdisk/0 job/0 q-vpc/results_simulation/SEM_NC_rho_diff.pdf',device = 'pdf',plot = g1,width = 16,units = 'cm',height = 15,dpi = 1200)

## Fit a nonlinear least square model for Morans'I and diff
plot.data$MI = colMeans(MI.mat)
plot(plot.data$MI, summary.diff[1,])
diff.model <- nls(diff50 ~ b * MI, data = plot.data,
                  start=list(b = -1))


summary(diff.model)
xx <- predict(diff.model)
cor(xx,plot.data$diff50)^2

plot.data$diff50.MI.fit <- xx
g2 = ggplot(data = plot.data, aes(x = MI)) +
  geom_linerange(aes(ymin=diff25, ymax=diff975),show.legend = T)+
  geom_line(aes(y=diff50.MI.fit,color="Fitted curve"),size=0.5) +
  geom_point(aes(y=diff50,shape="Diff with lower limit (2.5%) and upper limit (97.5%)"),color = 'black')+
  #geom_line(aes(y=diff25.fit,color="Lower limit (2.5%)"),linetype="dashed",size=0.1) +
  #geom_line(aes(y=diff975.fit,color="Upper limit (97.5%)"),linetype="dashed",size=0.1) +
  ylim(-0.9,0.01)+
  xlim(0.05,1)+
  scale_x_continuous(breaks = scales::breaks_width(0.2, offset = 0.05))+
  #geom_hline(yintercept = 0,color = 'red4',size = 0.3,lty = "dashed")+
  #geom_hline(yintercept = -1,color = 'red4',size = 0.3,lty = "dashed")+
  theme_bw()+
  theme(legend.position = c(0.4,0.3))+
  theme(legend.background = element_rect(fill = "transparent",colour = NA))

g2
ggplot2::ggsave('D:/BaiduSyncdisk/0 job/0 q-vpc/results_simulation/SEM_NC_MI_diff_mat.pdf',device = 'pdf',plot = g2,width = 16,units = 'cm',height = 15,dpi = 1200)







### an example of how to deal with 2 categorical variables#################################################################################
### the SAR process case
library(geodetector)
library(spdep)
library(spatialreg)

data("DiseaseData_shp")
str(DiseaseData_shp, max.level = 2)
head(DiseaseData_shp@data)


nc <- st_read(system.file("shapes/sids.shp", package="spData")[1], quiet=TRUE)
st_crs(nc) <- "+proj=longlat +datum=NAD27"
row.names(nc) <- as.character(nc$FIPSNO)
plot(st_geometry(nc))

# neighbourhood structure and weights matrix
nb_nc <- spdep::poly2nb(nc)
lisw_nc <- spdep::nb2listw(nb_nc, style = "W", zero.policy = TRUE)
plot(st_geometry(nc),border = 'grey')
plot( st_centroid(nc$geometry), add=TRUE )
plot(lisw_nc, st_centroid(nc$geometry), add= TRUE )
n <- nrow(nc)

b = c(1, 2.5, 1.5, 1.9, 2.3, 2.8, 1.7, 3, 3.2)
sigma <- 1.5 # again, this relates to the calculation of R2
rho <- seq(0.05, 0.95, by = 0.05)
Nsim <- 1000

# save results
diff.res <- matrix(nrow = Nsim, ncol = length(rho))
true.R2  <- matrix(nrow = Nsim, ncol = length(rho))
true.R2.pvalue <- matrix(nrow = Nsim, ncol = length(rho))
Q.mat  <- matrix(nrow = Nsim, ncol = length(rho))
Q.pvalue <- matrix(nrow = Nsim, ncol = length(rho))
MI.mat  <- matrix(nrow = Nsim, ncol = length(rho))
for (i in 1:length(rho)) {
  rho.temp <- rho[i]
  inv.temp <- invIrW(lisw_nc,rho.temp)

  # for each rho value, we do 1000 simulations
  for (j in 1:Nsim) {
    # generate x
    x1.temp <- sample(1:3,n,replace = TRUE)
    x2.temp <- sample(1:3,n,replace = TRUE)
    d.temp <- data.frame(x1=x1.temp, x2=x2.temp)
    x_mat.temp <- model.matrix( ~ factor(x1)*factor(x2), data = d.temp)

    mu = rnorm(n, 0, sqrt(sigma))
    d.temp$y <- x_mat.temp %*% b + inv.temp %*% (mu)

    ### The second thought: true R2 should be based on two parameters of rho and betas
    d.temp$y_tilde <- d.temp$y - rho.temp * lag.listw(lisw_nc,d.temp$y)
    # run an OLS model
    m.temp <- lm(y_tilde ~ x_mat.temp[,-1], data = d.temp)
    R2.true <- summary(m.temp)$r.squared

    mi = moran.test(as.vector(d.temp$y), lisw_nc)
    mi.num= as.numeric(unlist(mi$estimate[1]))


    #SST = var(inv.temp %*% (x_mat.temp %*% b)) + var(inv.temp %*% (mu))
    #SSE = var(inv.temp %*% (mu))
    #R2.true = 1- SSE/SST

    #x_mat.temp.1 <- x_mat.temp
    #x_mat.temp.1[,1] <- rep(1,times=n)
    #fit.y <- inv.temp %*% (x_mat.temp.1 %*% b)
    #R2.true <- cor(fit.y,d.temp$y)^2
    # Q-statistic
    d.temp$x1 <- ddd$x1
    d.temp$x2 <- ddd$x2
    Q.res <- interaction_detector("y",c("x1","x2"),tabledata = d.temp)
    Q.temp <- as.numeric(Q.res[1,3])
    #Q.p = as.numeric(Q.res[[1]][2])
    # save the difference
    MI.mat[j,i] = mi.num
    #Q.pvalue[j,i] = Q.p
    diff <- (Q.temp - R2.true) / R2.true * 100
    diff.res[j,i] <- diff
    true.R2[j,i]  <- R2.true
    Q.mat[j,i]    <- Q.temp
  }
  print(i)
  # the end of the j loop
}

#write.csv(diff.res,'D:/0 job/0 q-vpc/results_simulation/0 data.diff.res_inter.csv')
#write.csv(MI.mat,'D:/0 job/0 q-vpc/results_simulation/0 data.MI.mat_inter.csv')
#write.csv(true.R2,'D:/0 job/0 q-vpc/results_simulation/0 data.true.R2_inter.csv')
#write.csv(Q.mat,'D:/0 job/0 q-vpc/results_simulation/0 data.Q.mat_inter.csv')

diff.res = read.csv('D:/0 job/0 q-vpc/results_simulation/0 data.diff.res_inter.csv',header = T,sep = ',',row.names = 1)
MI.mat = read.csv('D:/0 job/0 q-vpc/results_simulation/0 data.MI.mat_inter.csv',header = T,sep = ',',row.names = 1)
true.R2 = read.csv('D:/0 job/0 q-vpc/results_simulation/0 data.true.R2_inter.csv',header = T,sep = ',',row.names = 1)
true.R2 = read.csv('D:/0 job/0 q-vpc/results_simulation/0 data.Q.mat_inter.csv',header = T,sep = ',',row.names = 1)

diff.res_act = diff.res

summary.diff <- apply(diff.res_act,2,quantile,probs=c(0.5,0.025,0.975),na.rm = T)
summary.r2 <- apply(true.R2,2,quantile,probs=c(0.5,0.025,0.975),na.rm = T)
#summary.diff <- apply(diff.res,2,quantile,probs=seq(0,1,0.05),na.rm = T)

plot(rho, summary.diff[1,])

library(ggplot2)

plot.data <- data.frame(rho = rho,diff50 = summary.diff[1,]/100,
                        diff25 = summary.diff[2,]/100,
                        diff975 = summary.diff[3,]/100)

## Fit a nonlinear least square model
diff.model <- nls(diff50 ~ b * rho^c, data = plot.data,
                  start=list(b = -0.5, c = 2))

summary(diff.model)
xx <- predict(diff.model)
cor(xx,plot.data$diff50)^2

plot.data$diff50.fit <- xx

## mapping
g3 = ggplot(data = plot.data, aes(x = rho)) +
  geom_linerange(aes(ymin=diff25, ymax=diff975),show.legend = T)+
  geom_line(aes(y=diff50.fit,color="Fitted curve"),size=0.5) +
  geom_point(aes(y=diff50,shape="Diff with lower limit (2.5%) and upper limit (97.5%)"),color = 'black')+
  #geom_line(aes(y=diff25.fit,color="Lower limit (2.5%)"),linetype="dashed",size=0.1) +
  #geom_line(aes(y=diff975.fit,color="Upper limit (97.5%)"),linetype="dashed",size=0.1) +
  ylim(-0.8,0.01)+
  xlim(0.05,1)+
  scale_x_continuous(breaks = scales::breaks_width(0.2, offset = 0.05))+
  #geom_hline(yintercept = 0,color = 'red4',size = 0.3,lty = "dashed")+
  #geom_hline(yintercept = -1,color = 'red4',size = 0.3,lty = "dashed")+
  theme_bw()+
  theme(legend.position = c(0.4,0.3))+
  theme(legend.background = element_rect(fill = "transparent",colour = NA))

g3

ggplot2::ggsave('D:/0 job/0 q-vpc/results_simulation/NC_inter_rho_diff_mat.pdf',device = 'pdf',plot = g3,width = 16,units = 'cm',height = 15,dpi = 1200)

## Fit a nonlinear least square model for Morans'I and diff
plot.data$MI = colMeans(MI.mat)
diff.model <- nls(diff50 ~ b * MI^c, data = plot.data,
                  start=list(b = -0.5, c = 2))

summary(diff.model)
xx <- predict(diff.model)
cor(xx,plot.data$diff50)^2

plot.data$diff50.MI.fit <- xx
g4 = ggplot(data = plot.data, aes(x = MI)) +
  geom_linerange(aes(ymin=diff25, ymax=diff975),show.legend = T)+
  geom_line(aes(y=diff50.MI.fit,color="Fitted curve"),size=0.5) +
  geom_point(aes(y=diff50,shape="Diff with lower limit (2.5%) and upper limit (97.5%)"),color = 'black')+
  #geom_line(aes(y=diff25.fit,color="Lower limit (2.5%)"),linetype="dashed",size=0.1) +
  #geom_line(aes(y=diff975.fit,color="Upper limit (97.5%)"),linetype="dashed",size=0.1) +
  ylim(-0.8,0.01)+
  xlim(0.05,1)+
  scale_x_continuous(breaks = scales::breaks_width(0.2, offset = 0.05))+
  #geom_hline(yintercept = 0,color = 'red4',size = 0.3,lty = "dashed")+
  #geom_hline(yintercept = -1,color = 'red4',size = 0.3,lty = "dashed")+
  theme_bw()+
  theme(legend.position = c(0.4,0.3))+
  theme(legend.background = element_rect(fill = "transparent",colour = NA))

g4
ggplot2::ggsave('D:/0 job/0 q-vpc/results_simulation/NC_inter_MI_diff_mat.pdf',device = 'pdf',plot = g4,width = 16,units = 'cm',height = 15,dpi = 1200)

#




