
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
b <- c(5,1.5, 1) # note that this yields an medium degree of explanation power
sigma <- 0.1 # again, this relates to the calculation of R2
rho <- seq(0.05, 0.95, by = 0.05)
Nsim <- 50

# save results
diff.res <- matrix(nrow = Nsim, ncol = length(rho))
true.R2  <- matrix(nrow = Nsim, ncol = length(rho))
true.R2.pvalue = matrix(nrow = Nsim, ncol = length(rho))
Q.mat  <- matrix(nrow = Nsim, ncol = length(rho))
Q.pvalue <- matrix(nrow = Nsim, ncol = length(rho))

for (i in 1:length(rho)) {
  rho.temp <- rho[i]
  inv.temp <- invIrW(lisw_nc,rho.temp)
  
  # for each rho value, we do 1000 simulations
  for (j in 1:Nsim) {
    # generate x
    x.temp <- sample(1:3,n,replace = TRUE)
    d.temp <- data.frame(y =rep(0,n), x = x.temp)
    x_mat.temp <- model.matrix( ~ factor(x.temp), data = d.temp)
    d.temp$y <- inv.temp %*% (x_mat.temp %*% b + rnorm(n, 0, sqrt(sigma)))
    
    
    if (T){
      sar.lag <- lagsarlm(y ~ x_mat.temp[,-1],data = d.temp, listw = lisw_nc)
      Wald1.p = summary(sar.lag)$Wald1$p.value
      
      fit.y = predict(sar.lag,listw = lisw_nc,pred.type = 'TC')
      #fit.y.lag = predict(sar.lag,listw = lisw_nc)

      
    }else{
      # model estimation 
      sar <- spautolm(y ~ x_mat.temp-1,data = d.temp, listw = lisw_nc, family = "SAR")
      summary(sar)
      fit.y <- sar$fit$fitted.values
      inv.temp_hat <- invIrW(lisw_nc,sar$lambda)
      mat.nc = listw2mat(lisw_nc)
      sar$lambda*mat.nc%*% d.temp$y
      error.sar = errorsarlm(y ~ x_mat.temp[,-1],data = d.temp, listw = lisw_nc)
      summary(error.sar)
      fit.y_hat = inv.temp_hat%*%( cbind(rep(1,n),x_mat.temp[,-1]) %*% sar$fit$coefficients)
      
    }
    R2.true <- 1 - crossprod(d.temp$y-fit.y) / crossprod(d.temp$y-mean(d.temp$y))
    
    # Q-statistic
    d.temp$f1 <- factor(d.temp$x)
    Q.res <- factor_detector("y","f1",tabledata = d.temp)
    Q.temp <- as.numeric(Q.res[[1]][1])
    Q.p = as.numeric(Q.res[[1]][2])
    
    # save the difference
    true.R2.pvalue[j,i] = Wald1.p
    Q.pvalue[j,i] = Q.p
    diff <- Q.temp - R2.true
    diff.res[j,i] <- diff
    true.R2[j,i]  <- R2.true
    Q.mat[j,i]    <- Q.temp
  }
  print(i)
 # the end of the j loop
}

QP = matrix(0,nrow = Nsim, ncol = length(rho))
RP = matrix(0,nrow = Nsim, ncol = length(rho))
QP[Q.pvalue>0.01] = NA
RP[true.R2.pvalue>0.01] = NA


diff.res = diff.res + QP + RP
summary.diff <- apply(diff.res,2,quantile,probs=c(0.5,0.025,0.975),na.rm = T)
plot(rho, summary.diff[1,])

library(ggplot2)

plot.data <- data.frame(rho = rho,diff50 = summary.diff[1,],
                        diff25 = summary.diff[2,],
                        diff975 = summary.diff[3,])

## Fit a nonlinear least square model
diff.model <- nls(diff50 ~ b * rho^c, data = plot.data, 
                  start=list(b = 0.5, c = -0.2))
summary(diff.model)
xx <- predict(diff.model)
cor(xx,plot.data$diff50)^2
plot.data$diff50.fit <- xx

ggplot(data = plot.data, aes(x = rho)) +
  geom_linerange(aes(ymin=diff25, ymax=diff975))+
  geom_point(aes(y=diff50)) +
  geom_line(aes(y=diff50.fit), col = "red") +
  theme_classic()
  


### an example of how to deal with 2 categorical variables
x1.temp <- sample(1:3,n,replace = TRUE)
x2.temp <- sample(1:3,n,replace = TRUE)
ddd <- data.frame(x1=x1.temp, x2=x2.temp)
x_mat.temp <- model.matrix( ~ factor(x1)*factor(x2)-1, data = ddd)
