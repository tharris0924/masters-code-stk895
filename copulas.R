library(ggplot2)
library(WVPlots)
library('btld')
dat1 <- data.frame(x=runif(1000))
ggplot(dat1, aes(x))+geom_histogram(aes(y=..density..),breaks=seq(0,1,0.1), position = "identity", lwd=0.2)
# qnorm is the icdf in a matter of speaking
trans_dat1 <- data.frame(norm_x=qnorm(dat1$x))
trans_dat2 <- cbind(dat1, trans_dat1)
trans_dat2
k<-ScatterHist(frame = trans_dat2, xvar = "x",
            yvar= "norm_x",
            title = 'Uniform Transformed to a Normal',
            contour = F,
            smoothmethod = 'none',
            point_color = "#006d2c",
            point_alpha = 0.3,
            hist_color = "Orange")
# as seen by the plot, the icdf streetchs the outer regions of the uniform to yield a normal

dat1 <- data.frame(Original=runif(1000))
ggplot(dat1, aes(Original))+geom_histogram(aes(y=..density..),breaks=seq(0,1,0.1), position = "identity", lwd=0.2)
# qnorm is the icdf in a matter of speaking
# this can be done for a beta distribution
trans_dat1 <- data.frame(Transformed=qbeta(dat1$Original, shape1 = 10, shape2 = 3))
trans_dat2 <- cbind(dat1, trans_dat1)
trans_dat2
stats::acf2AR()

ScatterHist(frame = trans_dat2,
            xvar = "Original",
            yvar= "Transformed",
            title = 'Uniform transformed to a Beta',
            contour = F,
            smoothmethod = 'none',
            point_color = "#006d2c",
            point_alpha = 0.3,
            hist_color = "Orange")
# or a gumbel
library(evd)
u<-dat1$Original
gumbel <- data.frame(transformed = qgumbel(u))
new <- cbind(gumbel,dat1)
str(new)

ScatterHist(frame = new,
            xvar = "Original",
            yvar= "transformed",
            title = 'Uniform transformed to a Gumbel',
            contour = F,
            smoothmethod = 'none',
            point_color = "#006d2c",
            point_alpha = 0.3,
            hist_color = "Orange")

# transform from distro to unif, apply inverse of icdf i.e. cdf
x_trans_trans <- pgumbel(gumbel)
colnames(gumbel)<- 'Original'
colnames(x_trans_trans)<-'Transformed'
new <- cbind(gumbel,x_trans_trans)
ScatterHist(frame = new,
            xvar = "Original",
            yvar= "Transformed",
            title = 'Gumbel transformed back to a Uniform',
            contour = F,
            smoothmethod = 'none',
            point_color = "#006d2c",
            point_alpha = 0.3,
            hist_color = "Orange")

# this is called the probability integral transform
# generate a gaussian copula
#step 1 generate a multivariate gaussian
### copula plots
library(btld)
library(copula)
library('RColorBrewer')
# Condtional Gaussian Copula Plots
library('dplyr')
n<-1000
norm.cop <- normalCopula(0.5, dim = 2) # intialize copula function
u <- rCopula(1000, norm.cop)
splom2(u)
u1 <- 0.05
U <- cCopula(cbind(u1, runif(1000)), copula = norm.cop, inverse = TRUE)
splom2(qnorm(U), cex=0.2)
plot(sort(U[,2]),sort(u[,2]),ylab = quote(U[2]), xlab= quote(U[1]), main="Conditional Gaussian Copula Plots")
## A large u_1\
png("gaussian1.png", width = 800, height = 500)
par(mfrow = c(1,1))

par(xpd = T, mar = par()$mar + c(0,0,0,7))
u1 <- 0.95
u2<-0.75
u3<-0.25
u4<-0.05
n<-1000
U <- cCopula(cbind(u1, runif(n)), copula = norm.cop, inverse = TRUE)
U2<-cCopula(cbind(u2,runif(n)), copula=norm.cop,inverse=TRUE)
U3 <- cCopula(cbind(u3, runif(n)), copula = norm.cop, inverse = TRUE)
U4<-cCopula(cbind(u4,runif(n)), copula=norm.cop,inverse=TRUE)

# plot(sort(U[,2]), ylab = quote(U[2]), xlab= quote(U[1]), main="Conditional Gaussian Copula Plots")
plot(sort(U[,2]),sort(u[,2]), ylab = quote(u[2]), xlab= quote(u[1]), main="Conditional Gaussian Copula Plots", type="l", lty=1, col="Purple", cex.lab=1.2)
lines(sort(U2[,2]),sort(u[,2]), col="Green", type="l", xlab = " ", ylab = " ", lty=2)
lines(sort(U3[,2]),sort(u[,2]), col="Red", type="l", xlab = " ", ylab = " ", lty=3)
lines(sort(U4[,2]),sort(u[,2]), col="Blue", type="l", xlab = " ", ylab = " ", lty=4)

legend(1.10, 0.25, c('U2|U1 = 0.95','U2|U1 =0.75','U2|U1=0.25','U2|U1=0.05'),
       col = c("Purple", "Green", "Red", "Blue"),
       text.col = "black", lty = c(1,2,3,4), cex=0.8,lwd=1)
par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()

png('pcopula.png', width = 800, height=500)

par(xpd = T, mar = par()$mar + c(0,0,0,7))
cols<-brewer.pal(9,'YlOrRd')
contour(norm.cop, pCopula, main='Gaussian Copula Contour Density Plots',col = cols, cex.lab=1.2)
abline(v=c(0.95,0.75, 0.25, 0.05),col=c( 'purple', 'green', 'red','blue'))

legend(1.10, 0.25, c('U2|U1 = 0.95','U2|U1 =0.75','U2|U1=0.25','U2|U1=0.05'),
       col = c("Purple", "Green", "Red", "Blue"),
       text.col = "black", lty = c(1,2,3,4), cex=0.8,lwd=1)
par(mar=c(5, 4, 4, 2) + 0.1)
dev.off()
library('ggplot2')
library('PerformanceAnalytics')

## one d-vector =^= 1-row matrix, works too :
dCopula(c(0.5, 0.5), norm.cop)
pCopula(c(0.5, 0.5), norm.cop)
plot(u)
plot(dCopula(u, norm.cop)) # return density functions
pCopula(u, norm.cop) # return cum-prob functions
persp (norm.cop, dCopula,main="3D Contour Surface Plot of Gaussian Copula", col="green",zlab = "Density")

ScatterHist(frame = mvn,
            xvar = "x1",
            yvar= "x2",
            title = 'Gaussian Multivariate',
            contour = T,
            smoothmethod = 'none',
            point_alpha = 0.3,
            point_color ='Red',
            hist_color = 'DarkOrange',
            contour_color = 'Green')

#step 2: uniformify the marginals using probability integral transform
x_unif_1 <- qbtld(mvn$x1, theta1 = 0.3,theta2 = 0.8,3,3)
x_unif_2 <- qbtld(mvn$x2, theta1 = 0.3,theta2 = 0.8,3,3)
mvbtld<- data.frame(y1=x_unif_1,y2=x_unif_2)
colnames(mvbtld)<-c('y1','y2')

ScatterHist(frame =mvbtld,
            xvar = "y1",
            yvar= "y2",
            title = 'MVBTLD',
            contour = T,
            smoothmethod = 'none',
            point_alpha = 0.3,
            point_color ='Red',
            hist_color = 'DarkOrange',
            contour_color = 'Green')
#step 3, transform into multivariate distribution using icdfs
m1 <-qgumbel(x_unif_1)
m2 <-qbeta(x_unif_2,shape1 = 10, shape2 = 2)
trans <- data.frame(MaximumRiverLevel=m1, ProbFlooding=m2)

ScatterHist(frame = trans,
            xvar = "MaximumRiverLevel",
            yvar= "ProbFlooding",
            title = 'Multivariate Distribution with Gumbel and Beta Marginals using Gaussian Copula',
            contour = T,
            smoothmethod = 'none',
            point_alpha = 0.3,
            point_color ='White',
            hist_color = 'Blue',
            contour_color = 'DarkBlue')


######################################################################
 qbltd <- function (x, theta1, theta2, alpha1, alpha2){

  alpha0<-(1-alpha1*theta1/2 -alpha2*(1-theta2)/2)/(theta2-theta1)
   print(alpha0)
  u_df <- data.frame(x=x)

  # specify valid bounds
  lower_df <- u_df[u_df$x<=alpha1*theta1/2,]
  middle_df<-u_df[(alpha1*theta1/2)<u_df$x & u_df$x<=(alpha1*theta1/2+alpha0*(theta2-theta1)),]
  upper_df <- u_df[u_df$x>(alpha1*theta1/2+alpha0*(theta2-theta1)),]

   # generate values
  lower_df<-sqrt((2*lower_df*theta1)/alpha1)
  middle_df<-(middle_df-(alpha1*theta1)/2)/alpha0 + theta1
  upper_df<-1-sqrt((2*(1-upper_df)*(1-theta2))/alpha2)

  df <- c(lower_df,middle_df,upper_df)
  qs <- data.frame(GenRvs=df)
   return(qs)

 }

rvs <-qbltd(runif(1000), theta1 = 0.3, theta2 = 0.8, alpha1 = 3, alpha2 =3)
plt <- ggplot(data=rvs, aes(x=GenRvs))
plt + geom_density(fill="Blue") + ggtitle('Density plot of Bimodal Triangular-Linked Distribution')
plt +geom_histogram(color="Black", fill="Blue", bins=30)+ggtitle("Histogram of Bimodal Triangular-Linked Distribution")

pbtld <- function (vec, alpha1, alpha2, theta1, theta2){
     u<-vec
    alpha0<-(1-alpha1*theta1/2 -alpha2*(1-theta2)/2)/(theta2-theta1)
    u_df <- data.frame(unif=u)

   # specify valid bounds
    lower_df <- u_df[u_df$unif<=theta1,]
    middle_df<-u_df[theta1<u_df$unif & u_df$unif<=theta2,]
    upper_df <- u_df[u_df$unif>theta2,]
    inputs<-c(lower_df,middle_df,upper_df)
  # generate values
    lower_cdf<-alpha1*lower_df^2/(2*theta1)
    middle_cdf<-alpha0*(middle_df-theta1)+alpha1*theta1/2
    upper_cdf<-1-(alpha2*(1-upper_df)^2)/(2*(1-theta2))

    cdf <- c(lower_cdf,middle_cdf,upper_cdf)
    density <- data.frame(values=inputs, Denisty=cdf)
    return(density)
  }
cdf <- pbtld(vec=rvs$GenRvs, theta1 = 0.3, theta2 = 0.8, alpha1 = 3, alpha2 =3)
# bivariate case
library(mvtnorm)
means <- c(0,0)
sigma <- matrix(c(1,0.5,0.5,1), ncol = 2)

x <- rmvnorm(n=100000, mean = means,sigma = sigma)
mvn<-data.frame(x)
colnames(mvn)<-c('x1','x2')

ScatterHist(frame = mvn,
            xvar = "x1",
            yvar= "x2",
            title = 'Gaussian Multivariate',
            contour = T,
            smoothmethod = 'none',
            point_alpha = 0.3,
            point_color ='Red',
            hist_color = 'DarkOrange',
            contour_color = 'Green')

#step 2: uniformify the marginals using probability integral transform
x_unif_1 <- pnorm(mvn$x1)
x_unif_2 <- pnorm(mvn$x2)
cop<- data.frame(y1=x_unif_1,y2=x_unif_2)

ScatterHist(frame = cop,
            xvar = "y1",
            yvar= "y2",
            title = 'Gaussian Copula',
            contour = T,
            smoothmethod = 'none',
            point_alpha = 0.3,
            point_color ='Red',
            hist_color = 'DarkOrange',
            contour_color = 'Green')

m1 <-qbltd(x=x_unif_1, alpha1 = 5, alpha2 = 1, theta1=0.3, theta2 = 0.7)
print(m1)
m2 <-qbltd(x =x_unif_2,alpha1 = 1, alpha2 = 5, theta1 = 0.3,theta2=0.7)
hist(m2$GenRvs)
trans <- cbind(m1, m2)
colnames(trans)<-c('x1', 'x2')

ScatterHist(frame = trans,
            xvar = "x1",
            yvar= "x2",
            title = 'Multivariate Bimodal Triangular-Linked Distribution using Gaussian Copula',
            contour = T,
            smoothmethod = 'none',
            point_alpha = 0.3,
            point_color ='Red',
            hist_color = 'Blue',
            contour_color = 'DarkBlue')
x=x_unif_1
library(rayshader)
library(ggplot2)

gg <- ggplot(data=mvbtld, aes(x=X1, y=X2))+geom_hex(bins=20, size=0) + scale_fill_viridis_c(option='C')+ggtitle('Hexagonal Plot of Multivariate Bimodal Triangular-Linked Distribution ')
gg
plot_gg(gg, width = 4, height = 4, scale = 300, multicore = T)
#####################################################################

library(mvtnorm)
means <- c(0,0,0)
sigma <- matrix(c(1,0.5, 0.5,0.5,1, 0.5, 0.5, 0.5, 1), ncol = 3)
print(sigma)

x <- rmvnorm(n=100000, mean = means,sigma = sigma)
mvn<-data.frame(x)
colnames(mvn)<-c('x1','x2', 'x3')
x_unif_1 <- pnorm(mvn$x1)
x_unif_2 <- pnorm(mvn$x2)
x_unif_3 <- pnorm(mvn$x3)
cop<- data.frame(y1=x_unif_1,y2=x_unif_2, y3=x_unif_3)

ggplot(data=cop, aes(x=y1, y=y2))+geom_hex(bins=20, size=0, color="Black") + scale_fill_viridis_c(option='C') + ggtitle("Gaussian Copula Hexagonal Plot")

p<-ggplot(data=cop, aes(x=y1, y=y2))+geom_hex(bins=20, size=0, color="Black") + scale_fill_viridis_c(option='C')
plot_gg(p, width = 4, height = 4, scale = 300, multicore = T)
show(p)

m1 <-qbltd(x=x_unif_1, alpha1 = 1, alpha2 = 5, theta1=0.3, theta2 = 0.7)
m2 <-qbltd(x =x_unif_2,alpha1 = 3, alpha2 = 3, theta1 = 0.3,theta2=0.7)
m3 <-qbltd(x =x_unif_3,alpha1 = 5, alpha2 = 1, theta1 = 0.3,theta2=0.7)

hist(m2$GenRvs)
trans <- cbind(m1, m2, m3)
colnames(trans)<-c('x1', 'x2','x3')
cor(trans)
gg<- ggplot(data=trans, aes(x=x2, y=x1))+geom_hex(bins=20, size=0, color="Black") + scale_fill_viridis_c(option='C')+xlab('x1')+ylab('y1')#+ggtitle('Multivariate Bimodal Triangular-Linked Distribution')#+xlab('x1')+ylab('y1')
gg
plot_gg(gg, width = 4, height = 4, scale = 300, multicore = T)
ScatterHist(frame = trans,
            xvar = "x1",
            yvar= "x2",
            title = 'Multivariate Bimodal Triangular-Linked Distribution using Gaussian Copula',
            contour = T,
            smoothmethod = 'none',
            point_alpha = 0.3,
            point_color ='Red',
            hist_color = 'Blue',
            contour_color = 'DarkBlue')

cor(trans, method = 'pearson')
options(warn = defaultW)
par(new=T)
plot(ecdf(trans$x1), xlim=c(0,1), ylim=c(0,1),
     main = "H-plot", col='Red')
lines(ecdf(trans$x2),col='Green', main=' ')
lines(ecdf(trans$x3),col='Blue', main=' ')
legend(0,0.8,c('BTLD(0.3,0.7,1,5)','BTLD(0.3,0.7,3,3)','BTLD(0.3,0.7,5,1)'))
plot(ecdf(trans$x3))
ecdf(trans$x1)

ecdf <- rep(1/100000,100000)
ecdf<- cumsum(ecdf)
plot(sort(trans$x1),ecdf,type = "l", lty="dashed",col="Blue", main = "H-Plot", xlab = "X", ylab="Cumulative Density")
par(new=TRUE)
lines(sort(trans$x2),ecdf, col="Red", type="l", xlab = " ", ylab = " ")
lines(sort(trans$x3),ecdf, col="Green", type="l", xlab = " ", ylab = " ")
legend(0.7, 0.5, c('BTLD(0.3,0.7,1,5)','BTLD(0.3,0.7,3,3)','BTLD(0.3,0.7,5,1)'), col = c("Blue", "Red", "Green"),
       text.col = "black", lty = c("dashed", "solid", "dotted"),
       merge = TRUE, bg = "gray90", text.font = "small")
trans$x1


rbtld <- function (size, alpha, theta){
  u<-runif(size)
  alpha0<-(1-(alpha[1]*theta[1]/2) -(alpha[2]*(1-theta[2]))/2)/(theta[2]-theta[1])
  A<-(alpha[1]*theta[1]/2)
  r<-which(u<A)
  x1<-sqrt((2*theta[1]*u[r])/alpha[1])
  C<-A+alpha0*(theta[2]-theta[1])
  r<-which(u>C)
  x3<-1-sqrt((2*(1-u[r])*(1-theta[2]))/alpha[2])
  ls<-A<u & u<C
  r<-which(ls)
  x2<-theta[1]+(u[r]-A)/alpha0
  x<-c(x1,x2,x3)
}

dbtld <- function (xvalues, alpha,theta){
  theta<-as.character(theta)
   alpha0<-(1-(alpha[1]*theta[1]/2) -(alpha[2]*(1-theta[2]))/2)/(theta[2]-theta[1])
   r <- which(xvalues<theta[1])
   x1<-alpha1*xvalues[r]/(theta[1])
   r<-which(xvalues>theta[2])
   x3<-(alpha2*(1-xvalues[r]))/(1-theta2)
   n <- length(xvalues)-length(x1)-length(x2)
   x2<-rep(alpha0, times = n)
   joined_pdf <- c(x1,x2,x3)
   return(joined_pdf)
}

lower_x<-0
upper_x<-1
width<-abs(upper_x-lower_x) # distance on x plane
x1<-runif(1000, lower_x, upper_x) # uniform random variables on x plane
z<-dbtld(x1,  alpha = c(5,1), theta =c(0.3, 0.7)) # function values giving values on z plane
height <- mean(z) # distance on z plane
area_mcmc <- width*height
area_mcmc

x<-rbtld(1000, c(2, 2), c(0.2500,  0.800))
d<-dbtld(sort(xvalues),c(2, 2), c(0.2500,  0.800))
hist(x, freq = F)
plot(d, type="l")

alpha<-c(3, 3)
theta<-c(0.2500,  0.800)

c0=(alpha[1]*theta[1]/2+alpha[2]*(1-theta[2])/2-1)/(theta[1]-theta[2])

# rbtld
alpha0<-(1-(alpha[1]*theta[1]/2) -(alpha[2]*(1-theta[2]))/2)/(theta[2]-theta[1])
u<-runif(1000)
A<-(alpha[1]*theta[1]/2)
r<-which(u<A)
x1<-sqrt((2*theta[1]*u[r])/alpha[1])
C<-A+c0*(theta[2]-theta[1])
# C<-(alpha[1]*theta[1])/2+alpha0*(theta[2]-theta[1]) #(this is an area)
r<-which(u>C)
x3<-1-sqrt((2*(1-u[r])*(1-theta[2]))/alpha[2])
ls<-A<u & u<C
r<-which(ls)
?which
x2<-theta[1]+(u[r]-A)/c0
x<-c(x1,x2,x3)


hist(rbtld(1000,alpha,theta))
hist(x, freq = FALSE)
library('btld')

hist(tlrnvd(1000,alpha,theta))

alpha0<-(1-(alpha[1]*theta[1]/2) -(alpha[2]*(1-theta[2]))/2)/(theta[2]-theta[1])
df <- data.frame(unif=xvalues)
lower_df <- df[df$unif<theta[1],]
n1<-length(lower_df)
upper_df <- df[df$unif>theta[2],]
n2<-length(upper_df)
rem<-1000-(n2+n1)
middle_df<-runif(rem,0.25,0.75)

lower_pdf<-alpha[1]*lower_df/(theta[1])
middle_pdf<-runif(rem,0.25,0.75)
upper_pdf<-(alpha[2]*(1-upper_df))/(1-theta[2])

joined_pdf <- c(lower_pdf,middle_pdf,upper_pdf)

plot(sort(xvalues),joined_pdf,type="l")
library('btld')
alpha<-c(3, 3)
theta<-c(0.2500,  0.800)
rvs<-rbtld(1000, alpha = alpha, theta = theta)
pdf<-dbtld(rvs,alpha,theta)
hist(rvs,probability = T, breaks = 50)
ecdf<-empcdf(rvs)
dens<-density(rvs, kernel="triangular")
?density
lines(dens$x,dens$y)
cbind(rvs,dens$y)
# how to find the maximum in a given interval
p1<-optimize(approxfun(dens$x,dens$y), interval = c(0,0.5), maximum=TRUE)
p2<-optimize(approxfun(dens$x,dens$y), interval = c(0.5,1), maximum=TRUE)
?optimize
ecdf<-empcdf(rvs)
sort.unit.gene<-sort(rvs)
ecdf<-ecdf[round(ecdf[,1], digits = 2)==round(p1$maximum, digits=2),]
ft1<-colMeans(ecdf)

ecdf<-empcdf(rvs)
ecdf<-ecdf[round(ecdf[,1], digits = 2)==round(p2$maximum, digits=2),]
ft2<-colMeans(ecdf)

find.mode<-function (df){
  intervals<-list(int1=c(0,0.5), int2=c(0.5,1))
  mode<-NULL
  for(i in intervals) {
    dens <- density(df)
    peak <- optimize(approxfun(dens$x, dens$y), maximum = TRUE, interval = i)
    ecdf <- empcdf(df)
    ecdf <- ecdf[round(ecdf[, 1], digits = 2) == round(peak$maximum, digits = 2),]
    mode <-rbind(mode,colMeans(ecdf))
  }
  return(mode)
}


find.mode(rvs)
library('goftest')
alpha<-btld_scales(c(0.1,0.5),c(0.4,0.9))
ks.test(rvs,"pbtld", alpha=alpha,theta=theta)
cvm.test(rvs,"pbtld",alpha=alpha,theta=theta,estimated=TRUE)
?cvm.test

rvs<-rtridist(1000,0.25)
ks.test(rvs,"ptridist", theta=0.25)
cvm.test(rvs,"ptridist",theta=0.25)

# multivariate normal mixture testing
library('btld')
alpha<-matrix(c(3,3,3,3,3,3),nrow=3,ncol=2, byrow=TRUE)
theta<-matrix(c(0.3,0.7,0.3,0.7,0.3,0.7),ncol=2,nrow=3,byrow=TRUE)
df<-rmvbtld(1000, alpha=alpha, theta=theta, sigma=0.5,dim=3)

scatter3D(df[,3], df[,2], df[,1], pch = 18,  theta = 120, phi = -40,
          main = "MVTL Scatterplot", xlab = "X1",
          ylab ="X2", zlab = "X3", ticktype="detailed",bty="g")
plot(df)

par(mfrow = c(1, 1))
# an illustration
alpha<-matrix(c(5,1,3,3,1,5),nrow=3,ncol=2, byrow=TRUE)
theta<-matrix(c(0.3,0.7,0.3,0.7,0.3,0.7),ncol=2,nrow=3,byrow=TRUE)
df<-rmvbtld(1000, alpha=alpha, theta=theta, sigma=0.5,dim=3)
scatter3D(df[,3], df[,2], df[,1], pch = 18,  theta = 20, phi = 40,
          main = "Simulated Example: MVTL Scatterplot All Marginals", xlab = "X1",
          ylab ="X2", zlab = "X3", ticktype="detailed",bty="g")

df<-as.data.frame(df)
colnames(df)<-c('x1','x2','x3')
par(mfrow=c(1,1))
ScatterHist(frame = df,
            xvar = "x1",
            yvar= "x3",
            title = 'MVTL pre-softmax normalization',
            contour = T,
            smoothmethod = 'none',
            point_alpha = 0.3,
            point_color ='Red',
            hist_color = 'Blue',
            contour_color = 'Green')

df<-btld::softmax(df, exp=F, base=1000)

ScatterHist(frame = df,
            xvar = "x1",
            yvar= "x3",
            title = 'MVTL post-softmax normalization',
            contour = T,
            smoothmethod = 'none',
            point_alpha = 0.3,
            point_color ='Red',
            hist_color = 'Blue',
            contour_color = 'Green')
library('btld')
alpha<-matrix(c(1,5,3,3,5,1),nrow=3,ncol=2, byrow=TRUE)
theta<-matrix(c(0.3,0.7,0.15,0.8,0.3,0.7),ncol=2,nrow=3,byrow=TRUE)
df<-rmvbtld(1000, alpha=alpha, theta=theta, sigma=0.5,dim=3)
library('plot3D')
df<-as.data.frame(df)
colnames(df)<-c('x1','x2','x3')
scatter3D(df[,3], df[,2], df[,1], pch = 18,  theta = 20, phi = 40,
          main = "Simulated Example: MVTL Scatterplot Identical Marginals ", xlab = latex2exp::TeX('$X_1$'),
          ylab =latex2exp::TeX('$X_2$'), zlab = latex2exp::TeX('$X_3$'), ticktype="detailed",bty="g")

x1<-rbtld(1000,alpha[1,],theta[1,])
x2<-rbtld(1000,alpha[2,],theta[2,])
par(mfrow=c(1,1))
ScatterHist(frame = df,
            xvar = "x1",
            yvar= "x2",
            title = 'Multivariate Triangular-Linked Distribution using Gaussian Copula with Identical Marginals',
            contour = T,
            smoothmethod = 'none',
            point_alpha = 0.3,
            point_color ='Red',
            hist_color = 'Blue',
            contour_color = 'DarkBlue')
ind<-data.frame(x1=x1,x2=x2)
ScatterHist(frame = ind,
            xvar = "x1",
            yvar= "x2",
            title = 'Independant Marginals',
            contour = T,
            smoothmethod = 'none',
            point_alpha = 0.3,
            point_color ='Red',
            hist_color = 'Green',
            contour_color = 'DarkGreen')

library('MixAll')
?plot3D
lab_x1<-latex2exp::TeX('$X_1$')
lab_x2<-latex2exp::TeX("$X_2$")
lab_x3<-latex2exp::TeX("$X_3$")
# 3D kde contours
freqz <- with(data.frame(x1,x2), MASS::kde2d(x1, x2, n = 50))
with(freqz, plot_ly(x = x, y = y, z = z, type = "surface", colors = "YlOrRd"))

freqz <- with(data.frame(df$x1,df$x2), MASS::kde2d(df$x1, df$x2, n = 50))
with(freqz, plot_ly(x = x, y = y, z = z, type = "surface", colors = "YlGnBu"))


m <- ggplot(df, aes(x = x1, y = x2)) +
  geom_point(fill="Black")

# contour lines
m + geom_density_2d_filled(alpha=1)

m <- ggplot(ind, aes(x = x1, y = x2)) +
  geom_point(fill="Black")
# contour lines
m + geom_density_2d_filled(alpha=1)
?ggplot2
library('ggplot2')
image2D(df[,1:2])
em_func <- function(input_dat, plot=F,...){
  dat <- input_dat
  em <- clusterDiagGaussian(data = dat, models = "gaussian_pk_sjk", criterion = "AIC",
                            strategy = clusterStrategy(nbTry = 5, nbInit = 2,
                                                       initMethod = "class",
                                                       initAlgo = "EM", nbInitIteration = 20,
                                                       initEpsilon = 0.01,
                                                       nbShortRun = 2, shortRunAlgo = "EM",
                                                       nbShortIteration = 100,
                                                       shortEpsilon = 1e-04, longRunAlgo = "EM",
                                                       nbLongIteration = 1000,
                                                       longEpsilon = 1e-07),...)
  if(plot==T){
    plot(em)
  }
  return(list(MU=em@component@mean, SIGMA=em@component@sigma,
              PI1=em@pk[1], PI2=em@pk[2], AIC=em@criterion))
}

ems<-em_func(df, plot = T, nbCluster=4)
xtable::xtable(ems$MU, digits = 4)

