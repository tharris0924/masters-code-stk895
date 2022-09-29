# Title     : Testing of btld
# Objective : Generate plots for jmva paper
# Created by: Tristan Harris
# Created on: 7/5/2021
library(ucminf)
library(caret)
library(btld)
library(copula)
library(WVPlots)
library(ggplot2)
library(MGLM)
library(plot3D)
library(latex2exp)
library(xtable)
library(LICORS)
library('mclust')
library('readr')

# intialize copula function
u <- rCopula(1000, norm.cop)
hist(u)# generate 100 rvs from the cop
getSigma(norm.cop) # get the matrix
mvn<-data.frame(x1=u[,1], x2=u[,2])

## one d-vector =^= 1-row matrix, works too :
dCopula(c(0.5, 0.5), norm.cop)
pCopula(c(0.5, 0.5), norm.cop)
plot(u)
dCopula(u, norm.cop) # return density functions
pCopula(u, norm.cop) # return cum-prob functions
persp (norm.cop, dCopula)
contour(norm.cop, pCopula)

ScatterHist(frame = mvn,
            xvar = "x1",
            yvar= "x2",
            title = 'Gaussian Copula',
            contour = T,
            smoothmethod = 'none',
            point_alpha = 0.3,
            point_color ='Red',
            hist_color = 'DarkOrange',
            contour_color = 'Green')

#step 2: uniformify the marginals using probability integral transform
x_unif_1 <- qbtld(mvn$x1,theta = c(0.3,0.7), alpha = c(3,3))
x_unif_2 <- qbtld(mvn$x2, theta = c(0.3,0.7),alpha = c(3,3))
mvbtld<- data.frame(y1=x_unif_1,y2=x_unif_2)
colnames(mvbtld)<-c('y1','y2')

ScatterHist(frame =mvbtld,
            xvar = "y1",
            yvar= "y2",
            title = 'MVTLD',
            contour = T,
            smoothmethod = 'none',
            point_alpha = 0.3,
            point_color ='Red',
            hist_color = 'DarkOrange',
            contour_color = 'Green')
library('copula')
# 3 dimensional case
clay.cop <-claytonCopula(2, dim = 3) # intialize copula function
u <- rCopula(1000, clay.cop)
tau(clay.cop)
cor(u)
library('btld')
# generate 100 rvs from the cop

m1 <-qbtld(x=u[,1],alpha = c(3,3), theta = c(0.3,0.7))
m2 <-qbtld(x =u[,2],alpha = c(3,3), theta = c(0.3,0.7))
m3 <-qbtld(x =u[,2],alpha = c(3,3), theta =  c(0.3,0.7))

alpha<- matrix(c(1,5,3,3,5,1),nrow = 3,ncol = 2, byrow = TRUE)
theta<-t(matrix(c(0.3,0.7,0.3,0.7,0.3,0.7), nrow = 3,ncol = 2, byrow = TRUE))

trans <- cbind(m1, m2, m3)
cor(trans)
trans <- data.frame(trans)/rowSums(trans)
colnames(trans)<-c('x1', 'x2','x3')

ScatterHist(frame =trans,
            xvar = "x1",
            yvar= "x3",
            title = 'MVBTLD',
            contour = T,
            smoothmethod = 'none',
            point_alpha = 0.3,
            point_color ='Blue',
            hist_color = 'Green',
            contour_color = 'Red')

ecdf <- rep(1/nrow(trans),nrow(trans))
ecdf<- cumsum(ecdf)
plot(sort(trans$x1),ecdf,type = "l", lty="dashed",col="Blue", main = "H-Plot", xlab = "X", ylab="Cumulative Density")
par(new=TRUE)
plot(sort(trans$x2),ecdf, col="Red", type="l", ann=F,axes = F)
par(new=TRUE)
plot(sort(trans$x3),ecdf, col="Green", type="l", ann=F,axes = F)
x <- seq(0,1,0.1)
abline(h=x,v=x, col = "lightgray", lty = 3) #add a grid

legend(0.8, 0.2, c('BTLD(0.3,0.7,1,5)','BTLD(0.3,0.7,3,3)','BTLD(0.3,0.7,5,1)'), col = c("Blue", "Red", "Green"),
       text.col = "black", lty = c("dashed", "solid", "dotted"),
       merge = TRUE, bg = "gray90")

# library('strucchange')
# breakpoints(ecdf~sort(trans$x2), breaks = 2)


theta<- matrix(c(0.3,0.7), nrow = 3, ncol=2, byrow = T)
alpha<-matrix(c(3,3), nrow = 3,ncol=2,byrow = T)
d<-rmvbtld(1000,theta,alpha,0.5,3)
par(mfrow=c(1,2))
?exp

d<-100^(d)/rowSums(100^(d))

d<-exp()
d<-d/rowSums(d)
d<-data.frame(d)
colnames(d)<-c('x1','x2','x3')

ScatterHist(frame =d,
            xvar = "x1",
            yvar= "x3",
            title = 'MVBTLD',
            contour = T,
            smoothmethod = 'none',
            point_alpha = 0.3,
            point_color ='Blue',
            hist_color = 'Green',
            contour_color = 'Red')

scatter3D(trans$x1, trans$x2, trans$x3, pch = 18,  theta = 40, phi = 20,
          main = "MVBTL", xlab = "X1",
          ylab ="X2", zlab = "X3",ticktype="detailed",bty="g")
scatter3D(d[,1], d[,2], d[,3], pch = 18,  theta = 40, phi = 20,
        main = "MVBTL", xlab = "X1",
         ylab ="X2", zlab = "X3", ticktype="detailed")

scatter3D(d[,1], d[,2], d[,3], pch = 18,  theta = 40, phi = 20,
        main = "MVBTL", xlab = "X1",
          ylab ="X2", zlab = "X3", ticktype="detailed",bty="g")

scalex1  <-btld::btld_scales(c(0.27,0.7),c(0.17, 0.23))
scalex2  <-btld::btld_scales(c(0.3,0.7),c(0.43, 0.53))
scalex3  <-btld::btld_scales(c(0.3,0.7),c(0.73, 0.82))



# compositional mvbtld
# b<-matrix(c(.05, .05,.95, .95),nrow=2,2)
# c<-matrix(c(4.6, 31.62, 25.23, 3.07),nrow=2,2)
# rho<-as.numeric(matrix(c(1.0000, 0.5, 0.5, 1.0000),nrow=2,2,byrow = T))
# alpha<- matrix(c(1,5,3,3,5,1),nrow = 3,ncol = 2, byrow = TRUE)
# theta<-t(matrix(c(0.3,0.7,0.3,0.7,0.3,0.7), nrow = 2,ncol = 3, byrow = TRUE))
# y<-rmvbtld(1000,theta,alpha,0.5,3)
# scales <- y/rowSums(y)
# comp<-1-rowSums(y)

compdf <- data.frame(y, comp)
compdf <- compdf[!compdf$comp<0,]
resampled <- compdf[sample(seq_len(nrow(compdf)), 30),]

alpha<-matrix(c(5,1,3,3,1,5),nrow=3,ncol=2, byrow=TRUE)
theta<-matrix(c(0.3,0.7,0.3,0.7, 0.3, 0.7),ncol=3,nrow=2,byrow=TRUE)
mn<-rmvbtld(1000, alpha=alpha, theta=theta, sigma=0.5,dim=3)
scatter3D(mn[,1], mn[,2], mn[,3], pch = 18,  theta = 40, phi = 20,
          main = "MVBTL", xlab = "X1",
          ylab ="X2", zlab = "X3", ticktype="detailed",bty="g")

#so as to ensure randomness
alpha<-matrix(c(1,5,3,3,5,1),nrow=3,ncol=2, byrow=TRUE)
theta<-matrix(c(0.3,0.7,0.3,0.7,0.3,0.7),ncol=3,nrow=2,byrow=TRUE)
mn<-mn_cmvbtld(500, size=20, alpha=alpha, theta=theta, sigma=0.5,dim=3, base = 1000)
ys<-rmvbtld(500, alpha=alpha, theta=theta, sigma=0.5,dim=3)
p<-1000^(ys)/rowSums(1000^(ys))
mn<-data.frame(mn)



scatter3D(mn[,1], mn[,2], mn[,3], pch = 18,  theta = 120, phi = 40,
          main = "CMVTL-MN Simulation", xlab = "X1",
          ylab ="X2", zlab = "X3", ticktype="detailed",bty="g", colkey=FALSE)


x_c <- cut(p[,1], 20)
y_c <- cut(p[,2], 20)

##  Calculate joint counts at cut levels:
z <- table(x_c, y_c)
yy=1:20
xx=1:20

##  Plot as a 3D histogram:
hist3D(yy,xx, z=z,border="black", theta=40, phi=40, axes=TRUE,label=TRUE, zlim=c(0,15), ticktype='detailed')

set.seed(555)

n <- 500
d <- 4
m <- 20
y <- rdirmn(n, m,c(5,2,5))

scatter3D(y[,1], y[,2], y[,3], pch = 18,  theta = 60, phi = 40,
          main = "CDM", xlab = "X1",
          ylab ="X2", zlab = "X3", ticktype="detailed",bty="g")
#
x_c <- cut(y[,1], 20)
y_c <- cut(y[,2], 20)
#
# ##  Calculate joint counts at cut levels:
z <- table(x_c, y_c)
yy=1:20
xx=1:20

?hist3D
br = seq(0,10,by=1)
ranges = paste(head(br,-1), br[-1], sep=" - ")
freq_said   = hist(counts$said_norm, breaks=br, include.lowest=TRUE, plot=FALSE)
freq_showed = hist(counts$showed_norm, breaks=br, include.lowest=TRUE, plot=FALSE)
freq_kenn=hist(counts$kennedy_norm, breaks=br, include.lowest=TRUE, plot=FALSE)
df<-data.frame(range = ranges, Said = freq_said$counts/352, Showed=freq_showed$counts/352, Kennedy = freq_kenn$counts/352)
df2<-t(df)
colnames(df2)<-c('0.05', '0.15','0.25','0.35','0.45','0.55','0.65','0.75','0.85','0.95')

plot(freq_kenn$mids,freq_kenn$density, xlab="Frequency", ylab="No of Documents (Density)")
lines()
hist(counts$kennedy_norm, prob=T)
lines(density(counts$showed_norm))
xtable(df2)
the<-matrix(c(0.1,0.9), ncol = 2,nrow=3,byrow = TRUE)
f<-matrix(c(0.1,0.4,0.7,0.85,0.9,0.95),ncol = 2,nrow=3,byrow = TRUE)
alpha<-c(18,0.63)
x<-sort(rbtld(1000,alpha,the[1,]))
dens<-dbtld(x,alpha,the[1,])

plot(x,dens, type="l")


lapply(red[1:3], FUN=hist)
names<-c('said', 'showed','Kennedy')
list <-lapply(1:ncol(red),
              function(col) ggplot2::qplot(red[[col]],
                                           geom = c("histogram"), bins=10,
                                            xlab=names[[col]]))
ggally_barDiag(red,mapping = aes(x=said), binwidth = 2,rescale = T)+
  geom_density(bw = 0.1, alpha = 0.3)
hist(red$said, freq = F, main="Histogram of said", xlab="Counts")
lines(density(red$said))
hist(red$showed, freq = F, main="Histogram of showed", xlab="Counts")
lines(density(red$showed))
hist(red$Kennedy, freq = F, main="Histogram of Kennedy", xlab="Counts")
lines(density(red$Kennedy))

cowplot::plot_grid(plotlist = list,ncol = 3)
?plot_grid
?qplot
?geom_density


ecdf <- rep(1/nrow(red),nrow(red))
ecdf<- cumsum(ecdf)
plot(sort(red$said),ecdf,type = "l", lty="dashed",col="Blue", main = "H-Plot of Observed Probabilities", xlab = "Y", ylab="Cumulative Density")
par(new=TRUE)
plot(sort(red$showed),ecdf, col="Red", type="l", ann=F,axes = F)
par(new=TRUE)
plot(sort(red$Kennedy),ecdf, col="Green", type="l", ann=F,axes = F)
 #add a grid


legend(0.8, 0.2, c('Said','Showed','Kennedy'), col = c("Blue", "Red", "Green"),
       text.col = "black", lty = c("dashed", "solid", "dotted"),
       merge = TRUE, bg = "gray90")

the<-matrix(c(0.05,0.95), ncol = 2,nrow=3,byrow = TRUE)
f<-matrix(c(0.1,0.4,0.7,0.85,0.9,0.95),ncol = 2,nrow=3,byrow = TRUE)

btld_scales(the[1,],f[1,])
r<-cor(red)
r<-r[upper.tri(r)]
norm.cop <- normalCopula(param = r, dim = 3,dispstr = "un")
alpha<-matrix(c(4,24,28,6,36,2),ncol = 2,nrow = 3,byrow = T)
cc<-rmvbtld(n=352,theta=the,alpha=alpha,sigma = r,dim=3,dispstr = "un", normalize = T)


cc<-cc
normalize(cc, byrow=T)
norm<-1/rowSums(cc)
normed<-data.frame(norm*cc)

scatter3D(normed[,1], normed[,2], normed[,3], pch = 18,  theta = 40, phi = 40,
          main = "CMVBTL", xlab = "X1",
          ylab ="X2", zlab = "X3", ticktype="detailed",bty="g")
scatter3D(red[,1], red[,2], red[,3], pch = 18,  theta = 40, phi = 40,
          main = "CMVBTL", xlab = "X1",
          ylab ="X2", zlab = "X3", ticktype="detailed",bty="g")

colnames(normed)<-c('x1', 'x2','x3')
btld::rmvbtld()

ScatterHist(frame =normed,
            xvar = "x2",
            yvar= "x1",
            title = 'MVBTLD',
            contour = T,
            smoothmethod = 'none',
            point_alpha = 0.3,
            point_color ='Blue',
            hist_color = 'Green',
            contour_color = 'Red')

orig<-cc
colnames(orig)<-c('x1', 'x2','x3')


ScatterHist(frame =orig,
            xvar = "x2",
            yvar= "x1",
            title = 'MVBTLD',
            contour = T,
            smoothmethod = 'none',
            point_alpha = 0.3,
            point_color ='Blue',
            hist_color = 'Green',
            contour_color = 'Red')



par(mfrow=c(1,3))

xx <- cut(counts$said_norm, 10)
yy<- cut(counts$showed_norm, 10)
#
# ##  Calculate joint counts at cut levels:
z <- table(xx, yy)
?table
##  Plot as a 3D histogram:
hist3D(z=z,border="black", theta=120, phi=20,label=TRUE, ticktype='detailed',bty="g", xlab="P(X='said')", ylab="P(X='showed')", zlab="Frequencies", main="Observed Probabilities", colkey=F)

# fit a compound dirichlet mulnomial and plot the frequency distributions of the probabilities
red<-matrix(c(counts$said_norm,counts$showed_norm,counts$kennedy_norm), ncol = 3,nrow=352)
MGLMfit(data=red, dist = "DM")
n<-352
m <- 10
y <- rdirmn(n, m,c(0.354036,0.066088,0.024166))


xx <- cut(y[,1], 10)
yy<- cut(y[,2], 10)
#
# ##  Calculate joint counts at cut levels:
z <- table(xx, yy)
lims<-seq(1,10,1)
?table
##  Plot as a 3D histogram:
hist3D(z=z,border="black", theta=120, phi=20,label=TRUE, ticktype='detailed',bty="g",xlab="P(X='said')", ylab="P(X='showed')",zlab=" ", main="Estimated CDM Probabilities", colkey=F)


dim<-3
n<-352
sigma<-r
size<-10
df <- matrix(nrow=352,ncol=3)
dens <- matrix(nrow=352,ncol=1)
df<-mn_cmvbtld(352, n=10, alpha = alpha,theta = the,sigma=sigma,dim =3, base=1000, dispstr = "un")
for(i in seq_len(n)){
   p <- rmvbtld(1, alpha = alpha,theta = the,sigma=sigma,dim = dim, dispstr = "un")
  # p<-exp(p)/rowSums(exp(p))
   p<-1000^(p)/rowSums(1000^(p))
  df[i,] <- t(rmultinom(n=1,size = size,prob=p))
  # dens[i,]<-dmultinom(t(df[i,]), size=size, prob=as.numeric(p))
}

p <- rmvbtld(352, alpha = alpha,theta = the,sigma=sigma,dim = dim, dispstr = "un")
p<-1000^(p)/rowSums(1000^(p))
xx <- cut(p[,1], 10)
yy<- cut(p[,2], 10)
#
# ##  Calculate joint counts at cut levels:
z <- table(xx, yy)
lims<-seq(1,10,1)

##  Plot as a 3D histogram:
hist3D( z=z,border="black", theta=120, phi=20,label=TRUE, ticktype='detailed',bty="g",xlab="P(X='said')", ylab="P(X='showed')", zlab=" ", main="Estimated CMVTL-MN Probabilities", colkey=F)



df<-mn_cmvbtld(352, n=10, alpha = alpha,theta = the,sigma=r,dim =3, base=1000, dispstr = "un")
scatter3D(df[,1], df[,2], df[,3], pch = 18,  theta = 120, phi = 40,
          main = "CMVBTL-MN", xlab = "Said",
          ylab ="Showed", zlab = "Kennedy", ticktype="detailed",bty="g", colkey = FALSE)



# read in corpus counts
counts<-read.csv('/home/tristan/Documents/CodingProjects/GIt/brown_normalized_counts.csv')
xx <- cut(counts$said_norm, 10)
yy<- cut(counts$showed_norm, 10)
# ##  Calculate joint counts at cut levels:
z <- table(xx, yy)

##  Plot as a 3D histogram:
hist3D(z=z,border="black", theta=120, phi=20,label=TRUE, ticktype='detailed',bty="g", xlab="P(X='said')", ylab="P(X='showed')", zlab=" ", main="Observed Probabilities", colkey = F)

# predicted CMN-MVTL
the<-matrix(c(0.05,0.95), ncol = 2,nrow=3,byrow = TRUE)
f<-matrix(c(0.1,0.4,0.7,0.85,0.9,0.95),ncol = 2,nrow=3,byrow = TRUE) # determined from btld_scales function
btld_scales(the[1,],f[1,])
r<-cor(red)
r<-r[upper.tri(r)]

mn<-mn_cmvbtld(352, size=10, alpha=alpha, theta=the, sigma=r,dim=3,base=1000, dispstr = "un")
##  Plot as a 3D histogram:
par(mfrow=c(1,1))
scatter3D(mn[,1], mn[,2], mn[,3], pch = 18,  theta = 120, phi = 40,
          main = "Predicted CMVTL-MN", xlab = "Said",
          ylab ="Showed", zlab = "Kennedy", ticktype="detailed",bty="g", colkey=FALSE)


# correct mn-mvbtl
dim<-3
n<-1000
sigma<-r
size<-20
df <- matrix(nrow=1000,ncol=3)
# dens <- matrix(nrow=352,ncol=1)
for(i in seq_len(n)){
  p <- rmvbtld(1, alpha = alpha,theta = theta,sigma=0.5,dim = dim)
  p<-exp(p)/rowSums(exp(p))
  df[i,] <- t(rmultinom(n=1,size = size,prob=p))
  # dens[i,]<-dmultinom(t(df[1,]), size=size, prob=as.numeric(p))
}

par(mfrow = c(1, 1))

scatter3D(df[,1], df[,2], df[,3], pch = 18,  theta = 120, phi = 40,
          main = "CMVBTL-MN", xlab = "X1",
          ylab ="X2", zlab = "X3", ticktype="detailed",bty="g", colkey = FALSE)


par(mfrow = c(1, 3))

scatter3D(red[,1], red[,2], red[,3], pch = 18,  theta = 120, phi = 40,
          main = "Observed Counts", xlab = "Said",
          ylab ="Showed", zlab = "Kennedy", ticktype="detailed",bty="g", colkey=FALSE)

scatter3D(y[,1], y[,2], y[,3], pch = 18,  theta = 120, phi = 40,
          main = "CDM", xlab = "Said",
          ylab ="Showed", zlab = "Kennedy", ticktype="detailed",bty="g", colkey = FALSE)


scatter3D(df[,1], df[,2], df[,3], pch = 18,  theta = 120, phi = 40,
          main = "CMVBTL-MN", xlab = "Said",
          ylab ="Showed", zlab = "Kennedy", ticktype="detailed",bty="g", colkey = FALSE)



p <- rmvbtld(352, alpha = alpha,theta = the,sigma=r,dim = dim, dispstr = "un")
not_soft<-data.frame(x1=p[,1],x2=p[,2])

ScatterHist(frame =not_soft,
            xvar = "x1",
            yvar= "x2",
            title = 'MVBTLD',
            contour = T,
            smoothmethod = 'none',
            point_alpha = 0.3,
            point_color ='Blue',
            hist_color = 'Green',
            contour_color = 'Red')


p<-1000^(p)/rowSums(1000^(p))
soft<-data.frame(x1=p[,1], x2=p[,2])

ScatterHist(frame =soft,
            xvar = "x1",
            yvar= "x2",
            title = 'MVBTLD (Normalized)',
            contour = T,
            smoothmethod = 'none',
            point_alpha = 0.3,
            point_color ='Blue',
            hist_color = 'Green',
            contour_color = 'Red')


help(faithful)
df<-datasets::faithful
plot(df$eruptions,df$waiting)
hist(data$eruptions, breaks = 30)
hist(df$waiting)




df2<-MASS::geyser

df2<- predict(pre, df2)
plot(df2$duration,df2$waiting)
hist(df2$duration, breaks = 30)
hist(df2$waiting, breaks = 30, freq = F)
btld::get_params_btld(df2$waiting)


H=bw.SJ(c(df2$waiting,df2$duration));
grid=seq(min(df2$duration)-3*H,max(df2$duration)+3*H,by=H/10)
gfungrid=c();

fn <- function(grid){
for(i in 1:length(grid))
  gfungrid[i]=( dnorm((df2[,1]-grid[i])/H)/H+dnorm((df2[,2]-grid[i])/H)/H )/n
###############################################
}

init <- c(0,0)
opt <- ucminf(init, gfungrid)$par

H=bw.SJ(c(Residual.1,Residual.2));
###############################################  
grid=seq(min(R.matrix)-3*H,max(R.matrix)+3*H,by=H/10)
gfungrid=c();
for(i in 1:length(grid))
  gfungrid[i]=( t(P[,1])%*%dnorm((R.matrix[,1]-grid[i])/H)/H+t(P[,2])%*%dnorm((R.matrix[,2]-grid[i])/H)/H )/n




fn <- function(x1,x2){
  d1 <- density(x1, n=298)
  d2 <- density(x2, n=298)
  return(c(d1,d2))
}

opt <- fn(df2$waiting, df2$duration)
max(opt$x)
max(opt$y)
wait <- opt$x
dur <- opt$y
first<-df2[which.max(df$eruptions),1]
second<-df2[which.max(dur),2]

tt <- c(1,2,3,2,1, 1, 2, 1)
idx<-which(diff(sign(diff(df$waiting)))==-2)+1

df<-df$waiting[idx]
idx<-which(diff(sign(diff(df)))==-2)+1
df[idx]


wait
which.max(wait)
which.max(dur)
nrow(df2$waiting)
length(df2$waiting)
which.peaks(df2$waiting)
which.peaks(df2$duration)

library(ggplot2)
bimod_lymphoma<-read.csv('explorer_download.csv')
Overall.Cond <-seq(0,85,5)
Freq <- as.numeric(bimod_lymphoma$Rate.per.100.000)
density<-Freq[!is.na(Freq)]

modeidxs<-which(diff(sign(diff(density)))==-2)+1
density[modeidxs]



dat<-data.frame(x=Overall.Cond,dens=density)
xs<-data.frame(x=Overall.Cond)
pre<-preProcess(dat,method="range")
df2<- predict(pre, dat)
library(pracma)

des<-dbtld(df2$x,c(5,5),c(0.25,0.9))


fx<-empcdf(dat$dens)

dat<-data.frame(df2,dens=dat$dens, pdf=fx)
dat
new<-dat[modeidxs,]
theta<-new$x
theta
new
alpha<-btld_scales(new$x,new$pdf.FnX)
btld

model<-dbtld(dat$x,alpha,theta)
dat<-data.frame(dat,model=model)
str(dat)
rm(g)
g<-ggplot(dat,aes(x = x,weight = dens.1)) +
  geom_histogram( bins = 30 )
g<-g+geom_histogram(aes(y=..density..))
g
g <- ggplot(df4, aes(x))
g <- g + geom_line(aes(y=model), colour="red")
g <- g + geom_line(aes(y=k30), colour="green")
g<-g+geom_line(aes(y=k40),colour="blue")

colnames(df)<-c('x', 'kNN_k=20', 'kNN_k=30', 'GMM_full','t','c','gmm_comp1','gmm_comp2')
g
help(geom_histogram)

lower_x<-0
upper_x<-1
width<-abs(upper_x-lower_x) # distance on x plane
x1<-runif(1000, lower_x, upper_x) # uniform random variables on x plane
z<-dbtld(x1,  alpha = , theta =c(0.3, 0.7)) # function values giving values on z plane
height <- mean(z) # distance on z plane
area_mcmc <- width*height
area_mcmc

library('mclust')
data(GvHD)
plot(GvHD.control)
plot(GvHD.pos)

CD4<-GvHD.control$CD4
obj<-hist(CD4)
density<-obj$density

modeidxs<-which(diff(sign(diff(density)))==-2)+1
mids<-obj$mids
mids[modeidxs]

unit.normalize<-function(x){
  xs.df<-data.frame(x)
  pre<-preProcess(xs.df,method="range")
  unit.df<- predict(pre, xs.df)
  return(unit.df)
}

find.thetas<-function(x){
  freq<- hist(x, include.lowest=TRUE, plot=FALSE,breaks=20)
  density<-freq$density
  modeidxs<-which(diff(sign(diff(density)))==-2)+1
  dens<-freq$dens
  mids<-freq$mids
  ms<-cbind(dens,mids)
  datr<-ms[modeidxs,]
  return(datr)
}
find.f_thetas<-function(x){
  unit<-unit.normalize(x)
  modes<-find.thetas(unit)
  ecdf<-empcdf(unit)
  rounded<-round(ecdf[,1],digits=2)
  mods<-modes[,1]
  f_thet<-matrix(nrow=2,ncol=2)
  for(i in seq_len(2)){
    idxs<-which(rounded==as.character(mods[i]))
    ids<-ecdf[idxs,]
    f_thet[i,]<-colMeans(ids)
  }
  return(f_thet)
}

find.alphas<-function(x){
  unit<-unit.normalize(x)
  f_thetas<-find.f_thetas(x)
  modes<-find.thetas(unit)
  alphas<-btld_scales(modes,f_thetas)
  return(alphas)
}
btld.params<-function(x){
  unit<-unit.normalize(x)
  f_thetas<-find.f_thetas(x)
  modes<-find.thetas(unit)
  alphas<-btld_scales(modes,f_thetas)
  parm.df<-data.frame(modes=modes,alphas=alphas)
  return(parm.df)
}
idxs<-which(rounded==as.character(mods[2]))
ids<-ecdf[idxs,]
f_thet2<-colMeans(ids)
f_thet<-rbind(f_thet1,f_thet2)

unit.CD4<-unit.normalize(GvHD.control$CD4)
alpha<-c(1.5,2.2)
modes<-c(0.15,0.45)
theta<-f_thet[,1]
d<-btld.params(unit.CD4)
x<-dbtld(sort(unit.CD4),alpha =alpha,theta = modes)
plot(sort(unit.CD4),x, type='l')
soted<-sort(unit.CD4)
hist(unit.CD4,freq=F,col="red", main="GvHD: CD4", xlab="x", breaks=30)
lines(x=sort(unit.CD4),y=x, col="green")
cdf<-pbtld(soted,alpha =alpha,theta = modes)
plot(soted,cdf, type="l", col="green",lty="dashed", main="KS-testing", xlab = 'x', ylab="Cumulative Probabilities")
ecdf<-empcdf(soted)
lines(soted,sort(ecdf[,2]),type="l",col="blue",lty="dotted")

legend(0.7, 0.2, c('BTLD(0.15,0.45,2.5,1.5)','ECDF'), col = c("Blue", "Green"),
       text.col = "black", lty = c("dashed", "dotted"),
       merge = TRUE, bg = "gray90")
# max(cdf-ecdf)
# kolsmir_test<-function(rvs){

kolsmir_test<-function(prop_cdf){
  library(KSgeneral)
  n<-length(cdf)
  i<-t(1:n)
  Dminus <- prop_cdf-(i-1)/n
  Dplus <- i/n - prop_cdf
  D <- max(c(Dminus,Dplus))
  print("KStest Stat:"); print(D)
  print(D*sqrt(n)>=1.63)
  pval<-cont_ks_c_cdf(D,n)
  print(paste0("p-value: ", pval))
}
kolsmir_test(cdf)

rm(dbtld)
latex2exp_examples()

hist(GvHD.control$CD3)

modes<-c(0.10000, 0.850000)
thets<-data.frame(find.thetas(soted))
thets[order(thets, decreasing=TRUE),]
alpha<-c(2,1.5)
modes<-c(0.075,0.875)
library(mixtools)
data(NOdata)
 parms<-btld.params(NOdata$NO, normalize=TRUE)
soted<-sort(as.matrix(unit.normalize(NOdata$NO)))
hist(soted)
x<-dbtld(sort(as.matrix(soted)),alpha =alpha,theta = modes)
plot(sort(as.matrix(soted)),x, type='l')
# soted<-sort(unit.CD4$x)
hist(soted,freq=F,col="red", main="NO2data: NO", xlab="x", breaks = 30, ylim = c(0,3))
lines(x=sort(soted),y=x, col="green")
cdf<-pbtld(sort(soted),alpha =alpha,theta = modes)
plot(sort(soted),cdf, type="l", col="green")
ecdf<-empcdf(soted)
lines(sort(soted),sort(ecdf[,2]),type="l",col="blue")
kolsmir_test(cdf)

hist(NOdata$Equivalence, breaks = 30, freq=F)

hist(GvHD.control$CD4)
hist(GvHD.control$CD8)
hist(GvHD.control$CD8b)
hist(GvHD.control$CD3, freq=F)


unit.CD4<-unit.normalize(NOdata$Equivalence)
alpha<-c(1.477273,1.25)
modes<-c(0.35,0.95)
theta<-f_thet[,1]
d<-btld.params(unit.CD4)
x<-dbtld(sort(unit.CD4),alpha =alpha,theta = modes)
plot(sort(unit.CD4),x, type='l')
soted<-sort(unit.CD4)
hist(unit.CD4,freq=F,col="red", main="GvHD: CD4", xlab="x")
lines(density(unit.CD4))
DD<-density(unit.CD4)
dens<-DD$y
length(dens)
modeidxs
modeidxs<-which(diff(sign(diff(dens)))==-2)+1
dens[modeidxs]
freq<- hist(unit.CD4, include.lowest=TRUE, plot=FALSE)
density<-freq$density
modeidxs<-which(diff(sign(diff(density)))==-2)+1
dens<-freq$dens
lines(x=sort(unit.CD4),y=x, col="green")
cdf<-pbtld(soted,alpha =alpha,theta = modes)
plot(soted,cdf, type="l", col="green",lty="dashed", main="KS-testing", xlab = 'x', ylab="Cumulative Probabilities")
ecdf<-empcdf(soted)
lines(soted,sort(ecdf[,2]),type="l",col="blue",lty="dotted")
kolsmir_test(cdf)
library('')
gene_data<-read_tsv('/home/tristan/Documents/DataSpellProjects/GIt/Bimodality_Genes/data/KIRC_FPKM_TP.tsv')
idx<-which(gene_data$symbol=="ZDHHC15")
uty<-gene_data[idx,]
logged<-log(uty$FPKM)
logged<-logged[!logged< -10]
hist(logged,freq=F,breaks=30)
unit.uty<-unit.normalize(logged)

hist(s,freq=F,breaks=30)
# softmax(uty$FPKM)
soft<-(10^unit.uty)/(10^(sum(unit.uty)))
btld.params(unit.uty)

samples<-sample(unit.uty, 100)
hist(samples,freq=F,breaks=30)
gene<-"ZDHHC15"

gene_data<-read_tsv('data/KIRC_FPKM_TP.tsv')
idx<-which(gene_data$symbol==gene)
gene<-gene_data[idx,]
hist(gene$FPKM,freq=F,breaks=30)
unit.gene<-unit.normalize(gene$FPKM)
hist(unit.gene,freq=F,breaks=30)
alpha<-c(4,1.5)
theta<-c(0.1,0.5)
sort.unit.gene<-sort(unit.gene)
x<-dbtld(sort.unit.gene,alpha =alpha,theta = theta)
hist(sort.unit.gene, freq = F,breaks = 30, col="blue", ylim = c(0,8))
lines(sort.unit.gene,x, type='l', col="red")
# plot cdf
cdf<-pbtld(sort.unit.gene,alpha =alpha,theta = theta)
plot(sort.unit.gene,cdf, type="l", col="green",lty="dashed", main="KS-testing", xlab = 'x', ylab="Cumulative Probabilities")
ecdf<-empcdf(sort.unit.gene)
lines(sort.unit.gene,sort(ecdf[,2]),type="l",col="blue",lty="dotted")
kolsmir_test(cdf)

# btld_scales(c(0.1,0.6),c(0.4,0.9))
alpha<-btld_scales(c(0.1,0.5),c(0.4,0.9))
ks.test(sort.unit.gene)
cvm.test(sort.unit.gene,"pbtld",alpha=alpha,theta=theta,estimated=TRUE)
?ad.test

# suppose we log it
logged<-log(unit.gene)
logged<-logged[!logged< -10]
unit.uty<-unit.normalize(logged)
ecdf<-empcdf(sort(unit.uty))
plot(sort(unit.uty),sort(ecdf[,2]),type="l",col="blue",lty="dotted")
hist(unit.uty,freq=F,breaks=30)
sort.unit.gene<-sort(unit.uty)
ecdf<-ecdf[round(ecdf[,1], digits = 2)==0.45,]
ft1<-colMeans(ecdf)

ecdf<-empcdf(sort(unit.uty))
ecdf<-ecdf[round(ecdf[,1], digits = 2)==0.85,]
ft2<-colMeans(ecdf)
alpha<-btldScales(c(ft1[1], ft2[1]),c(ft1[2],ft2[2]))
theta<-c(ft1[1], ft2[1])
p<-(cbind(alpha, theta))
# sort.unit.gene<-sort(unit.gene)
colnames(p)<-c('$\alpha$', '$\theta$')
xtable(p)
x<-dbtld(sort.unit.gene,alpha =alpha,theta = theta)
hist(sort.unit.gene, freq = F,breaks = 30, col="blue", ylim = c(0,8))
lines(sort.unit.gene,x, type='l', col="red")
# plot cdf
cdf<-pbtld(sort.unit.gene,alpha =alpha,theta = theta)
plot(sort.unit.gene,cdf, type="l", col="green",lty="dashed", main="KS-testing", xlab = 'x', ylab="Cumulative Probabilities")
ecdf<-empcdf(sort.unit.gene)
lines(sort.unit.gene,sort(ecdf[,2]),type="l",col="blue",lty="dotted")
kolsmir_test(cdf)

# softmax(uty$FPKM)
soft<-(10^unit.uty)/(10^(sum(unit.uty)))
btld.params(unit.uty)

samples<-sample(unit.uty, 100)
hist(samples,freq=F,breaks=30)


ref.labels <- c(rep("A", 45), rep("B" , 10), rep("C", 15), rep("D", 25), rep("E", 5))
predictions <- c(rep("A", 35), rep("E", 5), rep("D", 5),
                 rep("B", 9), rep("D", 1),
                 rep("C", 7), rep("B", 5), rep("C", 3),
                 rep("D", 23), rep("C", 2),
                 rep("E", 1), rep("A", 2), rep("B", 2))
df <- data.frame("Prediction" = predictions, "Reference" = ref.labels, stringsAsFactors=TRUE)