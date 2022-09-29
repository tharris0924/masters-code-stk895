library(ucminf)
library(caret)
library(btld)
library(copula)
library(WVPlots)
library(ggplot2)
library(mltest)
library(plot3D)
library(latex2exp)
library(xtable)
library(LICORS)
library('mclust')
library('readr')
library('MGLM')
# read in corpus counts
counts<-read.csv('/home/tristan/DataSpellProjects/GIt/brown_normalized_counts.csv')
red<-data.frame(said=counts$said_norm/10, showed=counts$showed_norm/10,Kennedy = counts$kennedy_norm/10)
xx <- cut(red[,1], 10)
yy<- cut(red[,2], 10)
# ##  Calculate joint counts at cut levels:
z <- table(xx, yy)
plot(counts)
##  Plot as a 3D histogram:
hist3D(z=z,border="black", theta=120, phi=20,label=TRUE, ticktype='detailed',bty="g", xlab="P(X='said')", ylab="P(X='showed')", zlab=" ", main="Observed Probabilities", colkey = F)

# fit CDM: note that the maximum likelihood fitting procedure orders parameter estimates alphabetically not according to
# the original dimensions of the data of "red". Also the data frame needs to be in counts not frequency Probabilities.
red<-data.frame(said=counts$said_norm,showed=counts$showed_norm,  Kennedy = counts$kennedy_norm)
fit<-MGLMfit(data=red, dist = "DM")
AIC(fit)
BIC(fit)
n<-352
m <- 10
y <- rdirmn(n, m,c(0.354036,0.066088, 0.024166))
loglik<-sum(ddirmn(y,c(0.354036,0.066088, 0.024166)))
mltest::ml_test(y,as.matrix(red))
xx <- cut(y[,1], 10)
yy<- cut(y[,2], 10)
#
# ##  Calculate joint counts at cut levels:
z <- table(xx, yy)
bic_bt(loglik,4,352)
##  Plot as a 3D histogram:
hist3D(z=z,border="black", theta=120, phi=20,label=TRUE, ticktype='detailed',bty="g",xlab="P(X='said')", ylab="P(X='showed')",zlab=" ", main="Predicted CDM Probabilities", colkey=F)


# predicted CMVTL-MN
the<-matrix(c(0.05,0.95), ncol = 2,nrow=3,byrow = TRUE)
f<-matrix(c(0.1,0.4,0.7,0.85,0.9,0.95),ncol = 2,nrow=3,byrow = TRUE) # determined from btld_scales function
btld_scales(the[2,],f[2,])
r<-cor(red)
r<-r[upper.tri(r)]

alpha<-matrix(c(4,24,28,6,36,2),ncol = 2,nrow = 3,byrow = T)
dmn(mn,p)
mn<-mn_cmvbtld(352, size=10, alpha=alpha, theta=the, sigma=r,dim=3,base=1000, dispstr = "un")
p<-cmpmvbtl(352,alpha,the,r,3,1000,dispstr = "un")
loglik<-sum(dmn(mn,p))

bic_bt(loglik,4,n)
aic_bt(loglik,4)
xx <- cut(mn[,1], 10)
yy<- cut(mn[,2], 10)
#
# ##  Calculate joint counts at cut levels:
z <- table(xx, yy)

mn_cmvbtld22<-function(n,size,alpha,theta,sigma,dim,base,...){
  df <- matrix(nrow=n,ncol=dim)
  loglik<-matrix(nrow=n,ncol=dim)
  for(i in seq_len(n)){
    p <- rmvbtld(1, alpha = alpha,theta = theta,sigma=sigma,dim = dim,...)
    p<-softmax(p,exp=F,base=base)
    mn<-rmultinom(n=1,size = size,prob=p)
    df[i,] <- t(mn)
    loglik[i,]<-p
  }
  loglik<-sum(dmn(df,loglik))
  ls<-list(mn=df,loglik=loglik)
  return(ls)
}

##  Plot as a 3D histogram:
hist3D( z=z,border="black", theta=120, phi=20,label=TRUE, ticktype='detailed',bty="g",xlab="P(X='said')", ylab="P(X='showed')", zlab=" ", main="Predicted CMVTL-MN Probabilities", colkey=F)
base<-1000
p <- t(rmvbtld(1, alpha = alpha,theta = the,sigma=r,dim = 3,dispstr = "un"))
p<-softmax(p,exp=F,base=base)
df <- t(rmultinom(n=1,size = 10,prob=p))
par(mfrow = c(1, 1))

cmpmvbtl<-function(n,alpha,theta,sigma,dim,base,...){
  comp <- matrix(nrow=n,ncol=dim)
  for(i in seq_len(n)){
  p <- rmvbtld(1, alpha = alpha,theta = theta,sigma=sigma,dim = dim,...)
  p<-softmax(p,exp=F,base=base)
  comp[i,] <- p
  }
  return(comp)
}



fit<-MGLMfit(data=red, dist = "DM")
AIC(fit)
BIC(fit)
scatter3D(red[,1], red[,2], red[,3], pch = 18,  theta = 120, phi = 40,
          main = "Observed Counts", xlab = "Said",
          ylab ="Showed", zlab = "Kennedy", ticktype="detailed",bty="g", colkey=FALSE)

scatter3D(y[,1], y[,2], y[,3], pch = 18,  theta = 120, phi = 40,
          main = "Predicted CDM", xlab = "Said",
          ylab ="Showed", zlab = "Kennedy", ticktype="detailed",bty="g", colkey = FALSE)

mn<-mn_cmvbtld(352, size=10, alpha=alpha, theta=the, sigma=r,dim=3,base=1000, dispstr = "un")
mn<-mn_cmvbtld22(352, size=10, alpha=alpha, theta=the, sigma=r,dim=3,base=1000, dispstr = "un")
scatter3D(mn[,1], mn[,2], mn[,3], pch = 18,  theta = 120, phi = 40,
          main = "Predicted CMVTL-MN", xlab = "Said",
          ylab ="Showed", zlab = "Kennedy", ticktype="detailed",bty="g", colkey = FALSE)
aic_bt<-function(loglik,k){
  return((-2*loglik)+(2*k))
}
bic_bt<-function(loglik,k,n){
  return((-2*loglik)+(k*log(n)))
}
aic_bt(mn$loglik, 6)
bic_bt(mn$loglik,6, 352)
# comparison for CDM
actuals<-as.matrix(red)
cdn<-mltest::ml_test(y,actuals)
ydf<-as.data.frame(y)
cols<-c('said','showed', 'Kennedy')
colnames(ydf)<-cols
ydf$type<-"DCM"


actuals_df<-as.data.frame(actuals)
dplyr::setdiff(actuals_df,ydf)
# comparison for CMVTL-MN
met_mtl<-mltest::ml_test(mn,actuals)

rnew<-red
rnew$type <- "Observed"
mn<-data.frame(mn$mn)
colnames(mn)<-cols
mn$type <- "CMVTL-MN"

df<-rbind(rbind(mn,rnew),ydf)
ggally_points(data = df, mapping=aes(x=showed,y=said,pch=type))
pal<-RColorBrewer::brewer.pal(3,"Dark2")
# said and showed: pred v observed
plot(mn$showed, mn$said, pch=1,col =pal[1], main="Showed and Said: Predicted vs Observed",xlab = "Counts of Showed",ylab="Counts of Said")
par(new=TRUE)
plot(rnew$showed,rnew$said, pch =20,col = pal[2], axes=F, xlab=" ", ylab = " ")
par(new=T)
plot(ydf$showed, ydf$said, pch=3, col="Black", axes=F, xlab= " ", ylab= " ")
legend(8,8,c('CMVTL-MN', "Observed","DCM"),pch=c(1,20,3),col=c(pal[1], pal[2],"Black"),box.col = "Gray")

plot(mn$showed, mn$Kennedy, pch=1,col =pal[1], main="Showed and Kennedy: Predicted vs Observed",xlab = "Counts of Showed",ylab="Counts of Kennedy")
par(new=TRUE)
plot(rnew$showed,rnew$Kennedy, pch =20,col = pal[2], axes=F, xlab=" ", ylab = " ")
par(new=T)
plot(ydf$showed, ydf$Kennedy, pch=3, col="Black", axes=F, xlab= " ", ylab= " ")
legend(8,8,c('CMVTL-MN', "Observed","DCM"),pch=c(1,20,3),col=c(pal[1], pal[2],"Black"),box.col = "Gray")

plot(mn$said, mn$Kennedy, pch=1,col =pal[1], main="Said and Kennedy: Predicted vs Observed",xlab = "Counts of Said",ylab="Counts of Kennedy")
par(new=TRUE)
plot(rnew$said,rnew$Kennedy, pch =20,col = pal[2], axes=F, xlab=" ", ylab = " ")
par(new=T)
plot(ydf$said, ydf$Kennedy, pch=3, col="Black", axes=F, xlab= " ", ylab= " ")
legend(8,8,c('CMVTL-MN', "Observed","DCM"),pch=c(1,20,3),col=c(pal[1], pal[2],"Black"),box.col = "Gray")
