
library('btld')
alpha<- matrix(c(1,5,3,3,5,1),nrow = 3,ncol = 2, byrow = TRUE)
theta<-matrix(c(0.3,0.7,0.3,0.7,0.3,0.7), nrow = 3,ncol = 2, byrow = TRUE)

X<-rmvbtld(1000, alpha=alpha, theta=theta, sigma=0.5,dim=3)
df<-as.data.frame(X)
colnames(df)<-c(quote(X[1]),quote(X[2]),quote(X[3]))

library('GGally')
ggpairs(df,upper=list(continuous="density"))
ggpairs(df,diag=list(continuous="barDiag"),upper=list(continuous="density"),title="Multivariate Bimodal Triangular-Linked Distribution", columnLabels =  c(latex2exp::TeX("$X_1$"),latex2exp::TeX("$X_2$"),latex2exp::TeX("$X[3]$")))
ggpairs(df,upper=list(continuous="density"))

p_(pm)

df<-as.data.frame(df)
colnames(df)<-c('x1','x2','x3')
scatter3D(df[,3], df[,2], df[,1], pch = 18,  theta = 20, phi = 40,
          main = "Simulated Example: MVTL Scatterplot Identical Marginals ", xlab = "X1",
          ylab ="X2", zlab = "X3", ticktype="detailed",bty="g")
df$type<-"MVTL"
x1<-as.matrix(rbtld(1000,alpha[1,],theta[1,]))
log_lik<-log(sum(dbtld(x1,alpha[1,],theta[1,])))
it<-mixtools::normalmixEM(x1, k=2)
cdf<-mixtools::compCDF(x1,it$posterior)
x3<-rbtld(1000,alpha[3,],theta[3,])
x2<-rbtld(1000,alpha[2,],theta[2,])
df2<-data.frame(x1=x1,x2=x2,x3=x3,type="Independant")
df_new<-rbind(df,df2)
ggally_density(df_new,ggplot2::aes(x=x1,y=x2, color=type))+xlab(quote(X[1]))+ylab(quote(X[2]))+ggplot2::theme(legend.title = element_text(size = rel(1.25)),legend.text = element_text(size = rel(1.25)), axis.title.y = element_text(size = rel(1.25)),axis.title.x = element_text(size = rel(1.25)))
ggally_density(df_new,ggplot2::aes(x=x3,y=x2, color= type ))+xlab(quote(X[3])) +ylab(quote(X[2]))+ggplot2::theme(legend.title = element_text(size = rel(1.25)),legend.text = element_text(size = rel(1.25)), axis.title.y = element_text(size = rel(1.25)),axis.title.x = element_text(size = rel(1.25)))
ggally_density(df_new,ggplot2::aes(x=x3,y=x1,color=type))+xlab(quote(X[3])) +ylab(quote(X[1]))+ggplot2::theme(legend.title = element_text(size = rel(1.25)),legend.text = element_text(size = rel(1.25)), axis.title.y = element_text(size = rel(1.25)),axis.title.x = element_text(size = rel(1.25)))

pm<-ggmatrix(plotList,3,1,data=df_new)

par(mfrow=c(1,1))
p_(pm)

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
# AIC tests
btld_bic<- -2*log(prod((dbtld(x1,alpha[1,],theta[1,]))))+4*log(1000)
btld_aic<- -2*log(prod((dbtld(x1,alpha[1,],theta[1,]))))+2*4
mix_aic<- -2*it$loglik+2*2*3
mix_bic <- -2*it$loglik+6*log(1000)
# mixture models
library('MixAll')
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

ems<-em_func(df, plot = T, nbCluster=2)
xtable::xtable(ems$MU, digits = 4)
# biostatistics (gene expression applications)
library('readr')
gene_data<-read_tsv("/home/tristan/DataSpellProjects/Bimodality_Genes/data/KIRC_FPKM_TP.tsv")

# ZDHHC15 and XRRA1 yield some interesting results
idx <- which(gene_data$symbol == "ZDHHC15")
gene<-gene_data[idx,]

logged <- log(gene$FPKM)
logged<-logged[!logged< -10]
unit.uty<-unit.normalize(logged)
ecdf<-empcdf(sort(unit.uty))
plot(sort(unit.uty),sort(ecdf[,2]),type="l",col="blue",lty="dotted")
hist(unit.uty,freq=F,breaks=30)
sort.unit.gene<-sort(unit.uty)
ecdf<-ecdf[round(ecdf[,1], digits = 2)==0.45,]
ft1<-colMeans(ecdf)
ecdf<-empcdf(sort(unit.uty))
ecdf<-ecdf[round(ecdf[,1], digits = 2)==0.80,]
ft2<-colMeans(ecdf)
alpha<-btldScales(c(ft1[1], ft2[1]),c(ft1[2],ft2[2]))
theta<-c(ft1[1], ft2[1])
# sort.unit.gene<-sort(unit.gene)
x<-dbtld(sort.unit.gene,alpha =alpha,theta = theta)
hist(unit.uty, freq = F,breaks = 30, col="blue", ylim = c(0,8))
lines(sort.unit.gene,x, type='l', col="red")
# plot cdf
cdf<-pbtld(sort.unit.gene,alpha =alpha,theta = theta)
plot(sort.unit.gene,cdf, type="l", col="green",lty="dashed", main="Goodness of Fit Testing", xlab = 'x', ylab="Cumulative Probabilities")
ecdf<-empcdf(sort.unit.gene)
lines(sort.unit.gene,sort(ecdf[,2]),type="l",col="blue",lty="dotted")
pnormmix <- function(x,mixture) {
  lambda <- mixture$lambda
  k <- length(lambda)
  pnorm.from.mix <- function(x,component) {
    lambda[component]*pnorm(x,mean=mixture$mu[component],
                            sd=mixture$sigma[component])
  }
  pnorms <- sapply(1:k,pnorm.from.mix,x=x)
  return(rowSums(pnorms))
}
gmm<-mixtools::normalmixEM(sort.unit.gene)
gmmCDF<-pnormmix(sort.unit.gene,gmm)
lines(sort.unit.gene,sort(gmmCDF),type="l",col="red",lty=4)
# softmax(uty$FPKM)
soft<-(10^unit.uty)/(10^(sum(unit.uty)))
btld.params(unit.uty)

samples<-sample(unit.uty, 100)
hist(samples,freq=F,breaks=30)

goftest::cvm.test(sort.unit.gene,"pnormmix",mixture=gmm,estimated=TRUE)
ks.test(sort.unit.gene,"pnormmix",mixture=gmm)
gmm<-mixtools::normalmixEM(sort.unit.gene)
aic_bt<-function(loglik,k){
  return((-2*loglik)+(2*k))
}
bic_bt<-function(loglik,k,n){
  return((-2*loglik)+(k*log(n)))
}
bic(gmm$loglik,3,510)
aic(gmm$loglik,3)
btlLoglik<-function(mat,alpha,theta,n){
  # laplace correction to avoid singularities
  density<-dbtld(mat,alpha,theta)+0.00001
  btld_bic<- -2*log(prod(density))+2*log(n)
  btld_aic<- -2*log(prod(density))+2*2
  return(list(bic=btld_bic,aic=btld_aic))

}
btlLoglik(sort.unit.gene,alpha,theta,510)

pnormmix <- function(x,mixture) {
  lambda <- mixture$lambda
  k <- length(lambda)
  pnorm.from.mix <- function(x,component) {
    lambda[component]*pnorm(x,mean=mixture$mu[component],
                            sd=mixture$sigma[component])
  }
  pnorms <- sapply(1:k,pnorm.from.mix,x=x)
  return(rowSums(pnorms))
}

gmmCDF<-pnormmix(sort.unit.gene,gmm)

# Theoretical CDF
# Empirical CDF
distinct.snoq <- sort(unique(sort.unit.gene))
tcdfs <- pnormmix(distinct.snoq,mixture=snoq.k2)
ecdfs <- ecdf(sort.unit.gene)(distinct.snoq)
plot(gmmCDF,ecdfs ,xlab="Theoretical CDF",ylab="Empirical CDF",xlim=c(0,1),
      col="Red", ylim=c(0,1),main="PP-Plot")
abline(0,1)
par(new=TRUE)
plot(cdf,ecdfs,axes = F, xlab=" ", ylab=" ", col="green", pch = "*")
legend(0.9,0.3,c("BTL","GMM"),pch=c("*","o"), col=c("Green","Red"))
plot(sort.unit.gene,cdf, type="l", col="green",lty=2, main="Goodness of Fit testing", xlab = 'x', ylab="Cumulative Probabilities")
ecdf<-empcdf(sort.unit.gene)
lines(sort.unit.gene,sort(ecdf[,2]),type="l",col="blue",lty=3)
gmmCDF<-pnormmix(sort.unit.gene,gmm)
lines(sort.unit.gene,sort(gmmCDF),type="l",col="red",lty=4)
legend(0.8,0.2,c("BTL","ECDF","GMM"),lty=c(2,3,4),col=c("Green", "Blue", "Red"))
