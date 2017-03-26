#####MME 2017

library(MASS)
library(mvtnorm)

car.f = read.csv2(file="Motorins_F_MOD.csv")
X     = as.matrix(read.csv2(file="design_matrix.csv"))
Z     = as.matrix(read.csv2(file="dummy_matrix.csv"))

cs    = c("A","B","C","D","E","F","G","H","I")
eg    = expand.grid(cs[1:5],cs[1:7],cs[1:7],cs[1:9])
er    = paste(eg[,1],eg[,2],eg[,3],eg[,4],sep="")
sr    = unique(with(car.f,paste(Kilometres,Zone,Bonus,Make,sep="")))
mi    = setdiff(er,sr)
mm    = mapply(function(i)substr(mi,i,i), 1:4)   #brakuj¹ce kombinacje liter

fX    = X[NULL,]
for(i in 1:nrow(mm))
  fX  = rbind(fX,as.numeric(c(1,mm[i,1]==cs[2:5],mm[i,3]==cs[2:7],mm[i,4]==cs[2:9])))

mz    = paste(mm[,2],mm[,3],sep="")
cz    = as.vector(t(outer(cs[1:7],cs[1:7],paste,sep="")))

fZ    = Z[NULL,]
for(i in 1:nrow(mm))
  fZ  = rbind(fZ,cz==mz[i])

fZB=paste(mm[,2],mm[,3],sep="")

glmm  = glmmPQL(fixed=Claims~Kilometres+Bonus+Make+offset(log(Insured)),
              random=~1|Zone_Bonus, family=poisson(link="log"), data=car.f) 


Beta  = as.matrix(fixef(glmm))
v     = as.matrix(ranef(glmm))

sigma2_v      = (0.218182)^2              
Y.fv          = exp(fX%*%Beta+fZ%*%v)     

ile=1000

BUF.g  = array(NA,dim=c(ile,nrow(fZ)))
BUF.gg = array(NA,dim=c(ile,nrow(fZ)))

set.seed(456)

for(i in 1:ile)
{
 v.g    = t(rmvnorm(n=1, mean=rep(0,nrow(v)), sigma=sigma2_v*diag(nrow(v))))
 Y.g    = rpois(n=nrow(car.f), lambda=exp(X%*%Beta+Z%*%v.g))
 Y.fv.g = rpois(n=nrow(fZ),    lambda=exp(fX%*%Beta+fZ%*%v.g))

 car.g = car.f
 car.g$Claims = Y.g 

 glmm.g = glmmPQL(fixed=Claims~Kilometres+Bonus+Make,
                random=~1|Zone_Bonus, family=poisson(link="log"), data=car.g) 

 Beta.gg  = as.matrix(fixef(glmm.g))
 v.gg     = as.matrix(ranef(glmm.g))

 Y.fv.gg  = exp(fX%*%Beta.gg+fZ%*%v.gg)

 BUF.g[i,]  = Y.fv.g   
 BUF.gg[i,] = Y.fv.gg   

 plot(Y.fv.g, Y.fv.gg, main=i)
}

ERR = abs(BUF.gg - BUF.g)
Q   = apply(ERR,2,quantile,probs=c(0.5,0.75,0.90,0.95))
MSE = sqrt(colMeans(ERR^2))
QM  = rbind(MSE,Q)

### FIGURE 2 ###
par(mar=c(3,4,3,2))
R=range(QM)
for(i in 1:ncol(QM))
{
  plot(QM[,i],xaxt="n",type="l",ylim=R,ylab="Measured accuracy",xlab="")
  par(new=T)
}
labele = c("RMSE",expression(q[0.5]),expression(q[0.75])
          ,expression(q[0.90]),expression(q[0.95]))

mtext(labele,at=1:5,side=1)

### FIGURE 1 ###
library(beanplot)

RR=c(0,3)
MM=c(-0.3,0.3)

par(mar=c(3,5,3,2))
par(mfrow=c(1,1))
beanplot(ERR[,14],ERR[,8],ERR[,19],ERR[,2],what=c(F,T,F,F),log="",col="gray",ylim=RR,
         ylab="True absolute errors and RMSE's",xlab="",xaxt="n")
lines(MM+1,c(MSE[14],MSE[14]),lwd=2,lty=2)
lines(MM+2,c(MSE[ 8],MSE[ 8]),lwd=2,lty=2)
lines(MM+3,c(MSE[19],MSE[19]),lwd=2,lty=2)
lines(MM+4,c(MSE[12],MSE[12]),lwd=2,lty=2)