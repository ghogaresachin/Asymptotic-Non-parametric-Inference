# Practical 1 : Consistency 1
#Q1 Cauchy Distribution
mu=5
sigma=1
eps=0.2
n=c(100,200,500,1000,5000)
est_prob1=est_prob2=0
MSE1=MSE2=0
for (i in 1:length(n)){
  x=matrix(c(rcauchy(n[i]*n[i],mu,sigma)),n[i],n[i])
  T_1=apply(x,2,mean)
  T_2=apply(x,2,median)
  est_prob1[i]=mean(abs(T_1-mu)<eps)
  est_prob2[i]=mean(abs(T_2-mu)<eps)
  MSE1[i]=mean((T_1-mu)^2)
  MSE2[i]=mean((T_2-mu)^2)
}
cbind(n,est_prob1,est_prob2,MSE1,MSE2)

# Q2 Uniform Distribution
a=0
theta=2
eps=0.2
n=c(5,200,500,1000,5000)
est_prob1=est_prob2=est_prob3=0
MSE1=MSE2=MSE3=0
for (i in 1:length(n)){
  x=matrix(c(runif(n[i]*n[i],a,theta)),n[i],n[i])
  T_1=apply(x,2,max)
  T_2=apply(x,2,mean)
  T_3=2*T_2
  est_prob1[i]=mean(abs(T_1-theta)<eps)
  est_prob2[i]=mean(abs(T_2-theta)<eps)
  est_prob3[i]=mean(abs(T_3-theta)<eps)
  MSE1[i]=mean((T_1-theta)^2)
  MSE2[i]=mean((T_2-theta)^2)
  MSE3[i]=mean((T_3-theta)^2)
}
cbind(n,est_prob1,est_prob2,est_prob3,MSE1,MSE2,MSE3)

#Q3) Normal Distribution
mu=3
sigma=1
eps=0.1
n=c(100,200,500,1000,5000)
est_prob1=est_prob2=0
for (i in 1:length(n)){
  x=matrix(c(rnorm(n[i]*n[i],mu,sigma)),n[i],n[i])
  T_1=apply(x,2,mean)
  T_2=apply(x,2,var)
  est_prob1[i]=mean(abs(T_1-mu)<eps)
  est_prob2[i]=mean(abs(T_2-sigma^2)<eps)
}
cbind(n,est_prob1,est_prob2)

#Q4) Emperical Distribution
mu=2.5
eps=0.5
n=c(100,200,500,1000,5000)
t=2
Fxt=1-exp(-t*mu)
est_prob=0
for(i in 1:length(n)){
  x=matrix(c(rexp(n[i]*n[i],1/mu)),n[i],n[i])
  y=x<=t
  T_1=apply(y,1,mean)
  est_prob[i]=mean(abs(T_1-Fxt)<eps)
}
cbind(n,est_prob)

#Q5) Poisson Distribution
lambda=3.5
eps=0.01
n=c(100,200,500,1000,5000)
est_prob1=est_prob2=0
for (i in 1:length(n)){
  x=matrix(c(rpois(n[i]*n[i],lambda)),n[i],n[i])
  x_bar=apply(x,2,mean)
  T_2=exp(-x_bar)
  s=apply(x,2,var)
  T_1=exp(-s)
  est_prob2[i]=mean(abs(T_2-exp(-lambda))<eps)
  est_prob1[i]=mean(abs(T_1-exp(-lambda))<eps)
}
cbind(n,est_prob2,est_prob1)

#Q6) Exponential Distribution(location parameter theta)
a=0
b=1
eps=0.2
theta=2
est_prob1=est_prob2=0
n=c(5,100,200,300,500,1000,5000)
for (i in 1:length(n)){
  x=matrix(c(runif(n[i]*n[i],a,b)),n[i],n[i])
  y=theta-log(1-x)       # or  y=theta+log(1/(1-x))  
  T_1=apply(y,2,min)
  T_2=apply(y,2,mean)-1
  est_prob1[i]=mean(abs(T_1-theta)<eps)
  est_prob2[i]=mean(abs(T_2-theta)<eps)
}
cbind(n,est_prob1,est_prob2)

#Practical 2: Moment Estimator and Consistency
#Q1 Binomial Distribution
rm(list=ls(all=TRUE))
m=5;
p=0.7;
eps=0.1;
n=c(50,200,500,1000,2000);
p1=p2=0
for (i in 1:length(n))
{
  x=matrix(c(rbinom(n[i]*n[i],m,p)),n[i],n[i]);
  m1=apply(x,2,mean)
  m2=apply(x,2,var)
  mhat=m1^2/(m1-m2)
  Phat=(m1-m2)/m1;
  p1[i]=mean(abs(Phat-p)<eps)
  p2[i]=mean(abs(mhat-m)<eps)
}
cbind(n,p1,p2)
par(mfrow=c(2,2))
hist(Phat)
qqnorm(Phat)
hist(mhat)
qqnorm(mhat)

#Q2 Gamma Distribution
alpha=2.5        # Shape Parameter
beta=1.5         # Scale Parameter
eps=0.1
est_prob1=est_prob2=est_prob3=est_prob4=0
est_alpha=est_beta=0
n=c(50,100,200,500,1000)
for(i in 1:length(n)){
  x=matrix(c(rgamma(n[i]*n[i],alpha,1/beta)),n[i],n[i])
  m_1=apply(x,1,mean)   
  m_2=apply(x,1,var)  
  T_1=m_1^2/m_2         # Estimator of alpha
  T_2=m_2/m_1           # Estimator of beta
  est_alpha[i]=mean(T_1)
  est_beta[i]=mean(T_2)
  est_prob3[i]=mean(abs(m_1-alpha*beta)<eps)     # Sample mean is consistent for population mean
  est_prob4[i]=mean(abs(m_2-alpha*beta^2)<eps)
  est_prob1[i]=mean(abs(T_1-alpha)<eps)
  est_prob2[i]=mean(abs(T_2-beta)<eps)
}  
cbind(n,est_alpha,est_beta,est_prob1,est_prob2,est_prob3,est_prob4)
par(mfrow=c(1,2))
hist(T_1)
hist(T_2)


#Q3 Beta distribution of first kind
rm=(list=ls(all=TRUE))
a=2.5         # Shape Parameter
b=1.5         # Scale Parameter
eps=0.1
est_prop1=0
est_prop2=0
est_a=est_b=0
n=c(50,100,200,500,1000)
for(i in 1:length(n))
{
  x=matrix(c(rbeta(n[i]*n[i],a,b)),n[i],n[i])
  m_1=apply(x,1,mean)   
  var=apply(x,1,var)
  m_2=var+m_1^2                             # Second raw moment  
  T_1=(m_1*(m_1-m_2))/(m_2-m_1^2)          # Estimator of a (shape parameter)
  T_2=(m_2-m_1)*(m_1-1)/(m_2-m_1^2)        # Estimator of b (scale parameter)
  est_a[i]=mean(T_1)
  est_b[i]=mean(T_2)
  est_prop1[i]=mean(abs(T_1-a)<eps)
  est_prop2[i]=mean(abs(T_2-b)<eps)
}  
cbind(n,est_a,est_b,est_prop1,est_prop2)
par(mfrow=c(1,2))
hist(T_1)
hist(T_2)
qqnorm(T_1)
qqnorm(T_2)

#Practical 3: Consistency and CAN
rm(list = ls())
#Q1)Normal Distribution
mu=5
sigma=1
n=c(20,30,100,300,500,1000)
est_prob1=est_prob2=est_prob3=0
eps=0.1
for(i in 1:length(n)){
  x=matrix(c(rnorm(n[i]*n[i],mu,sigma)),n[i],n[i])
  T1=apply(x,2,mean)
  T2=apply(x,2,var)
  T3=T1^2+T2^2
  est_prob1[i]=mean(abs(T1-mu)<eps)
  est_prob2[i]=mean(abs(T2-sigma^2)<eps)
  est_prob3[i]=mean(abs(T3-(mu^2+sigma^2))<eps)
}
cbind(n,est_prob1,est_prob2,est_prob3)
par(mfrow=c(2,3))
hist(T1)
qqnorm(T1)
hist(T2)
qqnorm(T2)
hist(T3)
qqnorm(T3)

rm(list=ls())
#Q2Poisson Distribution
lambda=1.5
n=c(20,30,100,300,500,1000)
est_prob1=est_prob2=0
eps=0.1
for(i in 1:length(n)){
  x=matrix(c(rpois(n[i]*n[i],lambda)),n[i],n[i])
  T1=apply(x,2,mean)
  T2=apply(x,2,var)
  est_prob1[i]=mean(abs(T1-lambda)<eps)
  est_prob2[i]=mean(abs(T2-lambda)<eps)
}
cbind(n,est_prob1,est_prob2)
par(mfrow=c(1,2))
hist(T1)
qqnorm(T1)
hist(T2)
qqnorm(T2)

#Q3)Laplace Distribution (When mu=0)
rm(list=ls())
theta=0.8
mu=0
n=c(20,30,100,300,500,1000)
est_prob1=est_prob2=est_prob3=0
eps=0.1
for (i in 1:length(n)){
  y=matrix(c(runif(n[i]*n[i],0,1)),n[i],n[i])
  x=ifelse(y<=0.5,mu+log(2*y)*theta,mu-log(2-2*y)*theta)
  T1=apply(x,2,mean)          
  m2=apply(x,2,var)/2
  T2=sqrt(m2)
  T3=apply(abs(x),2,mean)
  est_prob1[i]=mean(abs(T1-theta)<eps)
  est_prob2[i]=mean(abs(T2-theta)<eps)
  est_prob3[i]=mean(abs(T3-theta)<eps)
}
cbind(n,est_prob1,est_prob2,est_prob3)
par(mfrow=c(2,2))
qqnorm(T2)
hist(T3)
qqnorm(T3)


#Q3 Laplace Distribution
rm(list=ls())
theta=1.5
mu=0.5
n=c(20,30,100,300,500,1000)
est_prob1=est_prob2=est_prob3=est_prob4=0
eps=0.1
for (i in 1:length(n)){
  y=matrix(c(runif(n[i]*n[i],0,1)),n[i],n[i])
  x=ifelse(y<=0.5,mu+log(2*y)*theta,mu-log(2-2*y)*theta)
  T1=apply(x,2,mean)           # MOM
  T2=apply(x,2,median)         # MLE
  T3=sqrt(apply(x,2,var)/2)    # MOM
  T4=apply(abs(x-T2),2,mean)   # MLE
  est_prob1[i]=mean(abs(T1-mu)<eps)
  est_prob2[i]=mean(abs(T2-mu)<eps)
  est_prob3[i]=mean(abs(T3-theta)<eps)
  est_prob4[i]=mean(abs(T4-theta)<eps)
}
cbind(n,est_prob1,est_prob2,est_prob3,est_prob4)
par(mfrow=c(2,4))
hist(T1)
qqnorm(T2)
hist(T2)
qqnorm(T2)
hist(T3)
qqnorm(T3)
hist(T4)
qqnorm(T4)


#Q4 Exponential(location=mu,scale=sigma)
rm(list=ls())
sigma=1.5
mu=0.5
n=c(20,30,100,300,500,1000)
est_prob1=est_prob2=est_prob3=est_prob4=est_prob5=0
eps=0.5
for (i in 1:length(n)){
  y=matrix(c(runif(n[i]*n[i],0,1)),n[i],n[i])
  x=mu-sigma*log(1-y)
  T1=apply(x,1,mean)
  T2=apply(x,1,min)
  T3=apply(x,1,var)
  T4=apply(abs(x-T2),2,mean)
  m2=T1^2+T3
  T5=sqrt(m2/2)
  est_prob1[i]=mean(abs(T1-mu)<eps)     # T1 is not consistent for mu
  est_prob2[i]=mean(abs(T2-mu)<eps)
  est_prob3[i]=mean(abs(T3-sigma)<eps)  # T3 is not consistent for sigma 
  est_prob4[i]=mean(abs(T4-sigma)<eps)
  est_prob5[i]=mean(abs(T5-sigma))
}
cbind(n,est_prob1,est_prob2,est_prob3,est_prob4,est_prob5)
par(mfrow=c(2,2))
hist(T2)
qqnorm(T2)
hist(T3)
qqnorm(T3)
hist(T4)
qqnorm(T4)

# Practical 4:Invariance property of consistent CAN and percentile estimators
# Normal
rm(list=ls(all=TRUE))
mu=2.1;sigma=1.35;
mu2=mu^2; s2=sigma^4
eps=0.1;
V=sigma^2;
n=c(50,200,500,1000,2000,3000,4000);
p1=p2=p3=0
for (i in 1:length(n))
{
  x=matrix(c(rnorm(n[i]*n[i],mu,sigma)),n[i],n[i])
  m=apply(x,2,mean)
  T1=m^2;
  v=(apply(x,2,var))
  T2=v^2
  T3=T1+T2
  p1[i]=mean(abs(T1-mu2)<eps)
  p2[i]=mean(abs(T2-s2)<eps)
  p3[i]=mean(abs(T3-(mu2+s2))<eps)
}
cbind(n,p1,p2,p3)
par(mfrow=c(2,3))
hist(T1)
hist(T2)
hist(T3)
qqnorm(T1)
qqnorm(T2)
qqnorm(T3)


# Poisson
rm(list=ls(all=TRUE))
lam=2.4;
eps=0.01;
a=exp(-lam);
b=lam*exp(-lam)
n=c(50,200,500,1000,2000);
p1=p2=0
for (i in 1:length(n)){
  x=matrix(c(rpois(n[i]*n[i],lam)),n[i],n[i]);
  m=apply(x,2,mean);
  T1=exp(-m);
  T2=m*T1;
  p1[i]=mean(abs(T1-a)<eps)
  p2[i]=mean(abs(T2-b)<eps)
}
cbind(n,p1,p2)
par(mfrow=c(2,2))
hist(T1)
hist(T2)
qqnorm(T1)
qqnorm(T2)


# Binomial
rm(list=ls(all=TRUE))
m=1;
p1=0.35;q1=1-p1;
eps=0.01;
a=p1*(1-p1);
n=c(50,200,500,1000,2000);
p=matrix(rep(0,3*length(n)),ncol=length(n))
for (i in 1:length(n))
{
  x=matrix(c(rbinom(n[i]*n[i],m,p1)),n[i],n[i]);
  T1=1-(apply(x,2,mean));
  T2=apply(x,2,mean)
  T3=T2*T1;
  p[1,i]=mean(abs(T1-q1)<eps);
  p[2,i]=mean(abs(T2-p1)<eps);
  p[3,i]=mean(abs(T3-a)<eps);
}
par(mfrow=c(2,3))
hist(T1)
hist(T2)
hist(T3)
qqnorm(T1)
qqnorm(T2)
qqnorm(T3)

# Binomial(m,p)
rm(list=ls(all=TRUE))
p=0.5
m=10
n=c(5,20,50,200,500,1000,2000)
a=(1-p)^m
b=m*p*(1-p)^(m-1)
c=p*(1-p)
d=m*(m-1)/2*p^2*(1-p)^(m-2)
e=((1-p)/p)^2
eps=0.01
p1=p2=p3=p4=p5=0
for (i in 1:length(n)) {
  x=matrix(c(rbinom(n[i]*n[i],m,p)),n[i],n[i])
  m1=apply(x,1,mean)
  m2=apply(x,1,var)
  mhat=m1^2/(m1-m2)
  phat=(m1-m2)/m1
  T1=(1-phat)^mhat
  T2=mhat*phat*(1-phat)^(mhat-1)
  T3=phat*(1-phat)
  T4=mhat*(mhat-1)/2*phat^2*(1-phat)^(mhat-2)
  T5=((1-phat)/phat)^2
  p1[i]=mean(abs(T1-a)<eps)
  p2[i]=mean(abs(T2-b)<eps)
  p3[i]=mean(abs(T3-c)<eps)
  p4[i]=mean(abs(T4-d)<eps)
  p5[i]=mean(abs(T5-e)<eps)
}
cbind(n,p1,p2,p3,p4,p5)
par(mfrow=c(4,3))
hist(T1)
hist(T2)
hist(T3)
hist(T4)
hist(T5)
qqnorm(T1)
qqnorm(T2)
qqnorm(T3)
qqnorm(T4)
qqnorm(T5)

#Practical 5: Method of Scoring
# Scoring
#Q1) Normal Distribution
rm(list=ls(all=TRUE))
n=25; mu=4;sigma=1;
x=rnorm(n,mu,1)
u=c(median(x),rep(0,10))
v=c(median(x),rep(0,10))
I=1/sigma^2;
for (i in 1:10)
{
  d_1=n*(mean(x)-u[i]) #dlogl/du
  d2_1=-n #d^2logl/du^2
  u[i+1]=u[i]-d_1/d2_1 #using newton raphson method
  v[i+1]=u[i]+d_1/(n*I) #using fisher scoring method
}u
v

#Q2) Cauchy Distribution
n=25
x=rcauchy(n,0,1)
u=c(median(x),rep(0,10))  # Newton Raphson Method
v=c(median(x),rep(0,10))  # Fisher Scoring Method
I=1/2
for ( i in 1:10){
  z=x-u[i]
  y=x-v[i]
  d_1=2*sum(z/(1+z^2))                  # dlogl/du
  d2_1=2*sum((2*z^2-1-z^2)/(1+z^2)^2)   # d^2logl/du^2
  d_1_=2*sum(y/(1+y^2)) 
  u[i+1]=u[i]-d_1/d2_1                  # Newton Raphson Method
  v[i+1]=v[i]+d_1_/(n*I)                # Fisher Scoring Method
}
u
v

#Normal(theta,theta^2)
rm(list=ls(all=TRUE))
n=25; th=2.8;
x=rnorm(n,th,th)
u=c(mean(x),rep(0,10))
v=c(mean(x),rep(0,10))
for (i in 1:10)
{
  I=3/u[i]^2;
  d_1=(-n/u[i])+(1/u[i]^3)*sum(x^2)-(1/u[i]^2)*sum(x) #dlogl/du
  d2_1=n/u[i]^2-(3/u[i]^4)*sum(x^2)+(2/u[i]^3)*sum(x); #d^2logl/du^2
  u[i+1]=u[i]-d_1/d2_1 #using newton raphson method
  v[i+1]=u[i]+d_1/(n*I) #using fisher scoring method
}
u
v


# Multinomial Distribution
n1=1977
n2=906
n3=904
n4=32
n=n1+n2+n3+n4
th=c((1-4*n2/n),rep(0,9))
for ( i in 1:10){
  I=(2*th[i]+1)/(2*th[i]*(1-th[i])*(2+th[i]))   # dlogl/dtheta
  d_1=(n1/(2+th[i]))-(n2+n3)/(1-th[i])+n4/th[i] # I_x(theta)
  th[i+1]=th[i]+d_1/(n*I)                       # Method of scoring
}
th

# Cauchy(theta,1)
n=25
x=rcauchy(n,4,1)
u=c(median(x),rep(0,10))  # Newton Raphson Method
v=c(median(x),rep(0,10))  # Fisher Scoring Method
I=1/2
for ( i in 1:10){
  z=x-u[i]
  y=x-v[i]
  d_1=2*sum(z/(1+z^2))                  # dlogl/du
  d_1_=2*sum(y/(1+y^2)) 
  d2_1=2*sum((2*z^2-1-z^2)/(1+z^2)^2)   # d^2logl/du^2
  u[i+1]=u[i]-d_1/d2_1                  # Newton Raphson Method
  v[i+1]=v[i]+d_1_/(n*I)                # Fisher Scoring Method
}
u
v
z=qnorm(0.05)
LCI=u[11]-z*(sqrt(2/n))
UCI=u[11]+z*(sqrt(2/n))



n=25
mu=4
sigma=1
x=rnorm(25,mu,sigma)
u=c(median(x),rep(0,10))
v=c(median(x),rep(0,10))
I=1/sigma^2
for (i in 1:11){
  d1=n*(mean(x)-u[i])
  d2=-n
  u[i+1]=u[i]-d1/d2
  v[i+1]=v[i]+d1/(n*I)
}
u
v

# Practical 6: ACI
# Poisson
lambda=2
n=20
x=rpois(n,lambda)
xbar=mean(x)
z=qnorm(0.025,lower.tail = F)
# Pivotal Method
L=(2*xbar+(z^2)/n-sqrt((2*xbar*z^2)/n+z^2/n^2))/2
U=(2*xbar+(z^2)/n+sqrt((2*xbar*z^2)/n+z^2/n^2))/2
cbind(L,U)
# VST Method
L1=(sqrt(xbar)-z/(2*sqrt(n)))^2
U1=(sqrt(xbar)+z/(2*sqrt(n)))^2
cbind(L1,U1)

# Binomial
n=25
x=rbinom(n,1,prob = 0.4)
xbar=mean(x)
z=qnorm(0.025,lower.tail = F)
# ACI
l1=((2*n*xbar+z^2)-sqrt(4*z^2*n*xbar*(1-xbar)+z^2/2))/(2*z^2+n)
u1=((2*n*xbar+z^2)+sqrt(4*z^2*n*xbar*(1-xbar)+z^2/2))/(2*z^2+n)
cbind(l1,u1)

# Laplace
n=25
y=runif(n,0,1)
x=ifelse(y<=0.5,mu+log(2*y),mu-log(2-2*y))
xbar=mean(x)
med=median(x)
z=abs(qnorm(0.025,lower.tail = F))
# ACI based on mean
L=xbar-z*sqrt(2/n)
U=xbar+z*sqrt(2/n)
cbind(L,U)
# Median
L=med-z*sqrt(1/n)
U=med+z*sqrt(1/n)
cbind(L,U)

# Exponential Distribution
n=25
theta=3
x=rexp(n,1/theta)
xbar=mean(x)
#ACI using Pivotal Method
l_=xbar/(1+z/sqrt(n))
u_=xbar/(1-z/sqrt(n))
cbind(l_,u_)
# ACI using VST method
l_2=exp(log(xbar)-z*sqrt(1/n))
u_2=exp(log(xbar)+z*sqrt(1/n))
cbind(l_2,u_2)

