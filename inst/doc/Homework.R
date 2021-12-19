## ----echo=TRUE, fig.align='center', fig.height=5, fig.width=5-----------------
x = 1:10
y = runif(10)
plot(x,y,type = "p")

## ----echo=TRUE,results='hide'-------------------------------------------------
x = 1:20; y = x^2;
lmr = lm(y ~ x)
Co = summary(lmr)$coefficients
print(Co)

## ----echo=TRUE----------------------------------------------------------------
knitr::kable(Co,align = 'c')

## -----------------------------------------------------------------------------

for(i in c(0.5,1,2,5,10))
{
  sigma = i #参数sigma取不同的值
  N = 10000
  u = runif(N,min=0,max=1)
  x = sqrt(-2*sigma^2*log(1-u)) #模拟生成服从Rayleigh分布的随机变量
  hist(x,freq = F) #这里画出直方图，注意纵坐标显示的是频率，因为freq = F
  lines(density(x),col = "red") #这里画出密度函数，也可以观察出生成的随机变量是否具有Rayleigh分布的形式
}


## -----------------------------------------------------------------------------
for(p in c(0.1,0.2,0.4,0.5,0.6,0.75,0.8,0.9)) #题目中给定的p1=0.75也放在循环里了
{
    ar = matrix(0,1,1000)
  for(i in 1:1000) #生成1000个服从N(0,1)和N(3,1)混合分布的随机变量
  {
    u = runif(1,min=0,max=1)
    if(u<p)
      x = rnorm(1,0,1) else
        x = rnorm(1,3,1)
    ar[i] = x
  }
  hist(ar,freq = F,main = "样本直方图") #画出样本直方图
  lines(density(ar),col = "red") #画出样本密度函数
  plot.ecdf(ar,main = "样本经验分布函数") #画出样本经验函数
}


## -----------------------------------------------------------------------------
lamda = 6
alpha = 1
beta = 2 #明确参数的取值
X = matrix(0,1,1000)
for(i in 1:1000) #重复1000次模拟
{
  N = rpois(1,10*lamda) #随机生成N(10)
  y = rgamma(N,shape = alpha,scale = 1/beta) #随机生成Y1,Y2,...,YN,均服从Gamma(1,2)
  X[i] = sum(y) #Y1,Y2,...,YN的和为X(10)的一个模拟值
}
mean = mean(X)
var = sd(X)^2
cat("估计的均值和方差分别为",mean,var)

Mean = 10*lamda*alpha*(1/beta)
Var = 10*lamda*alpha*(alpha+1)*(1/beta^2)
cat("理论的均值和方差分别为",Mean,Var)

## -----------------------------------------------------------------------------
lamda = 10
alpha = 2
beta = 1
X = matrix(0,1,1000)
for(i in 1:1000)
{
  N = rpois(1,10*lamda)
  y = rgamma(N,shape = alpha,scale = 1/beta)
  X[i] = sum(y)
}
mean = mean(X)
var = sd(X)^2
cat("估计的均值和方差分别为",mean,var)

Mean = 10*lamda*alpha*(1/beta)
Var = 10*lamda*alpha*(alpha+1)*(1/beta^2)
cat("理论的均值和方差分别为",Mean,Var)

## -----------------------------------------------------------------------------
lamda = 1
alpha = 3
beta = 4
X = matrix(0,1,1000)
for(i in 1:1000)
{
  N = rpois(1,10*lamda)
  y = rgamma(N,shape = alpha,scale = 1/beta)
  X[i] = sum(y)
}
mean = mean(X)
var = sd(X)^2
cat("估计的均值和方差分别为",mean,var)

Mean = 10*lamda*alpha*(1/beta)
Var = 10*lamda*alpha*(alpha+1)*(1/beta^2)
cat("理论的均值和方差分别为",Mean,Var)

## -----------------------------------------------------------------------------
lamda = 2
alpha = 8
beta = 10
X = matrix(0,1,1000)
for(i in 1:1000)
{
  N = rpois(1,10*lamda)
  y = rgamma(N,shape = alpha,scale = 1/beta)
  X[i] = sum(y)
}
mean = mean(X)
var = sd(X)^2
cat("估计的均值和方差分别为",mean,var)

Mean = 10*lamda*alpha*(1/beta)
Var = 10*lamda*alpha*(alpha+1)*(1/beta^2)
cat("理论的均值和方差分别为",Mean,Var)

## -----------------------------------------------------------------------------
MC_cdf = function(x,a,b) {
  u = runif(5000,min = 0,max = x)
  f = gamma(a+b)/(gamma(a)*gamma(b))*u^(a-1)*(1-u)^(b-1) #这是贝塔分布的密度函数，用上一步生成的u来计算f(u)
  cdf = x*mean(f) #估计F(x)
}

x = seq(0.1,0.9,0.1)
for (i in x)
{
  MC_estimate = MC_cdf(i,3,3) #估计Beta(3,3)的cdf
  pbeta_value = pbeta(i,3,3)
  print(c(i,MC_estimate,pbeta_value)) #分别输出x的值,Monte Carlo估计值和pbeta返回的值
}

## -----------------------------------------------------------------------------
# 设计一个函数，可以利用对偶变量生成服从Rayleigh分布的随机变量
Rayleigh_Sample = function(sigma){
  u = runif(1,min=0,max=1)
  x1 = sqrt(-2*sigma^2*log(1-u))
  x2 = sqrt(-2*sigma^2*log(u))
  return((x1+x2)/2)
}
Rayleigh_Sample(0.5) # 取不同的参数sigma，得到生成的随机变量
Rayleigh_Sample(1)
Rayleigh_Sample(2)
Rayleigh_Sample(5)

## -----------------------------------------------------------------------------
# 计算使用不同方法（对偶变量、独立变量）下样本的方差及减少的比例
Rayleigh_Sample_Var = function(sigma,n,antithetic = TRUE){
  u = runif(n)
  if (antithetic){
    X1 = sqrt(-2*sigma^2*log(1-u))
    X2 = sqrt(-2*sigma^2*log(u))
    variance = (var(X1)+var(X2)+2*cov(X1,X2))/4 #计算(X1+X2)/2的样本方差，此时X1和X2负相关
  }
  else{
    v = runif(n)
    X1 = sqrt(-2*sigma^2*log(1-u))
    X2 = sqrt(-2*sigma^2*log(1-v))
    variance = (var(X1)+var(X2))/4 #计算(X1+X2)/2的样本方差，此时X1和X2独立
  }
variance
}

var1 = Rayleigh_Sample_Var(1,1000,antithetic = FALSE) #取参数sigma=1，样本量n=1000（都是任意的），计算样本方差
var2 = Rayleigh_Sample_Var(1,1000,antithetic = TRUE)
ratio = 100*(var1-var2)/var1 #计算方差减少的比例（百分比）
print(c(var1,var2,ratio)) #输出结果

## -----------------------------------------------------------------------------
x = seq(1,10,0.05)
g = (1/(sqrt(2*pi)))*x^2*exp(-x^2/2)
f1 = exp(-x) 
f2 = 1/x^2 #挑选g(x)的两个importance function

# 在同一张图上画出g,f1,f2的密度函数
plot(x,g,type = "l",ylab = "",ylim = c(0,0.5),lwd = 2)
lines(x, f1, lty = 2, lwd = 2, col = "red")
lines(x, f2, lty = 4, lwd = 2, col = "blue")
legend("topright",legend = c("g","f1","f2"),lty = c(1,2,4),col = c("black","red","blue"),lwd = 2,inset = 0.05)
# 在同一张图上画出g与f1,f2的比例
plot(x,g/f1,type = "l",ylab = "",ylim = c(0,2),lty = 2,lwd = 2,col = "red")
lines(x,g/f2,lty = 4,lwd = 2,col = "blue")
legend("topright",legend = c("f1","f2"),lty = c(2,4),col = c("red","blue"),lwd = 2,inset = 0.05)

n = 5000 
se = numeric(2)
g = function(x){
  (1/(sqrt(2*pi)))*x^2*exp(-x^2/2)*(x > 1)
}
# 用f1作为importance function
u = runif(n,0,1)
x = -log(1-u) #满足密度函数为f1
fg = g(x)/exp(-x)
se[1] = sd(fg)
# 用f2作为importance function
u = runif(n,0,1)
x = 1/(1-u) #满足密度函数为f2
fg = g(x)/(1/x^2)
se[2] = sd(fg)
print(round(se,4)) #输出估计量的标准差

## -----------------------------------------------------------------------------
g = function(x){
  (1/(sqrt(2*pi)))*x^2*exp(-x^2/2)*(x > 1)
}
n = 10000                 
u = runif(n,0,1)
x = 1/(1-u) #X的密度函数为f2 = 1/x^2
fg = g(x)/(1/x^2)
MC_est = mean(fg)
MC_est
Real_value = integrate(g,1,Inf)
Real_value

## -----------------------------------------------------------------------------
# Use a Monte Carlo experiment to estimate the coverage probability of the t-interval
set.seed(111)
n = 20
alpha = 0.05
m = 1000 #重复次数
LCL = numeric(m)
UCL = numeric(m)
for(i in 1:m)
{
  x = rchisq(n, 2)
  LCL[i] = mean(x) + qt(alpha/2,df=n-1)*sd(x)/sqrt(n) #注意qt中的lower.tail默认为TRUE,意为P(X ≤ x),所以这里是+号,计算UCL时为-号
  UCL[i] = mean(x) - qt(alpha/2,df=n-1)*sd(x)/sqrt(n)
}
#计算coverage probability(CP)
sum = 0
for (i in 1:1000) 
{
  if(LCL[i] <= 2 && UCL[i] >= 2) #表示χ2(2)的理论均值2落在CI内
    sum = sum + 1
}
CP = sum/1000
CP

# Example 6.4(Confidence interval for variance) 用例子中的方法估计方差的CP
UCL_var = numeric(m)
for(i in 1:m)
{
  y = rchisq(n, 2)
  UCL_var[i] = (n-1)*var(y)/qchisq(alpha,df = n-1)
}
sum_var = 0
for (i in 1:1000) 
{
  if(UCL_var[i] > 4) #表示χ2(2)的理论方差4落在CI内
    sum_var =sum_var + 1
}
CP_var = sum_var/1000
CP_var

## ----eval=FALSE---------------------------------------------------------------
#  n = c(10,20,50,100,200)
#  L = length(n)
#  p1_hat = numeric(L)
#  p2_hat = numeric(L)
#  p3_hat = numeric(L)
#  alpha = 0.05
#  mu0 = 1
#  m = 10000
#  p1 = numeric(m)
#  p2 = numeric(m)
#  p3 = numeric(m)
#  
#  # χ2(1)
#  for(i in 1:L)
#  {
#    num = n[i]
#    for(j in 1:m)
#    {
#      x = rchisq(num,1)
#      ttest1 = t.test(x,mu = mu0)
#      p1[j] = ttest1$p.value
#    }
#    p1_hat[i] = mean(p1 <= alpha)
#  }
#  cat("在抽样总体为χ2(1)时,当n分别为10,20,50,100,200,Type I error rate分别为",p1_hat,"\n")
#  
#  # U(0,2)
#  for(i in 1:L)
#  {
#    num = n[i]
#    for(j in 1:m)
#    {
#      x = runif(num,0,2)
#      ttest2 = t.test(x,mu = mu0)
#      p2[j] = ttest2$p.value
#    }
#    p2_hat[i] = mean(p2 <= alpha)
#  }
#  cat("在抽样总体为U(0,2)时,当n分别为10,20,50,100,200,Type I error rate分别为",p2_hat,"\n")
#  
#  # Exp(1)
#  for(i in 1:L)
#  {
#    num = n[i]
#    for(j in 1:m)
#    {
#      x = rexp(num,1)
#      ttest3 = t.test(x,mu = mu0)
#      p3[j] = ttest3$p.value
#    }
#    p3_hat[i] = mean(p3 <= alpha)
#  }
#  cat("在抽样总体为Exp(1)时,当n分别为10,20,50,100,200,Type I error rate分别为",p3_hat,"\n")
#  

## ----echo=FALSE---------------------------------------------------------------
cat("Question 1: What is the corresponding hypothesis test problem?","\n")
cat("Answer 1: The null hypothesis is the powers for two methods are same, while the alternative hypothesis is the powers for two methods are different.","\n")

cat("Question 2: What test should we use? Z-test, two-sample t-test, paired-t test or McNemar test? Why?","\n")
cat("Answer 2: Maybe we should use the McNemar test, because the McNemar test is used to determine whether the proportions of categories in two related groups significantly differ from each other, while the Z-test is used to determine whether two population means are different when the variances are known, the two-sample t-test is used to determine whether the means of two random samples of independent observations are equal or not, and the paired-t test is used to determine whether the mean difference between pairs of measurements is zero.","\n")
 
cat("Question 3: Please provide the least necessary information for hypothesis testing.","\n")
cat("Answer 3: If we use McNemar test, at least we need to know the two things, the first is the probability(or the number of incidents, denoted as n1) when Method 1 is positive while Method 2 is negative, denoted as p1, the second is the probability(or the number of incidents, denoted as n2) when Method 1 is negative while Method 2 is positive, denoted as p2. At this time, the null hypothesis is p1 = p2, and the McNemar test statistic (n1 - n2)^2/(n1 + n2) has a chi-squared distribution with 1 degree of freedom, if the result is significant, this provides sufficient evidence to reject the null hypothesis, which means that the marginal proportions are significantly different from the powers for two methods.","\n")

## ----eval=FALSE---------------------------------------------------------------
#  # Skewness test
#  # 不妨就使X和Y均从μ=(0,0,0)',∑=I(单位阵)的多元高斯分布N(μ,∑)中生成
#  library(MASS)
#  Mardia_test = function(data){
#    n = nrow(data)
#    d = ncol(data)
#    central = matrix(0,n,d)
#    for(i in 1:d){
#      central[,i] = data[,i]-mean(data[,i])
#    }
#    sigma_hat = t(central)%*%central/n
#    a = central%*%solve(sigma_hat)%*%t(central)
#    b = sum(colSums(a^{3}))/(n*n)
#    test = n*b/6
#    chi = qchisq(0.95,d*(d+1)*(d+2)/6)
#    as.integer(test > chi)
#  }
#  
#  set.seed(123)
#  mu = c(0,0,0)
#  sigma = matrix(c(1,0,0,0,1,0,0,0,1),nrow = 3,ncol = 3)
#  m = 10000
#  n = c(10, 20, 30, 50, 100, 500)
#  p.reject = numeric(length(n))
#  for(i in 1:length(n)){
#    p.reject[i] = mean(replicate(m, expr={
#      data = mvrnorm(n[i],mu,sigma) #产生服从多元正态分布的随机数
#      Mardia_test(data)
#    }))
#  }
#  print(p.reject)

## ----eval=FALSE---------------------------------------------------------------
#  # Power of the skewness test
#  # 不妨就取μ1=μ2=(0,0,0)',∑1=I(单位阵),∑2=diag(100,100,100)(对角阵),且contaminated multivariate normal distribution为(1−ε)N(μ1,∑1)+εN(μ2,∑2)
#  library(MASS)
#  set.seed(123)
#  mu1 = c(0,0,0)
#  mu2 = c(0,0,0)
#  sigma1 = matrix(c(1,0,0,0,1,0,0,0,1),nrow = 3,ncol = 3)
#  sigma2 = matrix(c(100,0,0,0,100,0,0,0,100),nrow = 3,ncol = 3)
#  sigma = list(sigma1,sigma2)
#  m = 2500
#  n = 30
#  epsilon = c(seq(0, 0.15, 0.01), seq(0.15, 1, 0.05))
#  N = length(epsilon)
#  pwr = numeric(N)
#  for (j in 1:N){
#    e = epsilon[j]
#    sktests = numeric(m)
#    for (i in 1:m){
#      index= sample(c(1, 2), replace = TRUE, size = n, prob = c(1-e, e))
#      data = matrix(0,nrow = n,ncol = 3)
#      for(t in 1:n){
#        if(index[t] == 1) data[t,] = mvrnorm(1,mu1,sigma1) #从N(μ1,∑1)中生成
#        else data[t,] = mvrnorm(1,mu2,sigma2) #从N(μ2,∑2)中生成
#      }
#      sktests[i] = Mardia_test(data) #利用上述Skewness test中的检验
#    }
#    pwr[j] = mean(sktests)
#  }
#  plot(epsilon, pwr, type = "b", xlab = bquote(epsilon), ylim = c(0,1))
#  abline(h = 0.1, lty = 3)
#  se = sqrt(pwr * (1-pwr) / m)
#  lines(epsilon, pwr+se, lty = 3)
#  lines(epsilon, pwr-se, lty = 3)

## -----------------------------------------------------------------------------
library(boot)
library(bootstrap)
B = 1000
n = nrow(scor)
lambda_hat = eigen((n-1)/n*cov(scor))$values #极大似然估计
theta_hat = lambda_hat[1] / sum(lambda_hat)


set.seed(111)
func = function(data,i){
  x = data[i,]
  lambda = eigen((n-1)/n*cov(x))$values
  theta = lambda[1] / sum(lambda)
  return(theta)
}

bootstrap_result = boot(
  data = cbind(scor$mec,scor$vec,scor$alg,scor$ana,scor$sta),
  statistic = func,R = B)
theta_b = bootstrap_result$t

bias_boot = mean(theta_b)-theta_hat
se_boot = sqrt(var(theta_b))
# 输出结果
cat('The Bootstrap estimate of bias is :',bias_boot,'\n')
cat('The Bootstrap estimate of standard error is :',se_boot,'\n')

## -----------------------------------------------------------------------------
theta_j = rep(0,n)
for(i in 1:n)
{
  x = scor[-i,]
  lambda = eigen((n-2)/(n-1)*cov(x))$values
  theta_j[i] = lambda[1] / sum(lambda) 
}
bias_jack = (n-1) * (mean(theta_j) - theta_hat)
se_jack = sqrt((n-1) * mean((theta_j - mean(theta_j))^2))
# 输出结果
cat('The Jackknife estimate of bias is :',bias_jack,'\n')
cat('The Jackknife estimate of standard error is :',se_jack,'\n')

## -----------------------------------------------------------------------------
boot.ci(bootstrap_result, conf = 0.95, type = c('perc','bca'))

## ----eval=FALSE---------------------------------------------------------------
#  # estimate the standard normal bootstrap confidence interval, the basic bootstrap confidence interval, and the percentile confidence interval
#  n = 50 #样本量
#  m = 5000
#  B = 1000
#  set.seed(123)
#  CI_norm1 = matrix(0,m,2)
#  CI_basic1 = matrix(0,m,2)
#  CI_perc1 = matrix(0,m,2) #用来存放 normal population 的结果
#  CI_norm2 = matrix(0,m,2)
#  CI_basic2 = matrix(0,m,2)
#  CI_perc2 = matrix(0,m,2) #用来存放 χ2(5) distribution 的结果
#  
#  skewness = function(data,i){
#    x = data[i]
#    xbar = mean(x)
#    m3 = mean((x - xbar)^3)
#    m2 = mean((x - xbar)^2)
#    return(m3 / m2^1.5)
#  }
#  
#  for(i in 1:m)
#  {
#    x1 = rnorm(n)
#    boot_result1 = boot(data = x1,statistic = skewness,R = B)
#    boot_CI1 = boot.ci(boot_result1,type=c('norm','basic','perc'))
#    CI_norm1[i,]<-boot_CI1$norm[2:3]
#    CI_basic1[i,]<-boot_CI1$basic[4:5]
#    CI_perc1[i,]<-boot_CI1$perc[4:5]
#  }
#  for(i in 1:m)
#  {
#    x2 = rchisq(n,5)
#    boot_result2 = boot(data = x2,statistic = skewness,R = B)
#    boot_CI2 = boot.ci(boot_result2,type=c('norm','basic','perc'))
#    CI_norm2[i,]<-boot_CI2$norm[2:3]
#    CI_basic2[i,]<-boot_CI2$basic[4:5]
#    CI_perc2[i,]<-boot_CI2$perc[4:5]
#  }
#  # use confidence interval to estimate the coverage probabilities
#  mu1 = 0 #N(0,1)的理论偏度为0
#  mu2 = sqrt(8/5) #χ2(5)的理论偏度为sqrt(8/5)
#  CP1_norm = mean(CI_norm1[,1]<=mu1 & CI_norm1[,2]>=mu1)
#  CP1_basic = mean(CI_basic1[,1]<=mu1 & CI_basic1[,2]>=mu1)
#  CP1_perc = mean(CI_perc1[,1]<=mu1 & CI_perc1[,2]>=mu1)
#  # 以上,计算三种CI下N(0,1)的CP
#  CP2_norm = mean(CI_norm2[,1]<=mu2 & CI_norm2[,2]>=mu2)
#  CP2_basic = mean(CI_basic2[,1]<=mu2 & CI_basic2[,2]>=mu2)
#  CP2_perc = mean(CI_perc2[,1]<=mu2 & CI_perc2[,2]>=mu2)
#  # 以上,计算三种CI下χ2(5)的CP
#  # 下面输出CP的结果
#  CP1 = cbind(CP1_norm,CP1_basic,CP1_perc)
#  CP2 = cbind(CP2_norm,CP2_basic,CP2_perc)
#  result = rbind(CP1,CP2)
#  colnames(result) = c("norm","basic","perc")
#  rownames(result) = c("N(0,1)","χ2(5)")
#  print(result)

## ----eval=FALSE---------------------------------------------------------------
#  # Find the proportion of times that the confidence intervals miss on the left, and the porportion of times that the confidence intervals miss on the right
#  left1 = cbind(mean(CI_norm1[,1]>mu1),mean(CI_basic1[,1]>mu1),mean(CI_perc1[,1]>mu1))
#  right1 = cbind(mean(CI_norm1[,2]<mu1),mean(CI_basic1[,2]<mu1),mean(CI_perc1[,2]<mu1))
#  # 以上,计算N(0,1)的miss on the left and right
#  left2 = cbind(mean(CI_norm2[,1]>mu2),mean(CI_basic2[,1]>mu2),mean(CI_perc2[,1]>mu2))
#  right2 = cbind(mean(CI_norm2[,2]<mu2),mean(CI_basic2[,2]<mu2),mean(CI_perc2[,2]<mu2))
#  # 以上,计算χ2(5)的miss on the left and right
#  # 下面输出miss的结果
#  miss = rbind(left1,right1,left2,right2)
#  colnames(miss) = c("norm","basic","perc")
#  rownames(miss) = c("left loss of N(0,1)","right loss of N(0,1)","left loss of χ2(5)","right loss of χ2(5)")
#  print(miss)

## -----------------------------------------------------------------------------
set.seed(1)
x = rnorm(20,1,2)
y = rt(20,df=1) # Generate two groups of random numbers
#null hypothesis: true rho is equal to 0
cor.test = cor.test(x,y,method = "spearman")
rho = cor.test$estimate
p = cor.test$p.value

R = 10000
z = c(x, y)
K = length(z)
rhos = numeric(R)

for (i in 1:R) {
  k = sample(1:K,size = K/2,replace = FALSE)
  x1 = z[k]
  y1 = z[-k]
  rhos[i] = cor(x1,y1,method = "spearman")
}

hist(rhos,breaks = 100)
p1 = mean(abs(rhos) > rho)
cat('The p-value reported by cor.test is:',p,'\n')
cat('The achieved significance level of the permutation test is:',p1,'\n')

## ----eval=FALSE---------------------------------------------------------------
#  library(RANN)
#  library(energy)
#  library(Ball)
#  library(boot)
#  set.seed(12345)
#  m = 50
#  k = 3
#  p = 2
#  n1 = 20
#  n2 = 20
#  N = c(n1,n2)
#  R = 1000
#  
#  Tn = function(z,ix,sizes,k) {
#    n1 = sizes[1]
#    n2 = sizes[2]
#    n = n1 + n2
#    if(is.vector(z)) z = data.frame(z,0)
#    z = z[ix, ]
#    NN = nn2(data = z,k = k + 1)
#    block1 = NN$nn.idx[1:n1,-1]
#    block2 = NN$nn.idx[(n1+1):n,-1]
#    i1 = sum(block1 < n1 + 0.5)
#    i2 = sum(block2 > n1 + 0.5)
#    return((i1 + i2) / (k * n))
#  }
#  
#  eqdist.nn = function(z,sizes,k) {
#    boot.obj = boot(data = z,statistic = Tn,R = R,sim = "permutation",sizes = sizes,k = k)
#    ts = c(boot.obj$t0,boot.obj$t)
#    p_value = mean(ts >= ts[1])
#    list(statistic = ts[1],p.value = p_value)
#  }
#  
#  p_values = matrix(NA,m,3)

## ----eval=FALSE---------------------------------------------------------------
#  ### Unequal variances and equal expectations
#  set.seed(12345)
#  for(i in 1:m){
#    x = matrix(rnorm(n1*p,0,1),ncol = p)
#    y = matrix(rnorm(n2*p,0,2),ncol = p)
#    z = rbind(x,y)
#    p_values[i,1] = eqdist.nn(z,N,k)$p.value #NN
#    p_values[i,2] = eqdist.etest(z,sizes = N,R = R)$p.value #energy
#    p_values[i,3] = bd.test(x = x,y = y,R = 1000,seed=i*12345)$p.value #Ball
#  }
#  alpha = 0.05
#  pow = colMeans(p_values < alpha)
#  data.frame(methods = c('NN','energy','Ball'),pow)

## ----eval=FALSE---------------------------------------------------------------
#  ### Unequal variances and unequal expectations
#  set.seed(12345)
#  for(i in 1:m){
#    x = matrix(rnorm(n1*p,0,1),ncol = p)
#    y = matrix(rnorm(n2*p,1,2),ncol = p)
#    z = rbind(x,y)
#    p_values[i,1] = eqdist.nn(z,N,k)$p.value #NN
#    p_values[i,2] = eqdist.etest(z,sizes = N,R = R)$p.value #energy
#    p_values[i,3] = bd.test(x = x,y = y,R = 1000,seed=i*12345)$p.value #Ball
#  }
#  alpha = 0.05
#  pow = colMeans(p_values < alpha)
#  data.frame(methods = c('NN','energy','Ball'),pow)

## ----eval=FALSE---------------------------------------------------------------
#  ### Non-normal distributions: t distribution with 1 df (heavy-tailed distribution), bimodel distribution (mixture of two normal distributions)
#  set.seed(12345)
#  for(i in 1:m){
#    x = matrix(rt(n1*p,df = 1),ncol = p)
#    y1 = rnorm(n2*p,1,1)
#    y2 = rnorm(n2*p,1,2)
#    w = rbinom(n2*p,1,0.5)
#    y = matrix(w*y1 + (1-w)*y2,ncol = p) #bimodel distribution
#    z = rbind(x,y)
#    p_values[i,1] = eqdist.nn(z,N,k)$p.value #NN
#    p_values[i,2] = eqdist.etest(z,sizes = N,R = R)$p.value #energy
#    p_values[i,3] = bd.test(x = x,y = y,R = 1000,seed=i*12345)$p.value #Ball
#  }
#  alpha = 0.05
#  pow = colMeans(p_values < alpha)
#  data.frame(methods = c('NN','energy','Ball'),pow)

## ----eval=FALSE---------------------------------------------------------------
#  ### Unbalanced samples
#  set.seed(12345)
#  n1 = 20
#  n2 = 50
#  N = c(n1,n2)
#  for(i in 1:m){
#    x = matrix(rt(n1*p,df = 1),ncol = p)
#    y = matrix(rnorm(n2*p,0,1),ncol = p)
#    z = rbind(x,y)
#    p_values[i,1] = eqdist.nn(z,N,k)$p.value #NN
#    p_values[i,2] = eqdist.etest(z,sizes = N,R = R)$p.value #energy
#    p_values[i,3] = bd.test(x = x,y = y,R = 1000,seed=i*12345)$p.value #Ball
#  }
#  alpha = 0.05
#  pow = colMeans(p_values < alpha)
#  data.frame(methods = c('NN','energy','Ball'),pow)

## -----------------------------------------------------------------------------
set.seed(123)
f = function(x,theta,eta) {
  return(1/(pi * theta * (1+((x-eta)/theta)^2)))
}

n = 10000
x = numeric(n)
u = runif(n)
theta = 1
eta = 0

x[1] = rnorm(1)
k = 0
for(i in 2:n) {
  xt = x[i-1]
  y = rnorm(1,mean=xt,sd=1.5) #choose normal distribution as proposal distribution
  num = f(y,theta,eta) * dnorm(xt,mean=y,sd=1.5)
  den = f(xt,theta,eta) * dnorm(y,mean=xt,sd=1.5)
  if(u[i] <= num/den) {
    x[i] = y
  }else {
    x[i] = xt #y is rejected
    k = k+1
  }
}

index = 1001:n
Estimated = quantile(x[index],seq(0,1,0.1))
Ture = qcauchy(seq(0,1,0.1))
data.frame(Estimated,Ture)

a = ppoints(100)
Qc = qcauchy(a)
Q = quantile(x, a)
qqplot(Qc, Q, xlim = c(-2,2), ylim = c(-2,2), xlab = "Standard Cauchy Quantiles", ylab = "Sample Quantiles")
hist(x[index], breaks = "scott", freq = FALSE)
lines(Qc, f(Qc, 1, 0))

## -----------------------------------------------------------------------------
n = 50
a = 30
b = 40
func = function(x, y) {
  gamma(n + 1) / (gamma(x + 1) * gamma(n - x + 1)) * y^(x + a - 1) * (1 - y)^(n - x + b - 1)
}
# generate the chain with target joint density f(x, y)
m = 10000
d = 2
x = matrix(0, nrow = m, ncol = d)
for (i in 2:m) {
  xt = x[i-1,] #initial (x,y)=(0,0)
  xt[1] = rbinom(1, n, xt[2])
  xt[2] = rbeta(1, xt[1] + a, n - xt[1] + b)
  x[i,] = xt
}
plot(x, cex = 0.5) #show the chain

## ----eval=FALSE---------------------------------------------------------------
#  # For 9.3
#  Gelman.Rubin = function(psi) {
#    psi = as.matrix(psi)
#    n = ncol(psi)
#    k = nrow(psi)
#    psi.means = rowMeans(psi)
#    B = n * var(psi.means)
#    psi.w = apply(psi, 1, "var")
#    W = mean(psi.w)
#    v.hat = W * (n-1)/n + (B/n)
#    r.hat = v.hat / W
#    return(r.hat)
#  }
#  
#  f = function(x,theta,eta) {
#    return(1/(pi * theta * (1+((x-eta)/theta)^2)))
#  }
#  theta = 1
#  eta = 0
#  
#  chain = function(n) {
#    x = numeric(n)
#    u = runif(n)
#    x[1] = rnorm(1)
#    k = 0
#    for(i in 2:n) {
#      xt = x[i-1]
#      y = rnorm(1,mean=xt,sd=1.5) #choose normal distribution as proposal distribution
#      num = f(y,theta,eta) * dnorm(xt,mean=y,sd=1.5)
#      den = f(xt,theta,eta) * dnorm(y,mean=xt,sd=1.5)
#      if(u[i] <= num/den) {
#        x[i] = y
#      }else {
#        x[i] = xt #y is rejected
#        k = k+1
#      }
#    }
#    return(x)
#  }
#  
#  k = 4
#  n = 10000
#  b = 1001
#  
#  #generate the chains
#  set.seed(12)
#  X = matrix(0, nrow=k, ncol=n)
#  for (i in 1:k)
#    X[i, ] = chain(n)
#  
#  #compute diagnostic statistics
#  psi = t(apply(X, 1, cumsum))
#  for (i in 1:nrow(psi))
#    psi[i,] = psi[i,] / (1:ncol(psi))
#  print(Gelman.Rubin(psi))
#  
#  #plot psi for the four chains
#  par(mfrow=c(2,2))
#  for (i in 1:k)
#    plot(psi[i, (b+1):n], type="l", xlab=i, ylab=bquote(psi))
#  par(mfrow=c(1,1)) #restore default
#  
#  #plot the sequence of R-hat statistics
#  rhat = rep(0, n)
#  for (j in (b+1):n)
#    rhat[j] = Gelman.Rubin(psi[,1:j])
#  plot(rhat[(b+1):n], type="l", xlab="",ylim = range(1,1.4), ylab="R")
#  abline(h=1.2, lty=2)

## ----eval=FALSE---------------------------------------------------------------
#  # For 9.8
#  Gelman.Rubin = function(psi) {
#    psi = as.matrix(psi)
#    n = ncol(psi)
#    k = nrow(psi)
#    psi.means = rowMeans(psi)
#    B = n * var(psi.means)
#    psi.w = apply(psi, 1, "var")
#    W = mean(psi.w)
#    v.hat = W * (n-1)/n + (B/n)
#    r.hat = v.hat / W
#    return(r.hat)
#  }
#  
#  chain = function(m) {
#    x = matrix(0, nrow = m, ncol = d)
#    for (i in 2:m) {
#      xt = x[i-1,] #initial (x,y)=(0,0)
#      xt[1] = rbinom(1, n, xt[2])
#      xt[2] = rbeta(1, xt[1] + a, n - xt[1] + b)
#      x[i,] = xt
#    }
#    return(x)
#  }
#  
#  k = 4
#  m = 10000
#  d = 2
#  n = 50
#  a = 30
#  b = 40
#  b1 = 1001
#  
#  #For two-dimensional variables, we monitor the convergence of X and Y respectively
#  #generate the chains
#  set.seed(123)
#  X = matrix(0, nrow=k, ncol=m) #For X
#  Y = matrix(0, nrow=k, ncol=m) #For Y
#  for (i in 1:k){
#    C = chain(m)
#    X[i, ] = C[ ,1]
#    Y[i, ] = C[ ,2]
#  }
#  
#  #compute diagnostic statistics
#  psi = t(apply(X, 1, cumsum))
#  for (i in 1:nrow(psi))
#    psi[i,] = psi[i,] / (1:ncol(psi))
#  print(Gelman.Rubin(psi)) #For X
#  
#  phi = t(apply(Y, 1, cumsum))
#  for (i in 1:nrow(phi))
#    phi[i,] = phi[i,] / (1:ncol(phi))
#  print(Gelman.Rubin(phi)) #For Y
#  
#  #plot psi/phi for the four chains
#  par(mfrow=c(2,2))
#  for (i in 1:k)
#    plot(psi[i, (b1+1):m], type="l", xlab=i, ylab=bquote(psi)) #For X
#  par(mfrow=c(1,1)) #restore default
#  
#  par(mfrow=c(2,2))
#  for (i in 1:k)
#    plot(phi[i, (b1+1):m], type="l", xlab=i, ylab=bquote(phi)) #For Y
#  par(mfrow=c(1,1)) #restore default
#  
#  #plot the sequence of R-hat statistics
#  rhat = rep(0, m)
#  for (j in (b1+1):m)
#    rhat[j] = Gelman.Rubin(psi[,1:j])
#  plot(rhat[(b1+1):m], type="l", xlab="",ylim = range(1,1.4), ylab="R")
#  abline(h=1.2, lty=2) #For X
#  
#  rhat1 = rep(0, m)
#  for (j in (b1+1):m)
#    rhat1[j] = Gelman.Rubin(phi[,1:j])
#  plot(rhat1[(b1+1):m], type="l", xlab="",ylim = range(1,1.4), ylab="R")
#  abline(h=1.2, lty=2) #For Y

## -----------------------------------------------------------------------------
a = c(1,2) #vector
d = length(a) #dimension

#compute the k-th term 
Term = function(a,k) {
  d = length(a)
  return((-1)^k * exp((2*k + 2)*log(norm(a,type = "2")) - lgamma(k + 1) - k*log(2) - log(2*k + 1) - log(2*k + 2) + lgamma((d + 1)/2) + lgamma(k + 3/2) - lgamma(k + d/2 + 1)))
}

#compute the sum
n = 1000
Sum = function(a) {
  sum(sapply(0:n, function(k) Term(a,k)))
}

cat("The sum is",Sum(a))

## -----------------------------------------------------------------------------
solve.equation = function(k) {
  #the integrand
  integrand = function(u,n) {
    (1 + u^2/(n-1))^(-n/2)
  }
  #the upper limit ck
  c_k = function(a,n) {
    sqrt(a^2 * n / (n + 1 - a^2))
  }
  #the expression
  expr = function(a,n) {
    integral = function(u) {
      integrand(u,n)
    }
    c = c_k(a,n-1)
    2/sqrt(pi*(n-1)) * exp(lgamma(n/2)-lgamma((n-1)/2)) * integrate(integral, lower = 0, upper = c)$value
  }
  #express the equation
  f = function(a) {
    left = expr(a,k)
    right = expr(a,k + 1)
    return(left - right)
  }
  #solve the equation
  r = uniroot(f, lower=1, upper=2)$root
  return(r)
}

result = sapply(c(4:25, 100, 500, 1000), function(k) {solve.equation(k)})
result

## -----------------------------------------------------------------------------
m = 10 #number of iterations
n = 10
lambda = numeric(m)
lambda[1] = 2 #initialization

for (i in 2:m) {
  E = 0.54 + 0.48 + 0.33 + 0.43 + 0.91 + 0.21 + 0.85 + 3 * (1 + 1/lambda[i-1])
  lambda[i] = n/E
}

print(1/lambda) #get the estimated expectation
plot(1/lambda)

## -----------------------------------------------------------------------------
trims <- c(0, 0.1, 0.2, 0.5)
x <- rcauchy(100)

lapply(trims, function(trim) mean(x, trim = trim))
lapply(trims, mean, x = x)

## -----------------------------------------------------------------------------
# For Exercise 3
formulas <- list(
        mpg ~ disp,
        mpg ~ I(1/disp),
        mpg ~ disp + wt,
        mpg ~ I(1/disp) + wt
)
# For loops
lo <- list(NULL,NULL,NULL,NULL)
for (i in 1:4) {
  lo[[i]] <- lm(formulas[[i]], data = mtcars)
}
lo
# For lapply()
la <- lapply(formulas, lm, data = mtcars)
la

# extract R^2
rsq <- function(mod) summary(mod)$r.squared
result <- sapply(la, rsq)
cat("提取结果为：",result)

## -----------------------------------------------------------------------------
# For Exercise 5
bootstraps <- lapply(1:10, function(i) {
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})
# For loops
lo2 <- list(NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL)
for (i in 1:10) {
  lo2[[i]] <- lm(mpg ~ disp, data = bootstraps[[i]])
}
lo2
# For lapply() without an anonymous function
la2 <- lapply(bootstraps, lm, formula = mpg ~ disp)
la2

# extract R^2
rsq <- function(mod) summary(mod)$r.squared
result2 <- sapply(la2, rsq)
cat("提取结果为：",result2)

## -----------------------------------------------------------------------------
library(bootstrap)
# (a) a numeric data frame "scor"
result <- vapply(scor, sd, numeric(1))
result
cat("各列标准差为：",result,"\n")
# (b) a mixed data frame "iris"
result2 <- vapply(iris[vapply(iris, is.numeric, logical(1))], sd, numeric(1))
result2
cat("各列标准差为：",result2)

## -----------------------------------------------------------------------------
# sapply与lapply在功能实现上其实是一致的,只是输出结果略有不同,所以利用mclapply去实现mcsapply
mcsapply <- function(X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
  FUN <- match.fun(FUN)
  answer <- parallel::mclapply(X = X, FUN = FUN, ...)
  if (USE.NAMES && is.character(X) && is.null(names(answer)))
      names(answer) <- X
  if (!identical(simplify, FALSE) && length(answer))
      simplify2array(answer, higher = (simplify == "array"))
  else answer
}

## -----------------------------------------------------------------------------
m = 10000
d = 2
# for R
GibbsR = function(n, a, b){
  x = matrix(0, nrow = m, ncol = d)
  for (i in 2:m) {
    xt = x[i-1,]
    xt[1] = rbinom(1, n, xt[2])
    xt[2] = rbeta(1, xt[1] + a, n - xt[1] + b)
    x[i,] = xt
  }
  return(x)
}

# for Rcpp
library(Rcpp)
cppFunction('NumericMatrix GibbsC(int n, int a, int b) {
NumericMatrix x(10000,2);
x(0,0) = 0;
x(0,1) = 0;
for (int i = 1; i < 10000; i++) {
int xt = rbinom(1, n, x(i-1,1))[0];
double yt = rbeta(1, xt + a, n - xt + b)[0];
x(i,0) = xt;
x(i,1) = yt;
};
return x;
}')

## -----------------------------------------------------------------------------
n = 50
a = 30
b = 40
sampleR = GibbsR(n, a, b)
sampleC = GibbsC(n, a, b)
# the variable x of the chain
qqplot(sampleR[,1], sampleC[,1], xlab = "generate by R", ylab = "generate by Rcpp", main = expression("Q-Q plot of the variable x"))
abline(a = 0, b = 1, col = 2, lwd = 2)
# the variable y of the chain
qqplot(sampleR[,2], sampleC[,2], xlab = "generate by R", ylab = "generate by Rcpp", main = expression("Q-Q plot of the variable y"))
abline(a = 0, b = 1, col = 2, lwd = 2)

