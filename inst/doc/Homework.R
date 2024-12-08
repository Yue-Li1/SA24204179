## ----echo=FALSE---------------------------------------------------------------
library(knitr)
data <- data.frame(
  Name=c("A","B","C","D","E"),
  Content=c(0.11,0.13,0.14,0.16,0.18),
  Strength=c(43,45,45,49,50)
)
kable(data,caption="A Table")

## ----echo=FALSE,fig.width=8, fig.height=8-------------------------------------
x <- c(18,20,22,24,26)
y <- c(26.86,28.35,28.75,28.87,29.75)
lm01 <- lm(x~y)
summary(lm01)$coef
par(mfrow=c(2,2))
plot(lm01)

## ----echo=FALSE,fig.width=8, fig.height=8-------------------------------------
n <- 10000
my_function <- function(n,sigma){
  u <- runif(n)
  return(sigma * sqrt(-2*log(u)))
}
par(mfrow=c(2,2))
for(sigma in 1:4){
  samples <- my_function(n,sigma)
  hist(samples,prob = TRUE,main=paste("Rayleigh Distribution(σ=",sigma,")"),
       xlab="x",ylab="Density")
y <- seq(0,max(samples),length=1000)
rayleigh_density <- (y / sigma^2) * exp(-y^2 / (2 * sigma^2))
lines(y,rayleigh_density)
}

## ----echo=FALSE---------------------------------------------------------------
set.seed(123)
n <- 1000
p1 <- 0.75
b <- rbinom(n,size = 1,prob = p1)
sample_data <- b*rnorm(n,mean=0,sd=1)+(1-b)*rnorm(n,mean=3,sd=1)
hist(sample_data,prob=TRUE,main=expression("混合分布的直方图(P1=0.75)"),
     xlab="x",ylab="Density")
y <- seq(min(sample_data),max(sample_data),length=1000)
lines(y,p1*dnorm(y,mean=0,sd=1)+(1-p1)*dnorm(y,mean=3,sd=1))

## ----echo=FALSE---------------------------------------------------------------
set.seed(1234)
n <- 1000
p1_data <- c(0.3,0.45,0.5,0.55,0.7)
par(mfrow=c(2,3))
for(p1 in p1_data){
  b <- rbinom(n,size = 1,prob = p1)
  sample_data <- b*rnorm(n,mean=0,sd=1)+(1-b)*rnorm(n,mean=3,sd=1)
  hist(sample_data,prob=TRUE,main=paste("混合分布的直方图(P1=",p1,")"),
     xlab="x",ylab="Density")
y <- seq(min(sample_data),max(sample_data),length=1000)
lines(y,p1*dnorm(y,mean=0,sd=1)+(1-p1)*dnorm(y,mean=3,sd=1))
}

## ----echo=FALSE---------------------------------------------------------------
set.seed(12345)
lambda <- 3
alpha <- 4
beta <- 2
t <- 10
n <- 10000
compound_function <- function(lambda,alpha,beta,t,n){
  result <- numeric(n)
  for(k in 1:n){
    N_t <- rpois(1,lambda*t)
    if(N_t > 0){
      Y_k <- rgamma(N_t,shape=alpha,rate=beta)
      result[k] <- sum(Y_k)
  }else{
    result[k] <- 0
  }
  }
  return(result)
}
X_t <- compound_function(lambda,alpha,beta,t,n)
N_t_values <- compound_function(lambda,alpha,beta,t,n)
estimated_mean <- mean(X_t)
estimated_variance <- var(X_t)
theoretical_mean <- lambda*t*alpha/beta
theoretical_variance <- lambda*t*(alpha * (1 + alpha) / beta^2)
cat("模拟的X(10)均值：",estimated_mean,"\n")
cat("模拟的X(10)方差：",estimated_variance,"\n")
cat("理论的X(10)均值：",theoretical_mean,"\n")
cat("理论的X(10)方差：",theoretical_variance,"\n")

## -----------------------------------------------------------------------------
m <- 100000
x_values <- seq(0.1,0.9,by=0.1)
mc_cdf <- numeric(length(x_values))
for(i in 1:length(x_values)){
  x <- runif(m,min=0,max=x_values[i])
  mc_cdf[i] <- mean(30*x^2*(1-x)^2)*x_values[i]
}
true_value <- pbeta(x_values,shape1=3,shape2=3)
result <- data.frame(x=x_values,MC_estimate=mc_cdf,True_value=true_value)
print(result)

## -----------------------------------------------------------------------------
n <- 10000
sigma <- 2
u <- runif(n)
x1 <- sigma*sqrt(-2*log(u))
x2 <- sigma*sqrt(-2*log(1-u))
y <- (x1+x2)/2
var_x1 <- var(x1)
var_x2 <- var(x2)
var_independent <- (var_x1+var_x2)/2
var_dependent <- var(y)
var_reduction <- (var_independent-var_dependent)/var_independent*100
print(paste("方差减少百分比:",round(var_reduction,2),"%"))

## -----------------------------------------------------------------------------
m <- 10000
theta.hat <- se <- numeric(2)
g <- function(x) {
  (x^2 / sqrt(2 * pi)) * exp(-x^2 / 2) * (x > 1)
}
#using f1
x <- rexp(m,rate=1/2)
fg <- g(x) / dexp(x,rate=1/2)
theta.hat[1] <- mean(fg)
se[1] <- sd(fg)
#using f2
y <- rnorm(m,mean = 3,sd=1)
fg <- g(y) / dnorm(y,mean=1,sd=1)
theta.hat[2] <- mean(fg)
se[2] <- sd(fg)
cat("使用 f1 的估计值:", theta.hat[1], "，标准误:", se[1], "\n")
cat("使用 f2 的估计值:", theta.hat[2], "，标准误:", se[2], "\n")

## -----------------------------------------------------------------------------
fast_sort <- function(x) {
  if(length(x) <= 1) {
    return(x)
  }
  p <- sample(x,1)
  left <- x[x < p]
  middle <- x[x == p]
  right <- x[x > p]
  return(c(fast_sort(left), middle, fast_sort(right)))
}
n <- c(1e4,2e4,4e4,6e4,8e4)
an <- 100
sort_time <- function(n,an){
  times <- numeric(an)
  for(i in 1:an){
    randomly_numbers <- sample(1:n,n)
    begin <- Sys.time()
    fast_sort(randomly_numbers)
    end <- Sys.time()
    times[i] <- as.numeric(difftime(end,begin,units = "secs"))
  }
  return(mean(times))
}
sort_times <- sapply(n,sort_time,an=an)
tn <- n * log(n)
fit <- lm(sort_times ~ tn)
plot(tn, sort_times,main="Scatter Plot and Regression Line",xlab = "n * log(n)",ylab="Calculate computation time")
abline(fit)

## -----------------------------------------------------------------------------
n <- 1000
m <- 1000
skewness_m <- numeric(m)
skewness <- function(x) {
  n <- length(x)
  mean_x <- mean(x)
  m2 <- mean((x - mean_x)^2)
  m3 <- mean((x - mean_x)^3)
  sqrt_b1 <- m3 / m2^(3/2)
  return(sqrt_b1)
}
for(i in 1:m){
  samples <- rnorm(n)
  skewness_m[i] <- skewness(samples)
}
q_estimate <- quantile(skewness_m,c(0.025,0.05,0.95,0.975))
se_estimate <- sd(skewness_m)/sqrt(m)
q_approximate <- qnorm(c(0.025,0.05,0.95,0.975),mean=0,sd=sqrt(6/n))
print(paste("Standard Error: ",se_estimate))
comparsion <- data.frame(Quantiles=c(0.025,0.05,0.95,0.975),
      MC_Experiment_Estimate=q_estimate,
      Large_Sample_Approximate=q_approximate)
print(comparsion)

## -----------------------------------------------------------------------------
set.seed(123)
n <- 100
m <- 1000
alpha <- 0.05
power_tests1 <- function(n,m,alpha){
  pearson_test <- spearman_test <- kendall_test <- numeric(m)
  for(i in 1:m){
    x <- rnorm(n)
    y <- 0.5*x + rnorm(n)
    pearson_test[i] <- cor.test(x, y, method = "pearson")$p.value
    spearman_test[i] <- cor.test(x, y, method = "spearman")$p.value
    kendall_test[i] <- cor.test(x, y, method = "kendall")$p.value
  }
  power_pearson <- mean(pearson_test < alpha)
  power_spearman <- mean(spearman_test < alpha)
  power_kendall <- mean(kendall_test < alpha)
  print(paste("Pearson Power: ", power_pearson))
  print(paste("Spearman Power: ", power_spearman))
  print(paste("Kendall Power: ", power_kendall))
}
power_tests1(n, m, alpha)

## -----------------------------------------------------------------------------
power_tests2 <- function(n,m,alpha){
  pearson_test <- spearman_test <- kendall_test <- numeric(m)
  for(i in 1:m){
    x <- rnorm(n)
    y <- sin(x) + x^2 + rnorm(n)
    pearson_test[i] <- cor.test(x, y, method = "pearson")$p.value
    spearman_test[i] <- cor.test(x, y, method = "spearman")$p.value
    kendall_test[i] <- cor.test(x, y, method = "kendall")$p.value
  }
  power_pearson <- mean(pearson_test < alpha)
  power_spearman <- mean(spearman_test < alpha)
  power_kendall <- mean(kendall_test < alpha)
  print(paste("Pearson Power: ", power_pearson))
  print(paste("Spearman Power: ", power_spearman))
  print(paste("Kendall Power: ", power_kendall))
}
power_tests2(n, m, alpha)

## -----------------------------------------------------------------------------
p1 <- 0.651
p2 <- 0.676
n1 <- 10000
n2 <- 10000
p_total <- (p1*n1+p2*n2)/(n1+n2)
Z_value <- (p1-p2)/(sqrt((p_total)*(1-p_total)*(1/n1+1/n2)))
P_value <- 2*(1-pnorm(abs(Z_value)))
if(P_value < 0.05){
  print(paste("p值：",P_value))
  print("拒绝原假设，两个方法的功效存在显著差异")
}else{
  print(paste("p值：",P_value))
  print("接受原假设，两个方法的功效没有显著差异")
}

## -----------------------------------------------------------------------------
N <- 1000
m <- 10000
alpha <- 0.1
results_table <- matrix(0, nrow=3, ncol=2)
colnames(results_table) <- c("Bonferroni correction", "B-H correction")
rownames(results_table) <- c("FWER", "FDR", "TPR")

for (i in 1:m) {
  pvalues <- c(runif(950),rbeta(50, 0.1, 1))
  adjusted_bonferroni <- p.adjust(pvalues, method = "bonferroni") < alpha
  adjusted_bh <- p.adjust(pvalues, method = "BH") < alpha
  results_table["FWER", "Bonferroni correction"] <- results_table["FWER", "Bonferroni correction"] + any(adjusted_bonferroni[1:950])
  results_table["FDR", "Bonferroni correction"] <- results_table["FDR", "Bonferroni correction"] + sum(adjusted_bonferroni[1:950]) / max(1, sum(adjusted_bonferroni))
  results_table["TPR", "Bonferroni correction"] <- results_table["TPR", "Bonferroni correction"] + sum(adjusted_bonferroni[951:1000]) / 50
  
  results_table["FWER", "B-H correction"] <- results_table["FWER", "B-H correction"] + any(adjusted_bh[1:950])
  results_table["FDR", "B-H correction"] <- results_table["FDR", "B-H correction"] + sum(adjusted_bh[1:950]) / max(1, sum(adjusted_bh))
  results_table["TPR", "B-H correction"] <- results_table["TPR", "B-H correction"] + sum(adjusted_bh[951:1000]) / 50
}
results_table <- results_table/m
print(results_table)

## -----------------------------------------------------------------------------
data <- c(3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487)
lambda_hat <- 1/mean(data)
B <- 1000
lambda_star <- numeric(B)
set.seed(1234)
for(b in 1:B) {
  x_star <- sample(data, size = length(data), replace = TRUE)
  lambda_star[b] <- 1/mean(x_star)
}
bias <- mean(lambda_star) - lambda_hat
standard_error <- sd(lambda_star)
cat("λ的最大似然估计为:", lambda_hat, "\n")
cat("自举法估计的偏差为:", bias, "\n")
cat("自举法估计的标准误差为:", standard_error, "\n")

## -----------------------------------------------------------------------------
library(boot)
data <- c(3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487)
function_star <- function(data,i){
  lambda_hat <- 1 / mean(data[i])
  return(1 / lambda_hat)
}
set.seed(123)
results <- boot(data,statistic = function_star, R=1000)
ci <- boot.ci(results, type = c("norm", "basic", "perc", "bca"))
print(ci)

## -----------------------------------------------------------------------------
library(bootstrap)
data1 <- scor
n <- nrow(data1)
Cov_matrix <- cov(data1)
Positive_eigenvalues <- eigen(Cov_matrix)$values
theta_hat <- Positive_eigenvalues[1] / sum(Positive_eigenvalues)

jack_estimates <- numeric(n)
for (i in 1:n) {
  jack_data <- data1[-i, ]
  jack_Cov <- cov(jack_data)
  jack_eigenvalues <- eigen(jack_Cov)$values
  jack_estimates[i] <- jack_eigenvalues[1] / sum(jack_eigenvalues)
}

jack_bias <- (n-1) * (mean(jack_estimates) - theta_hat)
jack_se <- sqrt((n-1) * mean((jack_estimates - mean(jack_estimates))^2))

cat("估计值 theta_hat:", theta_hat, "\n")
cat("Jackknife 偏差估计:", jack_bias, "\n")
cat("Jackknife 标准误估计:", jack_se, "\n")

## -----------------------------------------------------------------------------
library(DAAG)
data(ironslag)
magnetic <- ironslag$magnetic
chemical <- ironslag$chemical
n <- length(magnetic)
e1 <- e2 <- e3 <- e4 <- numeric(n)

for (k in 1:n) {
  y <- magnetic[-k]
  x <- chemical[-k]
  
  J1 <- lm(y~x)
  yhat1 <- J1$coef[1] + J1$coef[2] * chemical[k]
  e1[k] <- magnetic[k]- yhat1
  
  J2 <- lm(y ~ x + I(x^2))
  yhat2 <- J2$coef[1] + J2$coef[2] * chemical[k] + J2$coef[3] * chemical[k]^2
  e2[k] <- magnetic[k]- yhat2
  
  J3 <- lm(log(y) ~ x)
  logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[k]
  yhat3 <- exp(logyhat3)
  e3[k] <- magnetic[k]- yhat3
  
  J4 <- lm(y ~ x + I(x^2) + I(x^3))
  yhat4 <- J4$coef[1] + J4$coef[2] * chemical[k] + J4$coef[3] * chemical[k]^2 + J4$coef[4] * chemical[k]^3
  e4[k] <- magnetic[k] - yhat4
}

rmse1 <- sqrt(mean(e1^2))
rmse2 <- sqrt(mean(e2^2))
rmse3 <- sqrt(mean(e3^2))
rmse4 <- sqrt(mean(e4^2))
Adjusted_R2_1 <- summary(J1)$adj.r.squared
Adjusted_R2_2 <- summary(J2)$adj.r.squared
Adjusted_R2_3 <- summary(J3)$adj.r.squared
Adjusted_R2_4 <- summary(J4)$adj.r.squared

cat("模型1的RMSE和调整后的R²:", rmse1, Adjusted_R2_1, "\n")
cat("模型2的RMSE和调整后的R²:", rmse2, Adjusted_R2_2, "\n")
cat("模型3的RMSE和调整后的R²:", rmse3, Adjusted_R2_3, "\n")
cat("模型4的RMSE和调整后的R²:", rmse4, Adjusted_R2_4, "\n")

## -----------------------------------------------------------------------------
library(DAAG)
data(chickwts)
x <- sort(as.vector(chickwts$weight[chickwts$feed == "soybean"]))
y <- sort(as.vector(chickwts$weight[chickwts$feed == "linseed"]))
Cramer_von_Mises_test <- function(x,y){
  n <- length(x)
  m <- length(y)
  z <- c(x,y)
  ECDF_x <- ecdf(x)(z)
  ECDF_y <- ecdf(y)(z)
  Cramer_von_Mises_statistic <- sum((ECDF_x-ECDF_y)^2) * (n*m) / (m+n)^2
  return(Cramer_von_Mises_statistic)
}
W0 <- Cramer_von_Mises_test(x,y)

set.seed(123)
R <- 1000
W1 <- numeric(R)
for (k in 1:R) {
  z <- sample(c(x, y))
  W1[k] <- Cramer_von_Mises_test(z[1:length(x)],z[(length(x)+1):length(z)])
}
p_value <- mean(W1 >= W0)

cat("原始 Cramér-von Mises 统计量:", W0, "\n")
cat("p 值:", p_value, "\n")

## -----------------------------------------------------------------------------
set.seed(1234)
n <- 50
x <- rnorm(n)
y <- rnorm(n)
Spearman0 <- cor(x,y,method = "spearman")

R <- 1000
Spearman1 <- numeric(R)
for (i in 1:R){
  y_permute <- sample(y)
  Spearman1[i] <- cor(x,y_permute,method = "spearman")
}
p_value_permutation <- mean(abs(Spearman1) >= abs(Spearman0))
cor_test <- cor.test(x,y,method = "spearman")
p_value_cor_test <- cor_test$p.value

cat("置换检验的p值:", p_value_permutation, "\n")
cat("cor.test的p值:", p_value_cor_test, "\n")

## -----------------------------------------------------------------------------
set.seed(1234)
Metropolis_hastings <- function(m, sigma) {
  x <- numeric(m)
  x[1] <- 0
  for(i in 2:m) {
    #从正态分布生成候选样本
    y <- rnorm(1, mean = x[i-1], sd = sigma)
    f_y <- 1 / (pi * (1+y^2))
    f_x <- 1 / (pi * (1+x[i-1]^2))
    r <- f_y / f_x
    u <- runif(1)
    if (u < r) {
      x[i] <- y
    } else {
      x[i] <- x[i-1]
    }
  }
  return(x)
}
n <- 10000
sigma <- 1
samples <- Metropolis_hastings(n,sigma)
final_samples <- samples[-(1:1000)]
deciles_generated <- quantile(final_samples,probs = seq(0.1,0.9,by=0.1))
deciles_theoretical <- qcauchy(seq(0.1,0.9,by=0.1))

cat("生成样本的十分位数:\n")
print(deciles_generated)
cat("\n理论标准Cauchy分布的十分位数:\n")
print(deciles_theoretical)

## -----------------------------------------------------------------------------
set.seed(1234)
n <- 10
a <- 2
b <- 2
m <- 10000
Gibbs_sampler <- function(n,a,b,m) {
  X_Y <- matrix(NA,nrow = m,ncol = 2)
  x <- 0
  y <- runif(1)
  for(i in 1:m){
    x <- rbinom(1,n,y)
    y <- rbeta(1,x+a,n-x+b)
    X_Y[i, ] <- c(x,y)
  }
  return(X_Y)
}
Gibbs_samples <- Gibbs_sampler(n,a,b,m)
head(Gibbs_samples, 20)

## -----------------------------------------------------------------------------
set.seed(1234)
Metropolis_hastings <- function(m, sigma) {
  x <- numeric(m)
  x[1] <- 0
  for(i in 2:m) {
    y <- rnorm(1, mean = x[i-1], sd = sigma)
    f_y <- 1 / (pi * (1+y^2))
    f_x <- 1 / (pi * (1+x[i-1]^2))
    r <- f_y / f_x
    u <- runif(1)
    if (u < r) {
      x[i] <- y
    } else {
      x[i] <- x[i-1]
    }
  }
  return(x)
}
n <- 10000
sigma <- 1
chains_num <- 4
Gelman_rubin <- function(psi1){
  m <- length(psi1)
  n <- length(psi1[[1]])
  means_chain <- sapply(psi1,mean)
  B <- n * var(means_chain)
  W <- mean(sapply(psi1,var))
  R_hat <- ((B / W)/n + (n - 1) / n)
  return(R_hat)
}
psi1 <- lapply(1:chains_num, function(i)Metropolis_hastings(n, sigma))
R_hat1 <- Gelman_rubin(psi1)
cat("Gelman-Rubin R-hat:\n")
print(R_hat1)

if (all(R_hat1 < 1.2)) {
  cat("链已收敛到目标分布。\n")
} else {
  cat("链尚未收敛，请继续运行。\n")
}

## -----------------------------------------------------------------------------
set.seed(1234)
n <- 10
a <- 2
b <- 2
chains_num <- 4
m <- 10000
Gibbs_sampler <- function(n,a,b,m) {
  X_Y <- matrix(NA,nrow = m,ncol = 2)
  x <- 0
  y <- runif(1)
  for(i in 1:m){
    x <- rbinom(1,n,y)
    y <- rbeta(1,x+a,n-x+b)
    X_Y[i, ] <- c(x,y)
  }
  return(X_Y)
}
Gelman_rubin <- function(psi2){
  m <- length(psi2)
  n <- nrow(psi2[[1]])
  means_chain <- sapply(psi2,colMeans)
  B <- n * apply(means_chain,1,var)
  W <- mean(sapply(psi2,function(chain)apply(chain,2,var)))
  R_hat <- ((B / W)/n + (n - 1) / n)
  return(R_hat)
}
psi2 <- lapply(1:chains_num, function(i) Gibbs_sampler(n, a, b, m))

R_hat2 <- Gelman_rubin(psi2)
cat("Gelman-Rubin R-hat:\n")
print(R_hat2)
if (all(R_hat2 < 1.2)) {
  cat("链已收敛到目标分布。\n")
} else {
  cat("链尚未收敛，请继续运行。\n")
}

## -----------------------------------------------------------------------------
#(a)计算第k项
kth_term <- function(a,k,d){
  a_norm <- sqrt(sum(a^2))
  gamma_term <- gamma((d+1)/2) * gamma(k+3/2) / gamma(k+1+d/2)
  k_term <- a_norm^(2*k+2)*(-1)^k / (factorial(k)*2^k*(2*k+1)*(2*k+2)) *gamma_term
  return(k_term)
}
#(b)求和
sum_kth <- function(a,d,e=1e-10,k_max=1000){
  sum <- 0
  k <- 0
  repeat{
    new_term <- kth_term(a,k,d)
    sum <- sum + new_term
    if(abs(new_term) < e || k > k_max){
      break
    }
    k <- k+1
  }
  return(sum)
}
#(c)求和(向量a已知)
a <- c(1, 2)
d <- 2
result <- sum_kth(a,d)
cat("无穷和为:", result, "\n")

## -----------------------------------------------------------------------------
ck_function <- function(a,k){
  ck <- sqrt((a^2 *k)/(k+1-a^2))
  return(ck)
}
LHS <- function(a,k){
  c_k1 <- ck_function(a,k-1)
  integral1 <- integrate(function(u)(1+u^2/(k-1))^(-k/2),lower=0,upper=c_k1)$value
  LHS_value <- 2*gamma(k/2)*integral1/(sqrt(pi*(k-1))*gamma((k-1)/2))
  return(LHS_value)
}
RHS <- function(a,k){
  c_k2 <- ck_function(a,k)
  integral2 <- integrate(function(u)(1+u^2/k)^(-(k+1)/2),lower=0,upper=c_k2)$value
  RHS_value <- 2*gamma((k+1)/2)*integral2/(sqrt(pi*k)*gamma(k/2))
  return(RHS_value)
}

a_result1 <- function(k){
  equation1 <- function(a){
    LHS(a,k) - RHS(a,k)
  }
  solution1 <- uniroot(equation1,lower=0.1,upper=sqrt(k))$root
  return(solution1)
}

k <- 4
a_solution1 <- a_result1(k)
cat("当 k =", k, "时，求解的 a 值为:", a_solution1, "\n")

## -----------------------------------------------------------------------------
#11.4的结果
Sk_minus_1 <- function(a, k) {
  t_minus_1 <- sqrt(a^2 * (k - 1) / (k - a^2))
  1 - pt(t_minus_1, df = k - 1)
}
Sk <- function(a, k) {
  t <- sqrt(a^2 * k / (k + 1 - a^2))
  1 - pt(t, df = k)
}
a_result2 <- function(k){
  equation2 <- function(a){
    Sk_minus_1(a,k) - Sk(a,k)
  }
  solution2 <- uniroot(equation2,lower=0.001,upper=sqrt(k))$root
  return(solution2)
}
k <- 4
a_solution2 <- a_result2(k)
cat("当 k =", k, "时，求解的 a 值为:", a_solution2, "\n")

## -----------------------------------------------------------------------------
Y <- c(0.54, 0.48, 0.33, 0.43, 1.00, 1.00, 0.91, 1.00, 0.21, 0.85)
tau <- 1
lambda <- 0.6
e <- 1e-6
i_max <- 1000
for(i in 1:i_max){
  lambda0 <- lambda
  observed_Y <- Y[Y < tau]
  miss_Y <- Y[Y == tau]
  lambda <- (sum(observed_Y)+length(miss_Y)*(lambda+tau))/length(Y)
  if(abs(lambda-lambda0)<e){
    break
  }
}
lambda_mle <- mean(Y)
cat("EM算法估计的λ:",lambda,"\n")
cat("MLE估计的λ:",lambda_mle,"\n")

## -----------------------------------------------------------------------------
library(boot)
A1 <- rbind(c(2,1,1),c(1,-1,3))
b1 <- c(2,3)
a <- c(4,2,9)
result <- simplex(a=a,A1=A1,b1=b1,maxi=FALSE)
print(result)

## -----------------------------------------------------------------------------
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1/disp),
  mpg ~ disp + wt,
  mpg ~ I(1/disp) + wt
)
result_loop <- list()
for (i in 1:4) {
  result_loop[[i]] <- lm(formulas[[i]], data = mtcars)
}

## -----------------------------------------------------------------------------
formulas <- list(
  mpg ~ disp,
  mpg ~ I(1/disp),
  mpg ~ disp + wt,
  mpg ~ I(1/disp) + wt
)
result_lapply <- lapply(formulas, function(f)lm(f, data = mtcars))

## -----------------------------------------------------------------------------
bootstraps <- lapply(1:10,function(i){
  rows <- sample(1:nrow(mtcars), rep = TRUE)
  mtcars[rows, ]
})
results <- list()
for(i in 1:10){
  results[[i]] <- lm(mpg ~ disp,data=bootstraps[[i]])
}

## -----------------------------------------------------------------------------
rsq <- function(mod) summary(mod)$r.squared
R21 <- sapply(result_loop, rsq)
cat("练习3模型的R^2:\n")
print(R21)
R22 <- sapply(results, rsq)
cat("练习4模型的R^2:\n")
print(R22)

## -----------------------------------------------------------------------------
trials <- replicate(
  100,
  t.test(rpois(10, 10), rpois(7, 10)),
  simplify = FALSE
)
p_values <- sapply(trials, function(x) x$p.value)
print(p_values[1:5])

## -----------------------------------------------------------------------------
MVlapply <- function(..., f, simplify = TRUE) {
  list_result <- Map(f, ...)
  if (simplify) {
    vector_result <- vapply(list_result, FUN.VALUE = numeric(length(list_result[[1]])), FUN = identity)
    return(vector_result)
  } else {
    return(list_result)
  }
}

## -----------------------------------------------------------------------------
fast_chisq_test <- function(observed1, observed2) {
  if (length(x) != length(y)) stop("两个样本的长度必须一样")
  X_Y <- table(observed1, observed2)
  expected <- outer(rowSums(X_Y), colSums(X_Y), FUN = "*") / sum(X_Y)
  chisq_stat <- sum((X_Y - expected)^2 / expected)
  return(chisq_stat)
}

## -----------------------------------------------------------------------------
fast_table <- function(observed1, observed2) {
  observed3 <- unique(c(observed1, observed2))
  table_new <- matrix(0,nrow=length(observed3),ncol=length(observed3),dimnames = list(observed3, observed3))
  for (i in seq_along(observed1)) {
    table_new[observed1[i], observed2[i]] <- table_new[observed1[i], observed2[i]] + 1
  }
  return(table_new)
}

## -----------------------------------------------------------------------------
library(Rcpp)
cppFunction('
#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::export]]
NumericMatrix Gibbs(int n, double a, double b, int m){
  NumericMatrix X_Y(m,2);
  double x = 0;
  double y = R::runif(0,1);
  for(int i = 0; i < m; i++){
    x = R::rbinom(n,y);
    y = R::rbeta(x+a,n-x+b);
    X_Y(i,0) = x;
    X_Y(i,1) = y;
  }
  return(X_Y);
}')

n <- 10
a <- 2
b <- 2
m <- 1000

samples_gibbs <- Gibbs(n, a, b, m)
x_gibbs <- samples_gibbs[, 1]
y_gibbs <- samples_gibbs[, 2]

Corresponding_samples <- function(n, a, b, m) {
  x_random <- rbinom(m, n, 0.5)
  y_random <- rbeta(m, x_random+a, n-x_random+b)
  Random_samples <- cbind(x_random, y_random)
  return(Random_samples)
}
Random_samples <- Corresponding_samples(n, a, b, m)
x_random <- Random_samples[, 1]
y_random <- Random_samples[, 2]

par(mfrow = c(1, 2))
qqplot(x_gibbs, x_random, xlab = "Gibbs sampler x", ylab = "Random x")
qqplot(y_gibbs, y_random, xlab = "Gibbs sampler y", ylab = "Random y")

## -----------------------------------------------------------------------------
library(microbenchmark)
ts <- microbenchmark(
  Gibbs = Gibbs(n, a, b, m),
  Random = Corresponding_samples(n, a, b, m),
  times = 10
)
summary(ts)[,c(1,3,5,6)]

