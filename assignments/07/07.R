library(psych)
library(tidyverse)

# 9.19
file <- file.choose()
df <- read_tsv(
  file,
  col_names = c(
    "index_growth",
    "index_profit",
    "index_newsales",
    "test_creativity",
    "test_mechanical",
    "test_abstract",
    "test_mathematics"
  )
)

X <- df[ , 1:7]
Z <- t(apply(X, 1, function(x) (x - colMeans(X))/sqrt(diag(cov(X)))))

## a. Factor Analysis by MLE, no rotations
fa.mle.2 <- factanal(Z, 2, rotation = "none")
fa.mle.2
fa.mle.3 <- factanal(Z, 3, rotation = "none")
fa.mle.3

## b. Factor Analysis by MLE, varimax rotations
fa.mle.2 <- factanal(Z, 2, rotation = "varimax", scores = "Bartlett")
fa.mle.2
fa.mle.3 <- factanal(Z, 3, rotation = "varimax", scores = "Bartlett")
fa.mle.3

## c. 
1-fa.mle.2$uniquenesses
### communalities
apply(fa.mle.2$loadings^2,1,sum)
1-fa.mle.3$uniquenesses
fa.mle.2$uniquenesses
### specific variances, or uniqueness \Phi_ii == 1 - h_i^2
fa.mle.3$uniquenesses
### LL' + \Phi
R.2 <- fa.mle.2$loadings %*% t(fa.mle.2$loadings) + diag(fa.mle.2$uniquenesses)
R.3 <- fa.mle.3$loadings %*% t(fa.mle.3$loadings) + diag(fa.mle.3$uniquenesses)
R <- cov(Z)
round(R.2, 3)
round(R.3, 3)
round(R, 3)

## d.
# H_0: \Sigma = \mf{L}\mf{L}' + \Phi, H_a: !=
# ~ \mathcal{X}_{\frac{(p-m)^2-p-m}{2}}
alpha <- 0.01
n <- dim(X)[1]
p <- dim(X)[2]
m <- 2
bartlett_correction <- (n - 1 - (2*p + 4*m +5)/6)
bartlett_correction
bartlett_correction * log(det(R.2)/det(R))
df <- ((p-m)^2-p-m)/2
df
qchisq(1 - alpha, df)
m <- 3
bartlett_correction <- (n - 1 - (2*p + 4*m +5)/6)
bartlett_correction
bartlett_correction * log(det(R.2)/det(R))
df <- ((p-m)^2-p-m)/2
df
qchisq(1 - alpha, df)

## e.
x_new <- c(110, 98, 105, 15, 18, 12, 35) 
Mu <- colMeans(X)
Sigma <- sqrt(diag(cov(X)))
z_new <- (x_new - Mu) / Sigma

m = 2
L <- fa.mle.2$loadings
Phi <- diag(fa.mle.2$uniquenesses)
f_R2 <- t(L) %*% solve(R) %*% z_new
f_LS2 <- (diag(rep(1,m)) + solve(t(L) %*% solve(Phi) %*% L)) %*% f_R2
cbind(f_R2, f_LS2)


m = 3
L <- fa.mle.3$loadings
Phi <- diag(fa.mle.3$uniquenesses)
f_R3 <- t(L) %*% solve(R) %*% z_new
f_LS3 <- (diag(rep(1,m)) + solve(t(L) %*% solve(Phi) %*% L)) %*% f_R3
cbind(f_R3, f_LS3)

# 9.32
file <- file.choose()
df <- read_tsv(file, col_names = FALSE)
colnames(df) = c(
  "Breed",
  "SalePr",
  "YrHgt",
  "FtFrBody",
  "PrctFFB",
  "Frame",
  "BkFat",
  "SaleHt",
  "SaleWt"
)
X <- as.matrix(df[, 3:9])
S <- cov(X)
R <- cor(X)

## a.
## Perform FA on covariance matrix via principal component method
eigen <- eigen(S)
ggplot(data = NULL, mapping = aes(x = 1:7)) + 
  geom_line(aes(y = eigen$values)) +
  geom_point(aes(y = eigen$values), shape = "triangle", size = 3, colour = "red") +
  xlab("index") + 
  ylab("eigenvalue")

fa.pc <- principal(X, covar = TRUE, nfactors = 2, rotate = "none", scores = TRUE)
#fa.pc <- fa(X, nfactors = 2, covar = TRUE, fm = "pa", rotate = "none", score = "Bartlett")
qplot(fa.pc$scores[,1], fa.pc$scores[,2]) + xlab("factor_1") + ylab("factor_2")
fa.pc

fa.pc.rot <- principal(X, covar = TRUE, nfactors = 2, rotate = "varimax", scores = TRUE)
#fa.pc.rot <- fa(X, nfactors = 2, covar = TRUE, fm = "pa", rotate = "varimax", score = "Bartlett")
qplot(fa.pc.rot$scores[,1], fa.pc.rot$scores[,2]) + xlab("factor_1") + ylab("factor_2")
fa.pc.rot

fa.mle <- fa(X, nfactors = 2, covar = TRUE, fm = "ml", rotate = "varimax", score = "Bartlett")
qplot(fa.mle$scores[,1], fa.mle$scores[,2]) + xlab("factor_1") + ylab("factor_2")
fa.mle

## b.
eigen <- eigen(R)
ggplot(data = NULL, mapping = aes(x = 1:7)) + 
  geom_line(aes(y = eigen$values)) +
  geom_point(aes(y = eigen$values), shape = "triangle", size = 3, colour = "red") +
  xlab("index") + 
  ylab("eigenvalue")

fa.pc <- principal(X, covar = FALSE, nfactors = 3, rotate = "none", scores = TRUE)
qplot(fa.pc$scores[,1], fa.pc$scores[,2]) + xlab("factor_1") + ylab("factor_2")
fa.pc

fa.pc.rot <- principal(X, covar = FALSE, nfactors = 3, rotate = "varimax", scores = TRUE)
qplot(fa.pc.rot$scores[,1], fa.pc.rot$scores[,2]) + xlab("factor_1") + ylab("factor_2")
fa.pc.rot

fa.mle <- fa(X, nfactors = 2, covar = FALSE, fm = "ml", rotate = "varimax", score = "Bartlett")
qplot(fa.mle$scores[,1], fa.mle$scores[,2]) + xlab("factor_1") + ylab("factor_2")
fa.mle
