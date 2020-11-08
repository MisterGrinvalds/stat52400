library(tidyverse)

# 8.20
## Use T8-6.dat :: Men's
file <- file.choose()
df <- read_tsv(file, col_names = FALSE)

# a.
R <- cor(df[,2:9])
R
eigen(R)

# b.
X <- as.matrix(df[,2:9])
Mu <- colMeans(X)
Sigma <- sqrt(diag(cov(X)))
Z <- t(apply(X, 1, function(x) (x - Mu)/Sigma))
Rz <- cor(Z)
eigen <- eigen(Rz)
eigen
table <- cbind(
  t(t(eigen$vectors[,1] * sqrt(eigen$values[1]))),
  t(t(eigen$vectors[,2] * sqrt(eigen$values[2]))),
  eigen$values,
  cumsum(eigen$values) / sum(eigen$values)
)
colnames(table) <- c("rho_Y1X", "rho_Y2X", "eigenvalue", "proportion")
table

# d.
Zdf <- data.frame(pc1 = Z %*% (eigen$vectors[,1]))
sort_i <- order(Zdf, decreasing = TRUE)
cbind(df$X1[sort_i], Zdf[sort_i,])

# e.
## Use T1-9.dat :: Women's
file <- file.choose()
df <- read_tsv(file, col_names = FALSE)

X <- as.matrix(df[,2:8])
Mu <- colMeans(X)
Sigma <- sqrt(diag(cov(X)))
Z <- t(apply(X, 1, function(x) (x - Mu)/Sigma))
Rz <- cor(Z)
eigen <- eigen(Rz)
table <- cbind(
  t(t(eigen$vectors[,1] * sqrt(eigen$values[1]))),
  t(t(eigen$vectors[,2] * sqrt(eigen$values[2]))),
  eigen$values,
  cumsum(eigen$values) / sum(eigen$values)
)
colnames(table) <- c("rho_Y1X", "rho_Y2X", "eigenvalue", "proportion")
table
Zdf <- data.frame(pc1 = Z %*% (eigen$vectors[,1]))
sort_i <- order(Zdf, decreasing = TRUE)
cbind(df$X1[sort_i], Zdf[sort_i,])

# 8.21
## Use T8-6.dat :: Men's
file <- file.choose()
df <- read_tsv(file, col_names = FALSE)
df.speeds <- df %>% mutate(
  X2 = 100 / X2,
  X3 = 200 / X3,
  X4 = 400 / X4,
  X5 = 800 / X5 / 60,
  X6 = 1500 / X6 / 60,
  X7 = 5000 / X7 / 60,
  X8 = 10000 / X8 / 60,
  X9 = 42195 / X9 / 60
)
colnames(df.speeds) <- c(
  "nation",
  "100m",
  "200m",
  "400m",
  "800m",
  "1.5km",
  "5km",
  "10k",
  "marathon"
)
df.speeds

X <- as.matrix(df.speeds[,2:9])
S <- cov(X)
eigen <- eigen(S)
eigen
table <- cbind(
  t(t(eigen$vectors[,1] * sqrt(eigen$values[1]))),
  t(t(eigen$vectors[,2] * sqrt(eigen$values[2]))),
  eigen$values,
  cumsum(eigen$values) / sum(eigen$values)
)
colnames(table) <- c("rho_Y1X", "rho_Y2X", "eigenvalue", "proportion")
table

Sdf <- data.frame(pc1 = t(apply(X, 1, function(x) (x - colMeans(X)))) %*% (eigen$vectors[,1]))
sort_i <- order(Sdf, decreasing = FALSE)
cbind(df$X1[sort_i], Sdf[sort_i,])

# 8.22
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

# PCA on S
eigen <- eigen(S)
eigen
cbind(eigenvalue = eigen$values, total_proportion = cumsum(eigen$values) / sum(eigen$values))
cbind(colnames(df)[3:9], eigen$vectors[,2:3])
ggplot(data = NULL, mapping = aes(x = 1:7)) + 
  geom_line(aes(y = eigen$values)) +
  geom_point(aes(y = eigen$values), shape = "triangle", size = 3, colour = "red") +
  xlab("index") + 
  ylab("eigenvalue")
Y = X %*% eigen$vectors[,1:2]
ggplot(data = NULL, mapping = aes(y = Y[,1], x = Y[,2], colour = as.factor(df$Breed))) + 
  geom_point() +
  xlab("y2_hat") +
  ylab("y1_hat") +
  labs(colour = "breed")
ggplot(data = NULL, mapping = aes(sample = Y[,1])) +
  geom_qq() +
  geom_qq_line()

# PCA on R
eigen <- eigen(R)
eigen
cbind(eigenvalue = eigen$values, total_proportion = cumsum(eigen$values) / sum(eigen$values))
cbind(colnames(df)[3:9], eigen$vectors[,2:4])
ggplot(data = NULL, mapping = aes(x = 1:7)) + 
  geom_line(aes(y = eigen$values)) +
  geom_point(aes(y = eigen$values), shape = "triangle", size = 3, colour = "red") +
  xlab("index") + 
  ylab("eigenvalue")
Y = X %*% eigen$vectors[,1:2]
ggplot(data = NULL, mapping = aes(y = Y[,1], x = Y[,2], colour = as.factor(df$Breed))) + 
  geom_point() +
  xlab("y2_hat") +
  ylab("y1_hat") +
  labs(colour = "breed")
ggplot(data = NULL, mapping = aes(sample = Y[,1])) +
  geom_qq() +
  geom_qq_line()
