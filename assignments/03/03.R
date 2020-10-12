library(plotly)
library(ppcc)
library(tidyverse)

# 4.28
file <- file.choose()
df <- read_delim(
  delim = "  ",
  file,
  col_names = c(
    "wind",
    "solar_radiation",
    "co",
    "no",
    "no2",
    "o3",
    "hc"
  )
)
df <- as.tibble(apply(df, 2, function(x) as.integer(str_trim(x, side = c("both")))))
ggplot(data = df, mapping = aes(sample = solar_radiation)) + geom_qq() + geom_qq_line()
x_sorted <- sort(df$solar_radiation)
q <- qnorm((1:length(x_sorted) - 0.5) / length(x_sorted))
SSx <- sum((x_sorted - mean(x_sorted))^2)
SSq <- sum((q - mean(q))^2)
SSxq <- sum((x_sorted - mean(x_sorted)) * (q - mean(q)))
r_statistic <- SSxq / sqrt(SSx * SSq)

# 4.29
## univariate normality and transformations
boxcox <- function(x, lambdas = seq(from = -1, to = 1, by = 0.05), variable = "x"){
  likelihoods <- c()
  n = length(x)
  for (lambda in lambdas){
    if (lambda != 0){x_lambda = (x^lambda - 1)/lambda} else{x_lambda <- log(x)}
    xbar = mean(x_lambda)
    likelihood <- (-n/2) * log((1/n)*sum((x_lambda-xbar)^2)) + (lambda-1)*sum(log(x))
    likelihoods <- c(likelihoods, likelihood)
  }
  plot <- ggplot(data = NULL, mapping = aes(x = lambdas, y = likelihoods)) + geom_line()
  print(plot)
  return(tibble(lambdas, likelihoods))
}


run_diagnostics <- function(df, vars){
  for (var in vars){
    print(ggplot(data = df, mapping = aes_string(x = var)) + geom_histogram(bins = 10))
    print(ggplot(data = df, mapping = aes_string(sample = var)) + geom_qq() + geom_qq_line())
    print(shapiro.test(df[[var]]))
    print(ppccTest(df[[var]]))
    print(boxcox(df[[var]]) %>% arrange(desc(likelihoods)) %>% print(n = 1))
  }
}


### test for univariate normality
vars = list("no2", "o3")
run_diagnostics(df, vars)
plot <- qplot(df$no2, df$o3) + xlab("NO_2") + ylab("O_3")
plot

### transform the data
df_new <- tibble(no2_new = (df$no2^-0.1 - 1)/-0.1, o3_new = (df$o3^0.25 - 1)/0.25)
vars = list("no2_new", "o3_new")
run_diagnostics(df_new, vars)
plot <- qplot(df_new$no2_new, df_new$o3_new) + xlab("(NO_2^-0.01 - 1)/-0.01") + ylab("(O_3^0.025 - 1)/0.025")
plot

## a. statistical distances
### untransformed
x_c <- cbind(df$no2 - mean(df$no2), df$o3 - mean(df$o3))
E <- cov(x_c)
distance <- diag(x_c %*% solve(E) %*% t(x_c))
sort(distance)

### transformed
x_c_new <- cbind(df_new$no2_new - mean(df_new$no2_new), df_new$o3_new - mean(df_new$o3_new))
E <- cov(x_c_new)
distance_new <- diag(x_c_new %*% solve(E) %*% t(x_c_new))
sort(distance_new)

## b. confidence region
pathify_vector_coordinates <- function(center, vector_coords){
  path <- center
  for (i in 1:dim(vector_coords)[1]){
    path <- rbind(path, vector_coords[i, ], center)
  }
  return(path)
}


confidence_ellipse <- function(X1, X2, alpha = 0.5, thetas = seq(0, 2*pi, length.out=200)){
  Mu1 <- mean(X1)
  Mu2 <- mean(X2)
  E <- cov(cbind(X1, X2))
  eigen <- eigen(E)
  ellipse <- sqrt(qchisq(alpha, 2)) * cbind(sqrt(eigen$values[1]) * cos(thetas), sqrt(eigen$values[2]) * sin(thetas)) %*% t(eigen$vectors)
  ellipse <- sweep(ellipse, 2, c(Mu1, Mu2), "+")
  basis <- sqrt(qchisq(alpha, 2)) * eigen$vectors %*% diag(sqrt(eigen$values))
  x_1 <- rbind(Mu1 + basis[1, ], Mu1 - basis[1, ])
  x_2 <- rbind(Mu2 + basis[2, ], Mu2 - basis[2, ])
  vectors <- cbind(c(x_1), c(x_2))
  center <- c(Mu1, Mu2) 
  return(list(ellipse = ellipse, vector_coords = vectors, center = center, path = pathify_vector_coordinates(center, vectors)))
}


plot_ellipse <- function(plot, ellipse, Mu_size = 3, Mu_shape = 4, Mu_col = "red", eigen_linetype = "dotted", eigen_col = "blue", ...){
  plot <- plot + geom_path(mapping = aes(x = ellipse$ellipse[,1], y = ellipse$ellipse[,2])) + 
    geom_point(mapping = aes(x = ellipse$center[1], y = ellipse$center[2]), size = Mu_size, shape = Mu_shape, col = Mu_col) +
    geom_path(mapping = aes(x = ellipse$path[,1], y = ellipse$path[,2]), linetype = eigen_linetype, col = eigen_col)
  return(plot)
}


### untransformed
plot <- qplot(df$no2, df$o3) + xlab("NO_2") + ylab("O_3")
ellipse <- confidence_ellipse(df$no2, df$o3)
plot_ellipse(plot, ellipse)
sum(ifelse(distance > qchisq(0.5, 2), 0, 1))
### transformed
plot <- qplot(df_new$no2_new, df_new$o3_new) + xlab("(NO_2^-0.01 - 1)/-0.01") + ylab("(O_3^0.025 - 1)/0.025")
ellipse <- confidence_ellipse(df_new$no2_new, df_new$o3_new)
plot_ellipse(plot, ellipse)
sum(ifelse(distance_new > qchisq(0.5, 2), 0, 1))

## c. chi-square plot
### untransformed
distance_sort <- sort(distance)
q <- qchisq((1:length(distance_sort) - 0.5)/length(distance_sort), df = 2)
qplot(q, distance_sort) + xlab("theoretical quantile") + ylab("sample_quantile")

### transformed
distance_sort <- sort(distance_new)
q <- qchisq((1:length(distance_sort) - 0.5)/length(distance_sort), df = 2)
qplot(q, distance_sort) + xlab("theoretical quantile") + ylab("sample_quantile")

# 4.30
X1 <- c(1, 2, 3, 3, 4, 5, 6, 8, 9, 11)
X2 <- c(18.95, 19.00, 17.95, 15.54, 14.00, 12.95, 8.94, 7.49, 6.00, 3.99)
df <- tibble(X1 = X1, X2 = X2)

## a.
ggplot(data = df, mapping = aes(sample = X1)) + geom_qq() + geom_qq_line()
boxcox(X1, lambdas = seq(-1, 1, 0.01)) %>% arrange(desc(likelihoods))
X1_new <- (X1^0.37 - 1)/0.37
ggplot(data = df, mapping = aes(sample = X1_new)) + geom_qq() + geom_qq_line()

## b.
ggplot(data = df, mapping = aes(sample = X2)) + geom_qq() + geom_qq_line()
boxcox(X2, lambdas = seq(-1, 2, 0.01)) %>% arrange(desc(likelihoods))
X2_new <- (X1^0.94 - 1)/0.94
ggplot(data = df, mapping = aes(sample = X2_new)) + geom_qq() + geom_qq_line()

## c.
boxcox <- function(X, p = dim(X)[2], lambdas = rep(list(seq(from = -1, to = 1, by = 0.1)), p)){
  n = dim(X)[1]
  lambdas <- expand.grid(lambdas)
  colnames(lambdas) <- colnames(X)
  likelihoods <- c()
  for (i in 1:dim(lambdas)[1]){
    X_temp <- c()
    log_sum_total <- c()
    for (j in 1:p){
      if(lambdas[i,j] != 0){
        X_j_new <- (X[,j] ^ lambdas[i,j] - 1)/lambdas[i,j]
      } else{
        X_j_new <- log(X[,j])
      }
      X_temp <- cbind(X_temp, X_j_new)
      log_sum <- (lambdas[i,j] - 1) * sum(log(X[,j]))
      log_sum_total <- c(log_sum_total, log_sum)
    }
    likelihood <- (-n/2) * log(det(cov(X_temp))) + sum(log_sum_total)
    likelihoods <- c(likelihoods, likelihood)
  }
  return(as.tibble(cbind(lambdas, likelihoods)))
}

seq_ <- seq(-2, 2, 0.01)
boxcox(cbind(X1), lambdas = seq_) %>% arrange(desc(likelihoods))
boxcox(cbind(X2), lambdas = seq_) %>% arrange(desc(likelihoods))
seq_ <- seq(-3, 3, 0.05)
bxcx <- boxcox(cbind(X1, X2), lambdas = list(seq_, seq_))
bxcx %>% arrange(desc(likelihoods))
contour(seq_, seq_, matrix(bxcx$likelihoods, nrow = length(seq_)))
bxcx_3D <- plot_ly(z=matrix(bxcx$likelihoods, ncol = length(seq_)))
bxcx_3D %>% add_surface()        
        
