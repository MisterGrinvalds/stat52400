library(ppcc)
library(tidyverse)

# 5.20
## Read Data
file <- file.choose()
df <- read_delim(
  file,
  col_names = FALSE,
  delim = "  "
) %>% 
  mutate(
    X1 = as.integer(trimws(X1)), 
    X2 = as.integer(trimws(X2))
  )

## Ellipse Functions
pathify_vector_coordinates <- function(center, vector_coords){
  path <- center
  for (i in 1:dim(vector_coords)[1]){
    path <- rbind(path, vector_coords[i, ], center)
  }
  return(path)
}


confidence_ellipse <- function(X1, X2, alpha = 0.5, thetas = seq(0, 2*pi, length.out=200), region = FALSE){
  Mu1 <- mean(X1)
  Mu2 <- mean(X2)
  E <- cov(cbind(X1, X2))
  eigen <- eigen(E)
  if (!region){
    threshold <- qchisq(alpha, 2)
  } else{
    n <- length(X1) # dim(X)[1]
    p <- 2          # dim(X)[2]
    threshold <- (p*(n-1))/(n*(n-p))*qf(p = 1-alpha, df1 = p, df2 = n-p)
  }
  ellipse <- sqrt(threshold) * cbind(sqrt(eigen$values[1]) * cos(thetas), sqrt(eigen$values[2]) * sin(thetas)) %*% t(eigen$vectors)
  ellipse <- sweep(ellipse, 2, c(Mu1, Mu2), "+")
  basis <- sqrt(threshold) * eigen$vectors %*% diag(sqrt(eigen$values))
  x_1 <- rbind(Mu1 + basis[1, ], Mu1 - basis[1, ])
  x_2 <- rbind(Mu2 + basis[2, ], Mu2 - basis[2, ])
  vectors <- cbind(c(x_1), c(x_2))
  center <- c(Mu1, Mu2) 
  return(
    list(
      ellipse = ellipse, 
      vector_coords = vectors, 
      center = center, 
      path = pathify_vector_coordinates(center, vectors), 
      threshold = threshold
    )
  )
}


plot_ellipse <- function(plot, ellipse, Mu_size = 3, Mu_shape = 4, Mu_col = "red", eigen_linetype = "dotted", eigen_col = "blue", ...){
  plot <- plot + geom_path(mapping = aes(x = ellipse$ellipse[,1], y = ellipse$ellipse[,2])) + 
    geom_point(mapping = aes(x = ellipse$center[1], y = ellipse$center[2]), size = Mu_size, shape = Mu_shape, col = Mu_col) +
    geom_path(mapping = aes(x = ellipse$path[,1], y = ellipse$path[,2]), linetype = eigen_linetype, col = eigen_col)
  return(plot)
}


## Check if c(190, 275) is in Confidence Region
ellipse <- confidence_ellipse(df$X1, df$X2, alpha = 0.05, region = TRUE)
plot <- qplot(df$X1, df$X2) + xlab("Tail Length") + ylab("Wing Length")
plot_ellipse(plot, ellipse) + geom_point(mapping = aes(x = 190, y = 275), size = 3, shape = 3)
mu <- c(mean(df$X1), mean(df$X2))
E <- cov(cbind(df$X1, df$X2))
X <- c(190, 275)
distance <- t(X - mu) %*% solve(E) %*% c(X - mu) 
distance
ellipse$threshold
distance < ellipse$threshold

## CIs
intervals <- function(X, alpha = 0.05, method = "bonferroni"){
  intervals <- c()
  n <- dim(X)[1]
  p <- dim(X)[2]
  if (method %in% c("bonferroni", "b")){
    threshold <- qt(1-alpha/(2*p), n-1)
  } else if (method %in% c("hotelling", "h")){
    threshold <- sqrt((p*(n-1))/(n-p)*qf(p = 1-alpha, df1 = p, df2 = n-p))
  }
  for (j in 1:p){
    mu <- mean(X[,j])
    MOE <- threshold * sqrt(var(X[,j])/n)
    print(mu)
    print(threshold^2)
    print(var(X[,j])/n)
    print(MOE)
    LCL <- mu - MOE
    UCL <- mu + MOE
    intervals <- rbind(intervals, c(LCL, UCL))
  }
  colnames(intervals) <- c("LCL", "UCL")
  rownames(intervals) <- lapply(1:p, FUN = function(x) paste("X", x, sep = ""))
  return(intervals)
}


intervals(cbind(df$X1, df$X2))
intervals(cbind(df$X1, df$X2), method = "h")

## Check for Normality
for (var in c("X1", "X2")){
  print(ggplot(data = df, mapping = aes_string(sample = var)) + geom_qq() + geom_qq_line())
  print(ggplot(data = df, mapping = aes_string(x = var)) + geom_histogram(bins = 10))
  print(shapiro.test(df[[var]]))
  print(ppccTest(df[[var]]))
  print(boxcox(as.matrix(df[[var]]), p = 1) %>% arrange(desc(likelihoods)) %>% print(n = 1))
}
ellipse <- confidence_ellipse(df$X1, df$X2, alpha = 0.5, region = FALSE)
plot <- qplot(df$X1, df$X2) + xlab("Tail Length") + ylab("Wing Length")
plot_ellipse(plot, ellipse)
distance <- diag(as.matrix(df) %*% solve(cov(df)) %*% t(df))
distance_sort <- sort(distance)
q <- qchisq((1:length(distance_sort) - 0.5)/length(distance_sort), df = 2)
qplot(q, distance_sort) + xlab("theoretical quantiles") + ylab("sample quantiles")


### transform? ###
boxcox <- function(X, p = dim(X)[2], lambdas = rep(list(seq(from = -1, to = 1, by = 0.1)), p)){
  n = dim(X)[1]
  lambdas <- expand.grid(lambdas)
  #if(!is.null(colnames(X))){colnames(lambdas) <- colnames(X)}
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


seq_ <- seq(-2, 4, 0.1)
bxcx <- boxcox(cbind(df$X1, df$X2), lambdas = list(seq_, seq_))
bxcx %>% arrange(desc(likelihoods))
contour(seq_, seq_, matrix(bxcx$likelihoods, nrow = length(seq_)))
X1_new <- (df$X1^0.6-1)/0.6
X2_new <- (df$X2^2.6-1)/2.6
ellipse <- confidence_ellipse(X1_new, X2_new, alpha = 0.5, region = FALSE)
plot <- qplot(X1_new, X2_new) + xlab("Tail Length") + ylab("Wing Length")
plot_ellipse(plot, ellipse)
X <- cbind(X1_new, X2_new)
distance <- diag(X) %*% solve(cov(X)) %*% t(X)
distance_sort <- sort(distance)
q <- qchisq((1:length(distance_sort) - 0.5)/length(distance_sort), df = 2)
qplot(q, distance_sort) + xlab("theoretical quantile") + ylab("sample_quantile")
