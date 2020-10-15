library(ppcc)
library(tidyverse)

# 1.
E <- rbind(c(3, -3, 0), c(-3, 12, 0), c(0, 0, 6))

# 2.
## a.
MUx <- c(0, -1, 2)
Ex <- rbind(c(2, -2, 1), c(-2, 5, -2), c(1, -2, 1))
A <- rbind(c(1, 0, -1), c(0, 1, -1))
A %*% MUx
A %*% Ex %*% t(A)

## b.
A <- rbind(c(1, 0.5, 0), c(0, 0, 1))
A %*% MUx
A %*% Ex %*% t(A)

## c.
A <- rbind(c(1, 0, 0), c(0, 1, 2))
A %*% MUx
A %*% Ex %*% t(A)

# 4.
file <- file.choose()
df <- read_tsv(file)

## normality and outliers
boxcox <- function(X, p = dim(X)[2], lambdas = rep(list(seq(from = -3, to = 3, by = 0.1)), p), title = NULL){
  n = dim(X)[1]
  lambdas <- expand.grid(lambdas)
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
  if (p == 1){print(ggplot(data = NULL, mapping = aes(x = lambdas$Var1, y = likelihoods)) + geom_line() + xlab(title))}
  return(as.tibble(cbind(lambdas, likelihoods)))
}


normal_diagnostics <- function(df, vars){
  for (var in vars){
    print(ggplot(data = df, mapping = aes_string(x = var)) + geom_histogram(bins = 10))
    print(ggplot(data = df, mapping = aes_string(sample = var)) + geom_qq() + geom_qq_line())
    print(shapiro.test(df[[var]]))
    print(ppccTest(df[[var]]))
    print(boxcox(as.matrix(df[[var]]), title = var) %>% arrange(desc(likelihoods)) %>% print(n = 1))
  }
}


is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}


vars <- c("Length", "Width")
normal_diagnostics(df.new, vars)
for (var in vars){
  new_var = paste("outlier_", var, sep = "")
  df <- df %>% mutate(!!new_var := ifelse(is_outlier(df[[var]]), df[[var]], as.numeric(NA)))
  df[[new_var]][c(which(!is.na(df[[new_var]])))] = c(which(!is.na(df[[new_var]])))
  print(ggplot(data = df, mapping = aes_string(y = var)) + 
          geom_boxplot(width = 0.2) + 
          theme(aspect.ratio = 5/2) + 
          geom_text(aes_string(x = 0, label = new_var), na.rm = TRUE, hjust = -0.3)
  )
  new_var = paste("std_", var, sep = "")
  df <- df %>% mutate(!!new_var := (df[[var]] - mean(df[[var]]))/sd(df[[var]]))
}
x_center <- t(apply(cbind(df$Length, df$Width), 1, function(x) x - colMeans(cbind(df$Length, df$Width))))
df <- df %>% mutate(
  distance = diag(x_center %*% solve(cov(cbind(df$Length, df$Width))) %*% t(x_center)),
  greq_chisq = ifelse(distance > qchisq(0.99, 2), distance, as.numeric(NA))
)
df %>% select(Length, std_Length, Width, std_Width)%>% filter(abs(std_Width) > 2 | abs(std_Length) > 2)
df$greq_chisq[c(which(!is.na(df$greq_chisq)))] = c(which(!is.na(df$greq_chisq)))
qplot(Length, Width, data = df) + 
  geom_text(aes(x = Length, label = greq_chisq), na.rm = TRUE, hjust = -0.3)
df.new <- df[-42,]

### transformation
seq_ <- seq(-5, 5, 0.125)
bxcx <- boxcox(X = cbind(df.new$Length, df.new$Width), lambdas = rep(list(seq_), 2))
bxcx %>% arrange(desc(likelihoods))
contour(seq_, seq_, matrix(bxcx$likelihoods, nrow = length(seq_)))
df.red$Width <- (df.new$Width^-1.25 -1) / -1.25
plot <- qplot(Length, Width, data = df.red)
plot

## Confidence Intervals
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
    LCL <- mu - MOE
    UCL <- mu + MOE
    intervals <- rbind(intervals, c(LCL, UCL))
  }
  colnames(intervals) <- c("LCL", "UCL")
  rownames(intervals) <- lapply(1:p, FUN = function(x) paste("X", x, sep = ""))
  return(intervals)
}


CIT2 <- intervals(cbind(df.red$Length, df.red$Width), method = "h")
CIBF <- intervals(cbind(df.red$Length, df.red$Width))
CIT2 # T^2 Confidence Intervals
CIBF # Bonferroni Confidence Intervals

## Confidence Ellipse
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
  plot <- plot + geom_path(data = NULL, mapping = aes(x = ellipse$ellipse[,1], y = ellipse$ellipse[,2])) + 
    geom_point(data = NULL, mapping = aes(x = ellipse$center[1], y = ellipse$center[2]), size = Mu_size, shape = Mu_shape, col = Mu_col) +
    geom_path(data = NULL, mapping = aes(x = ellipse$path[,1], y = ellipse$path[,2]), linetype = eigen_linetype, col = eigen_col)
  return(plot)
}


plot <- qplot(df.red$Length, df.red$Width) + xlab("Length") + ylab("Width")
ellipse <- confidence_ellipse(df.red$Length, df.red$Width, alpha = 0.05, region = TRUE)
plot_ellipse(plot, ellipse) + 
  geom_vline(aes(xintercept = c(CIT2[1,1], CIT2[1,2])), linetype = 2, color = "darkgreen") +
  geom_hline(aes(yintercept = c(CIT2[2,1], CIT2[2,2])), linetype = 2, color = "darkgreen") +
  geom_vline(aes(xintercept = c(CIBF[1,1], CIBF[1,2])), linetype = 2, color = "green") +
  geom_hline(aes(yintercept = c(CIBF[2,1], CIBF[2,2])), linetype = 2, color = "green")
  

# 5.
eigen(rbind(c(1, -1), c(-1, 1)))