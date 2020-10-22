library(codingMatrices)
#library(fractional)

# 6.5
## a.
### manual way
c_0 <- c(1, 0, 0)
c_1 <- c(1, -1, 0)
c_2 <- c(1, 0, -1)
c_3 <- c(0, 1, -1) # redundant
C <- rbind(c_1, c_2)

### automated way
mean_contrasts(contr.treatment(3))

### sample data:
n <- 40
p <- 3
X_bar <- c(46.1, 57.3, 50.4)
S <- rbind(
  c(101.3, 63.0, 71.0),
  c(63.00, 80.2, 55.6),
  c(71.00, 55.6, 97.4)
)

### test statistics
Tsq <- n * t(C %*% X_bar) %*% solve(C %*% S %*% t(C)) %*% (C %*% X_bar)

### critical value
alpha = 0.05
F_star <- qf(1-alpha, p-1, n-p+1) * ((n-1)*(p-1)) / (n-p+1)

## b.
C <- rbind(c_1, c_2, c_3)
X_diff <- C %*% X_bar
MOE <- sqrt(F_star * diag(C %*% S %*% t(C)) / n)
SCI <- cbind(X_diff - MOE, X_diff + MOE)
colnames(SCI) = c("LCL", "UCL")
SCI
C

# 6.8
wilksLambdaTest <- function(p, g, N, B, T, W, alpha = 0.05, largeSample = FALSE){

  ratio <- function(wilksLambda, root = FALSE, ln = FALSE){
    if(root) wilksLambda = sqrt(wilksLambda)
    ratio <- (1-wilksLambda)/wilksLambda
    ifelse(ln, return(log(wilksLambda)), return(ratio))
  }
  
  
  wilksLambda <- det(W)/det(B+W)
  if(largeSample == FALSE){
    if(p == 1 && g >= 2){
      obs <- (N-g)/(g-1) * ratio(wilksLambda)
      ref <- qf(1-alpha, g-1, N-g)
    } else if (p == 2 && g >= 2){
      obs <- (N-g-1)/(g-1) * ratio(wilksLambda, root = TRUE)
      ref <- qf(1-alpha, 2*(g-1), 2*(N-g-1))
    } else if (p >= 2 && g == 2){
      obs <- (N-p-1)/p * ratio(wilksLambda)
      ref <- qf(1-alpha, p, N-p-1)
    } else if (p >= 2 && g == 3){
      obs <- (N-p-2)/p * ratio(wilksLambda, root = TRUE)
      ref <- qf(1-alpha, 2*p, 2*(N-p-2))
    } else{
      obs <- -1*(N-1-(p+g)/2) * ratio(wilksLambda, ln = TRUE)
      ref <-qchisq(1-alpha, p*(g-1))
    }    
  } else{
    obs <- -1*(N-1-(p+g)/2) * ratio(wilksLambda, ln = TRUE)
    ref <-qchisq(1-alpha, p*(g-1))
  }
  return(
    list(
      alpha = alpha,
      wilksLambda = wilksLambda,
      Fobs = obs, 
      Fstar = ref
    )
  )
}


MANOVA <- function(X, ...){
  p <- dim(X[[1]])[2]
  g <- length(X)
  T <- matrix(rep(0, p*p),nrow = p, ncol = p)
  B <- T
  W <- T
  .X <- c()
  for(x in X){.X <- rbind(.X, x)}
  Xbar <- colMeans(.X)
  xbar <- list()
  for(x in X){
    n <- dim(x)[1]
    .xbar <- colMeans(x)
    xbar <- append(xbar, list(.xbar))
    B <- B + n * (.xbar - Xbar) %*% t(.xbar - Xbar)
    for(j in 1:n){
      T <- T + (x[j, ] - Xbar) %*% t(x[j, ] - Xbar)
      W <- W + (x[j, ] - .xbar) %*% t(x[j, ] - .xbar)
    }
  }
  N <- dim(.X)[1]
  df <- list(T = N-1, B = g-1, W = N-g)
  specs <- list(p = p, g = g, N = N)
  results <- wilksLambdaTest(p, g, N, B, T, W, ...)
  wilksLambda = det(W) / (det(B+W))
  return(
    list(
      T = T, 
      B = B, 
      W = W, 
      Xbar = Xbar, 
      xbar = xbar, 
      X = X,
      df = df,
      results = results,
      specs = specs
    )
  )
}


X_1 <- cbind(c(6, 5, 8, 4, 7), c(7, 9, 6, 9, 9))
X_2 <- cbind(c(3, 1, 2), c(3, 6, 3))
X_3 <- cbind(c(2, 5, 3, 2), c(3, 1, 1, 3))
X <- list(X_1, X_2, X_3)
MANOVA(X)

## c.
MANOVA(X, alpha = 0.01)
MANOVA(X, alpha = 0.01, largeSample = TRUE)
