library(GGally)
library(mvtnorm)
library(plotly)
library(ppcc)
library(tidyverse)

# Part I
file <- file.choose()
df <- read_tsv(file)
df$species <- as.factor(df$species)

## Clean Data
diagnosticsBoxCox <- function(
  df,
  factors,
  lambdas = list(seq(from = -3, to = 3, by = 0.1)), 
  title = NULL
){
  
  boxCox <- function(df, lambdas, title){
    X <- as.matrix(df)
    n <- dim(X)[1]
    p <- dim(X)[2]
    lambdas <- rep(lambdas, p) 
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
  
  results <- list()
  results[["seq"]] <- unlist(lambdas)
  results[["title"]] <- title
  if (factors %in% colnames(df)){
    results[["all"]] <- boxCox(df %>% select(-factors), lambdas, title)
    for (level in unique(df[[factors]])){
      results[[paste(as.character(level))]] <- boxCox(filterCategory(df, factors, level), lambdas, title)
    }  
  } else{
    results[["all"]] <- boxCox(df, lambdas, title)
  }
  return(results)
}


diagnosticsQQChiSq <- function(
  df, 
  factors, 
  alpha = 0.01, 
  opacity = 0.5
){
  df.tmp <- df %>% select(-factors)
  X.tmp <- apply(df.tmp, 2, function(x) x - mean(x))
  n <- dim(X.tmp)[1]
  p <- dim(X.tmp)[2]
  Z <-apply(X.tmp, 2, function(x) x / sd(x))
  df.tmp <- df.tmp %>% mutate(
    distance = diag(X.tmp %*% solve(cov(df.tmp)) %*% t(X.tmp)),
    greq_chisq = ifelse(distance > qchisq(1-alpha, p), distance, as.numeric(NA))
  )
  df.tmp$greq_chisq[c(which(!is.na(df.tmp$greq_chisq)))] = c(which(!is.na(df.tmp$greq_chisq)))
  df.tmp[[factors]] <- df[[factors]]
  df.tmp <- df.tmp %>% arrange(distance) %>% mutate(quantile = qchisq((1:n - 0.5)/n, p)) 
  plot <- ggplot(data = df.tmp, aes_string(colour = factors, alpha = opacity)) + 
    geom_point(aes(x = quantile,  y = distance)) +
    xlab("chi-square quantile") + 
    ylab("generalized squared distance") +
    geom_hline(yintercept = qchisq(1 - alpha, p), linetype = "dashed") +
    geom_text(aes(label = greq_chisq, x = quantile,  y = distance), hjust = -0.3, na.rm = TRUE, show.legend = FALSE)
  return(list(plot = plot, df = df, Z = Z))
}


diagnosticsUnivariateNormal <- function(
  df,
  factors
){
  
  diagnostics <- function(df, vars){
    diagnostics <- list()
    for (var in vars){
      diagnostics[[var]][["histogram"]] <- ggplot(data = df, mapping = aes_string(x = var)) + geom_histogram(bins = 10)
      diagnostics[[var]][["QQ"]] <- ggplot(data = df, mapping = aes_string(sample = var)) + geom_qq() + geom_qq_line()
      diagnostics[[var]][["shapiroWilk"]] <- shapiro.test(df[[var]])
      diagnostics[[var]][["correlation"]] <- ppccTest(df[[var]])
    }
    return(diagnostics)
  }
  
  vars <- colnames(df)[which(!(colnames(df)) %in% factors)] 
  results <- list()
  if (factors %in% colnames(df)){
    results[["all"]] <- diagnostics(df %>% select(-factors), vars)
    for (level in unique(df[[factors]])){
      results[[paste(as.character(level))]] <- diagnostics(filterCategory(df, factors, level), vars)
    }  
  } else{
    results[["all"]] <- diagnostics(df, vars)
  }
  return(results)
}


diagnosticsZ <- function(
  df, 
  factors, 
  threshold = 2.5
){
  df.tmp <- df %>% select(-factors)
  Z <- as.tibble(apply(df.tmp, 2, function(x) (x - mean(x)) / sd(x))) %>% 
    filter_all(any_vars(abs(.) > threshold))
  Z <- list("all" = Z)
  for (level in unique(df[[factors]])){
    df.tmp <- filterCategory(df, factors = factors, level = level)
    Z[[level]] <- as.tibble(apply(df.tmp, 2, function(x) (x - mean(x)) / sd(x))) %>% 
      filter_all(any_vars(abs(.) > threshold))
  }
  return(Z)
}


filterCategory <- function(
  df, 
  factors, 
  level
){
  expression <- paste(factors, '==', level)
  df.tmp <- df %>% 
    filter(rlang::eval_tidy(rlang::parse_expr(expression))) %>% 
    select(-factors)
  return(df.tmp)
}


markOutlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}


plotBoxCox <- function(
  results,
  xlab = NULL,
  ylab = NULL
){
  p <- dim(results$all)[2] - 1
  plots <- list()
  groups <- names(results)[which(!(names(results) %in% "seq"))]
  if (p > 2){
    print("Cannot plot, too many dimensions.")
  } else if (p == 1) {
    for (group in groups){
      plots[[as.character(group)]] <- ggplot(data = NULL, mapping = aes(x = results[[group]]$Var1, y = results[[group]]$likelihoods)) + 
        geom_line() + 
        ggtitle(paste("Box-Cox Results for Group: ", as.character(group), sep =""))
    }
  } else if (p == 2) {
    for (group in groups){
      plot <- ggplot(data = results[[group]], mapping = aes(x = Var1, y = Var2, z = likelihoods)) + 
        geom_contour_filled() + 
        ggtitle(paste("Box-Cox Results for Group: ", as.character(group), sep =""))
      if(!is.null(xlab)) {plot <- plot + xlab(xlab)}
      if(!is.null(ylab)) {plot <- plot + ylab(ylab)}
      plots[[as.character(group)]] <- plot
    }
  }
  return(plots)
} 

df.red <- df %>% select(-X3, -X4)
diagnosticsUnivariateNormal(df.red, "species")
diagnosticsZ(df.red, "species", threshold = 3)
ggpairs(df.red, aes(colour = species, alpha = 0.5))
box_cox <-diagnosticsBoxCox(df.red, "species")
plotBoxCox(box_cox)
diagnosticsQQChiSq(df.red, "species", alpha = 0.005)

## 1. MANOVA
testBoxM <- function(
  S,
  g, 
  n, 
  p, 
  alpha = 0.05
){
  n <- unlist(n)
  M <- (sum(n) - g) * log(det(S[["pooled"]])) - sum(n * log(unlist(lapply(S[1:g], function(x) det(x)))))
  u <- (sum(n^-1) - 1/(sum(n)-g))*(2*p^2+3*p-1)/(6*(p+1)*(g-1))
  C <- M * u
  v <- 0.5*p*(p+1)*(g-1)
  chisq_star <- qchisq(1 - alpha, v)
  return(
    list(M = M, u = u, C = C, v = v, chisq_star = chisq_star)
  )
}


testPillai <- function(
  B, 
  W, 
  df, 
  p, 
  alpha = 0.05
){
  pillaiTrace <- sum(diag(B * solve(B + W)))
  s <- min(df$B, p)
  n <- (df$W - p - 1)/2
  m <-  (abs(df$B - p) - 1)/2
  Fobs <- (2*n+s+1)/(2*m+s+1) * pillaiTrace/(p-pillaiTrace)
  Fstar <- qf(1 - alpha, s*(2*m+s+1), s*(2*n+s+1))
  return(
    list(pillaiTrace = pillaiTrace, Fobs = Fobs, Fstar = Fstar)
  )
}


testWilksLambda <- function(
  B, 
  T, 
  W, 
  g, 
  N, 
  p, 
  alpha = 0.05, 
  largeSample = FALSE
){
  
  ratio <- function(wilksLambda, root = FALSE, ln = FALSE){
    if(root) wilksLambda = sqrt(wilksLambda)
    ratio <- (1-wilksLambda)/wilksLambda
    ifelse(ln, return(log(wilksLambda)), return(ratio))
  }
  
  
  wilksLambda <- det(W)/det(B+W)
  if(largeSample == FALSE){
    if(p == 1 && g >= 2){
      coef <- (N-g)/(g-1)
      obs <-  ratio(wilksLambda)
      ref <- qf(1-alpha, g-1, N-g)
    } else if (p == 2 && g >= 2){
      coef <- (N-g-1)/(g-1)
      obs <- ratio(wilksLambda, root = TRUE)
      ref <- qf(1-alpha, 2*(g-1), 2*(N-g-1))
    } else if (p >= 2 && g == 2){
      coef <- (N-p-1)/p
      obs <- ratio(wilksLambda)
      ref <- qf(1-alpha, p, N-p-1)
    } else if (p >= 2 && g == 3){
      coef <- (N-p-2)/p
      obs <- ratio(wilksLambda, root = TRUE)
      ref <- qf(1-alpha, 2*p, 2*(N-p-2))
    } else{
      coef <- -1*(N-1-(p+g)/2)
      obs <- ratio(wilksLambda, ln = TRUE)
      ref <- qchisq(1-alpha, p*(g-1))
    }    
  } else{
    coef <- -1*(N-1-(p+g)/2)
    obs <- ratio(wilksLambda, ln = TRUE)
    ref <-qchisq(1-alpha, p*(g-1))
  }
  return(
    list(
      alpha = alpha,
      wilksLambda = wilksLambda,
      ratio = obs,
      coef = coef,
      Fobs = obs * coef,
      Fstar = ref
    )
  )
}


MANOVA <- function(
  df, 
  factors, 
  alpha = 0.05, 
  ...
){
  n <- list()
  N <- dim(df)[1]
  S <- list()
  .X <- df %>% select(-factors)
  X <- list()
  Xbar <- colMeans(.X)
  xbar <- list()
  for (level in unique(df[[factors]])){
    X[[level]] <- as.matrix(filterCategory(df, factors, level))
    n[[level]] <- dim(X[[level]])[1]
    S[[level]] <- cov(X[[level]])
    xbar[[level]] <- colMeans(filterCategory(df, factors, level))
  }
  p <- dim(X[[1]])[2]
  g <- length(X)
  S[["S"]]  <- cov(df %>% select(-factors))
  S[["pooled"]] <- matrix(data = 0, ncol = p, nrow = p) 
  for (group in 1:g){S[["pooled"]] <- S[["pooled"]] + (dim(X[[group]])[1] - 1) * S[[group]]}
  S[["pooled"]] <- S[["pooled"]] /(N-g)
  
  T <- matrix(rep(0, p*p), nrow = p, ncol = p)
  B <- matrix(rep(0, p*p), nrow = p, ncol = p)
  W <- matrix(rep(0, p*p), nrow = p, ncol = p)
  for(x in X){
    .n <- dim(x)[1]
    .xbar <- colMeans(x)
    B <- B + .n * (.xbar - Xbar) %*% t(.xbar - Xbar)
    for(j in 1:.n){
      T <- T + (x[j, ] - Xbar) %*% t(x[j, ] - Xbar)
      W <- W + (x[j, ] - .xbar) %*% t(x[j, ] - .xbar)
    }
  }
  
  df <- list(
    T = N-1, 
    B = g-1, 
    W = N-g
  )
  
  specs <- list(p = p, g = g, N = N, n = n)
  pillai <- testPillai(B, W, df, p, alpha)
  wilks <- testWilksLambda(B, T, W, g, N, p, alpha)
  boxM <- testBoxM(S, g, n, p, alpha)
  return(
    list(
      T = T, 
      B = B, 
      W = W,
      Xbar = Xbar, 
      xbar = xbar, 
      X = X,
      S = S,
      df = df,
      boxM = boxM,
      pillai = pillai,
      wilks = wilks,
      specs = specs
    )
  )
}
  
df.red <- df %>% select(-X3, -X4)
manova <- MANOVA(df.red, "species")
manova

## 2. describe the test and assumptions

## 3. H_0: Mu_1 == Mu_2
compareMeanVectors <- function(
  manova, 
  i, 
  j, 
  alpha = 0.05,
  delta_0 = c(0, 0),
  equal_S = TRUE,
  confint_method = "h"
){
  g <- manova$specs$g
  n <- manova$specs$n
  p <- manova$specs$p
  S <- manova$S
  delta <- t(manova$xbar[[i]] - manova$xbar[[j]]) - delta_0
  
  if (equal_S){
    .S <- (1/n[[i]] + 1/n[[j]]) * S$pooled
    obs <- delta %*% solve(.S) %*% t(delta)
    star <- qf(1 - alpha, p, n[[i]] + n[[j]] - p - 1)
    coef <- (n[[i]] + n[[j]] - 2) * p / (n[[i]] + n[[j]] - p - 1)
    c_2 <- coef * star
    test <- list(
        Tsq_obs = obs,
        F_star = star,
        coef = coef,
        critical_value = c_2
    )
  } else{
    .S <- 1/n[[i]]*S[[i]] + 1/n[[j]]*S[[j]]
    obs <- delta %*% solve(.S) %*% t(delta)
    star <- qchisq(1-alpha, p)
    c_2 <- star
    test <- list(
      Tsq_obs = obs,
      chisq_star = star
    )
  }
  
  intervals = list()
  if ("h" %in% confint_method){
    MOE <- sqrt(c_2) * sqrt(diag(.S))
    intervals[["h"]] <- list(
      center = delta,
      c_2 = c_2,
      MOE = MOE,
      intervals = cbind(LCL = c(delta - MOE), UCL = c(delta + MOE))
    )
  } 
  if ("b" %in% confint_method) {
    .c_2 <- qt(1-alpha/(2*p), sum(n[[i]] + n[[j]]) - g)
    MOE <- sqrt(.c_2) * sqrt(diag(.S))
    intervals[["b"]] <- list(
      center = delta,
      c_2 = .c_2,
      MOE = MOE,
      intervals = cbind(LCL = c(delta - MOE), UCL = c(delta + MOE))
    )
  }
  
  region <- list(
    center = delta,
    c_2 = c_2,
    .S = .S,
    MOE = c_2 * .S
  )
  
  return(
    list(
      alpha = alpha,
      intervals = intervals,
      manova = manova,
      methods = confint_method,
      region = region,
      test = test
    )
  )
}

df.red <- df %>% filter(species != 3) %>% select(-X3, -X4)
manova <- MANOVA(df.red, "species")
manova
bartlett.test(X1 + X2 ~ species, data = df.red)
comparison <- compareMeanVectors(manova, 1, 2, confint_method = c("b", "h"))
comparison

## 4.

makeConfidenceEllipse <- function(
  mean_comparison,
  thetas = seq(0, 2*pi, length.out=200)
){
  g <- mean_comparison$manova$specs$g
  n <- mean_comparison$manova$specs$n
  p <- mean_comparison$manova$specs$p
  .S <- mean_comparison$region$.S
  center <- mean_comparison$region$center
  c_2 <- mean_comparison$region$c_2
  eigen <- eigen(.S)
  
  ellipse <- sqrt(c_2) * cbind(sqrt(eigen$values[1]) * cos(thetas), sqrt(eigen$values[2]) * sin(thetas)) %*% t(eigen$vectors)
  ellipse <- sweep(ellipse, 2, c(center[[1]], center[[2]]), "+")
  basis <- sqrt(c_2) * eigen$vectors %*% diag(sqrt(eigen$values))
  x_1 <- rbind(center[[1]] + basis[1, ], center[[1]] - basis[1, ])
  x_2 <- rbind(center[[2]] + basis[2, ], center[[2]] - basis[2, ])
  vectors <- cbind(c(x_1), c(x_2))
  center <- c(center[[1]], center[[2]]) 
  return(
    list(
      ellipse = ellipse, 
      vector_coords = vectors, 
      center = center, 
      path = pathifyVectorCoordinates(center, vectors)
    )
  )
}


pathifyVectorCoordinates <- function(
  center, 
  vector_coords
){
  path <- center
  for (i in 1:dim(vector_coords)[1]){
    path <- rbind(path, vector_coords[i, ], center)
  }
  return(path)
}


plotConfidenceEllipse <- function(
  ellipse, 
  plot, 
  Mu_size = 3, 
  Mu_shape = 4, 
  Mu_col = "red", 
  eigen_linetype = "dotted", 
  eigen_col = "blue", ...
){
  plot <- plot + 
    geom_path(mapping = aes(x = ellipse$ellipse[,1], y = ellipse$ellipse[,2])) + 
    geom_point(mapping = aes(x = ellipse$center[1], y = ellipse$center[2]), size = Mu_size, shape = Mu_shape, col = Mu_col) +
    geom_path(mapping = aes(x = ellipse$path[,1], y = ellipse$path[,2]), linetype = eigen_linetype, col = eigen_col)
  return(plot)
}


plotConfidenceIntervals <- function(
  mean_comparison,
  plot,
  colours = c("green", "darkgreen")
){
  colour_i <- 1
  for (i in mean_comparison$methods){
    intervals <- mean_comparison$intervals[[i]]$intervals
    plot <- plot + 
      geom_vline(xintercept = c(intervals[1,1], intervals[1,2]), linetype = 2, color = colours[colour_i]) +
      geom_hline(yintercept = c(intervals[2,1], intervals[2,2]), linetype = 2, color = colours[colour_i])
    colour_i <- colour_i + 1
  }
  return(plot)
}


ellipse <- makeConfidenceEllipse(comparison)
plot <- ggplot(data = NULL) + xlab("X1") + ylab("X2")
plot <- plotConfidenceEllipse(ellipse, plot)
plot <- plotConfidenceIntervals(comparison, plot)
plot


# Part II.
ggpairs(df)

## 5.
pcaReport <- function(eigen){
  print(
    cbind(
      eigenvalue = eigen$values, 
      total_proportion = cumsum(eigen$values) / sum(eigen$values)
    )
  )
}


screePlot <- function(
  eigen,
  icon_color = "red",
  icon_shape = "triangle",
  icon_size = 3,
  title = NULL
){
  plot <- ggplot(data = NULL, mapping = aes(x = 1:4)) + 
    geom_line(aes(y = eigen$values)) +
    geom_point(aes(y = eigen$values), shape = icon_shape, size = icon_size, colour = icon_color) +
    xlab("index") + 
    ylab("eigenvalue") +
    ggtitle(title)
  return(plot)
}

df.red <- df %>% select(-X1, -X2)
diagnosticsUnivariateNormal(df.red, "species")
box_cox <- diagnosticsBoxCox(df.red, "species")
plotBoxCox(box_cox, xlab = "X3", ylab = "X4")
ggpairs(df, aes(colour = species, alpha = 0.5))
diagnosticsQQChiSq(df, "species", alpha = 0.005)

X <- as.matrix(df[, 1:4])
S <- cov(X)
eigen.S <- eigen(S)
pcaReport(eigen.S)
screePlot(eigen.S)
U <- eigen.S$vectors[, 1:2]
X.pca <- as.tibble(cbind(X %*% U))
X.pca$species <- df$species
ggplot(X.pca, aes(V1, V2, colour = species, alpha = 0.5)) +
  geom_point()
diagnosticsQQChiSq(X.pca, "species", alpha = 0.005)
diagnosticsUnivariateNormal(X.pca, "species")
U.prime <- eigen.S$vectors[, 3:4]
X.pca.prime <- as.tibble(cbind(X %*% U.prime))
colnames(X.pca.prime) <- c("V3", "V4")
X.pca.prime %>% filter_all(any_vars(abs(.) > 1))
X.pca.prime$species <- df$species
ggplot(X.pca.prime, aes(V3, V4, colour = species, alpha = 0.5)) +
  geom_point()

#3d?
U <- eigen.S$vectors[,1:3]
X.pca <- as.tibble(cbind(X %*% U))
X.pca$species <- df$species
fig <- plot_ly(X.pca, x = ~V1, y = ~V2, z = ~V3, color = ~species, colors = c('red', 'green', 'blue'))
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'V1'),
                                   yaxis = list(title = 'V2'),
                                   zaxis = list(title = 'V3')))
fig

# Correlation version, not needed
Z <- t(t(apply(X, 2, function(x) (x - mean(x))/sd(x))))
R <- cov(Z)
eigen.R <- eigen(R)
pcaReport(eigen.R)
screePlot(eigen.R)
U <- eigen.R$vectors[,1:2]
Z.pca <- as.tibble(Z %*% U)
Z.pca$species <- df$species
ggplot(Z.pca, aes(V1, V2, colour = species, alpha = 0.5)) +
  geom_point()

# 6.
manova <- MANOVA(X.pca, "species")
bartlett.test(V1 + V2 ~ species, X.pca)

## a.
priors <- unlist(manova$specs$n) / manova$specs$N
xbars <- manova$xbar
Spool <- manova$S$pooled
for (i in 1:3){
  X.pca[[paste("d_", as.character(i), sep="")]] <- c(t(xbars[[i]]) %*% solve(Spool) %*% t(X.pca[,1:2]) - rep(0.5 * t(xbars[[i]]) %*% solve(Spool) %*% xbars[[i]] + log(priors[i]), 150))
}
class <- c()
for (row in 1:150){
  x <- unlist(X.pca[row, 4:6])
  class <- c(class, order(x, decreasing = TRUE)[1])
}
X.pca <- X.pca %>% mutate(class = class, correct = species == class)
X.pca %>% print(n=150)
sum(!X.pca$correct)

## b.
table <- matrix(0, nrow = 3, ncol = 3)
df.tmp <- X.pca %>% select(species, class)
for (i in 1:150){
  cell <- c(unlist(df.tmp[i,]))
  table[cell[1], cell[2]] <- table[cell[1], cell[2]] + 1
}
colnames(table) <- c("setosa", "versicolor", "virginica")
rownames(table) <- c("setosa", "versicolor", "virginica")
table
6/150

## c.
new_obs <- c(6.7, 2.5, 5.8, 1.8)  %*% U
d <- c()
for (i in 1:3){
  d <- c(d, t(xbars[[i]]) %*% solve(Spool) %*% t(new_obs) - 0.5 * t(xbars[[i]]) %*% solve(Spool) %*% xbars[[i]] + log(priors[i]))
}
d
ggplot(X.pca, aes(V1, V2, colour = species, alpha = 0.5)) +
  geom_point() +
  geom_point(x = new_obs[1], y = new_obs[2], shape = "triangle", size = 4, show.legend = FALSE)
