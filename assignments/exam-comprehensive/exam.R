library(tidyverse)
library(car)
library(caret)
library(GGally)
library(heplots)
library(MASS)
#library(mvtnorm)
library(ppcc)
library(rstatix)

# read data
file <- file.choose()
df <- read_tsv(
  file,
  col_types = list(col_double(), col_double(), col_factor())
)
df <- df %>% dplyr::select(Angle, Width, Species)


############################
# x. normality diagnostics #
############################

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
    results[["all"]] <- boxCox(df %>% dplyr::select(-factors), lambdas, title)
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
  X <- as.matrix(df %>% dplyr::select(-factors))
  X.c <- apply(X, 2, function(x) x - mean(x))
  n <- dim(X)[1]
  p <- dim(X)[2]
  d <- diag(X.c %*% solve(cov(X.c)) %*% t(X.c))
  threshold <- qchisq(1 - alpha, p)
  df <- df %>% mutate(
    distance = d,
    greq_chisq = ifelse(distance > qchisq(1-alpha, p), distance, as.numeric(NA))
  )
  df$greq_chisq[c(which(!is.na(df$greq_chisq)))] = c(which(!is.na(df$greq_chisq)))
  df[[factors]] <- df[[factors]]
  df <- df %>% arrange(distance) %>% mutate(quantile = qchisq((1:n - 0.5)/n, p)) 
  plot <- ggplot(data = df, aes_string(colour = factors, alpha = opacity)) + 
    geom_point(aes(x = quantile,  y = distance)) +
    geom_abline(slope = 1, intercept = 0, type = "dotted") +
    coord_equal() + 
    xlab("chi-square quantile") + 
    ylab("generalized squared distance") +
    geom_hline(yintercept = threshold, linetype = "dashed") +
    geom_text(aes(label = greq_chisq, x = quantile,  y = distance), hjust = -0.3, na.rm = TRUE, show.legend = FALSE)
  return(
    list(
      plot = plot, 
      df = df, 
      specs = list(
        alpha = alpha,
        n = n,
        p = p,
        threshold = threshold
      )
    )
  )
}


diagnosticsUnivariateNormal <- function(
  df,
  factors
){
  
  diagnostics <- function(df, level, vars){
    diagnostics <- list()
    for (var in vars){
      diagnostics[[var]][["histogram"]] <- ggplot(data = df, mapping = aes_string(x = var)) + geom_histogram(bins = 10) + ggtitle(as.character(level))
      diagnostics[[var]][["QQ"]] <- ggplot(data = df, mapping = aes_string(sample = var)) + geom_qq() + geom_qq_line() + ggtitle(as.character(level))
      diagnostics[[var]][["shapiroWilk"]] <- shapiro.test(df[[var]])
      diagnostics[[var]][["correlation"]] <- ppccTest(df[[var]])
    }
    return(diagnostics)
  }
  
  vars <- colnames(df)[which(!(colnames(df)) %in% factors)] 
  results <- list()
  if (factors %in% colnames(df)){
    results[["all"]] <- diagnostics(df %>% dplyr::select(-factors), "all groups", vars)
    for (level in unique(df[[factors]])){
      results[[paste(as.character(level))]] <- diagnostics(filterCategory(df, factors, level), level, vars)
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
  df.tmp <- df %>% dplyr::select(-factors)
  Z <- as.tibble(apply(df.tmp, 2, function(x) (x - mean(x)) / sd(x))) %>% 
    dplyr::filter_all(any_vars(abs(.) > threshold))
  Z <- list("all" = Z)
  for (level in unique(df[[factors]])){
    df.tmp <- filterCategory(df, factors = factors, level = level)
    Z[[level]] <- as.tibble(apply(df.tmp, 2, function(x) (x - mean(x)) / sd(x))) %>% 
      dplyr::filter_all(any_vars(abs(.) > threshold))
  }
  return(Z)
}


filterCategory <- function(
  df, 
  factors, 
  level
){
  expression <- paste(factors, '==', '"', level, '"', sep ='')
  df.tmp <- df %>% 
    dplyr::filter(rlang::eval_tidy(rlang::parse_expr(expression))) %>% 
    dplyr::select(-factors)
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


diagnosticsUnivariateNormal(df, "Species")
diagnosticsZ(df, "Species", threshold = 3)
ggpairs(df, aes(colour = Species, alpha = 0.5), legend = 1)
box_cox <-diagnosticsBoxCox(df, "Species")
plotBoxCox(box_cox) # maybe transform Angle, V2?
diagnosticsQQChiSq(df, "Species", alpha = 0.01) #overspread in center, underspread in outer region

df.transformed <- df %>% mutate(Angle = Angle^2.6)
diagnosticsUnivariateNormal(df.transformed, "Species")
ggpairs(df.transformed, aes(colour = Species, alpha = 0.5))
box_cox <-diagnosticsBoxCox(df.transformed, "Species")
plotBoxCox(box_cox)
diagnosticsQQChiSq(df.transformed, "Species", alpha = 0.01) #overspread in center, underspread in outer region

#rstatix::box_m(df[,1:2], as.vector(as.matrix(df[,3])))
heplots::boxM(cbind(Angle, Width) ~ Species, data = df)
bartlett.test(Angle + Width ~ Species, data = df)
leveneTest(Angle + Width ~ Species, data = df)

#rstatix::box_m(df.transformed[,1:2], as.vector(as.matrix(df.transformed[,3])))
heplots::boxM(cbind(Angle, Width) ~ Species, data = df.transformed)
bartlett.test(Angle + Width ~ Species, data = df.transformed)
leveneTest(Angle + Width ~ Species, data = df.transformed)


###############
# a-b. MANOVA #
###############

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
  .X <- df %>% dplyr::select(-factors)
  X <- list()
  Residuals <- list()
  Treatments <- list()
  Xbar <- colMeans(.X)
  xbar <- list()
  for (level in unique(df[[factors]])){
    X[[level]] <- as.matrix(filterCategory(df, factors, level))
    n[[level]] <- dim(X[[level]])[1]
    S[[level]] <- cov(X[[level]])
    xbar[[level]] <- t(colMeans(filterCategory(df, factors, level)))
    Treatments[[level]] <- Xbar - xbar[[level]]
    Residuals[[level]] <- t(apply(X[[level]], 1, function(x) x - xbar[[level]]))
  }
  p <- dim(X[[1]])[2]
  g <- length(X)
  Residuals <- as.tibble(do.call(rbind, Residuals))
  colnames(Residuals) <- colnames(.X)
  Residuals <- list(
    Residuals = Residuals, 
    Standardized = as.tibble(apply(Residuals, 2, function(x) (x-mean(x))/sd(x)))
  )
  S[["S"]]  <- cov(df %>% dplyr::select(-factors))
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
  return(
    list(
      Residuals = Residuals,
      Treatments = Treatments,
      T = T, 
      B = B, 
      W = W,
      Xbar = Xbar, 
      xbar = xbar, 
      X = X,
      S = S,
      df = df,
      pillai = pillai,
      wilks = wilks,
      specs = specs
    )
  )
}


manova <- MANOVA(df, "Species")
manova
qplot(Angle, Width, data = manova$Residuals$Standardized, colour = df$Species) + labs(colour = "Species")

# Same in cars!
lm <- lm(cbind(Width, Angle) ~ Species, data = df, contrasts = list(Species=contr.sum))
manova.cars <- Manova(lm)
summary(manova.cars)

############################
# c. Compare Mu_1 and Mu_2 #
############################

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
  delta <- manova$xbar[[i]] - manova$xbar[[j]] - delta_0
  
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


df.red <- df %>% dplyr::filter(Species != "Hep")
df.red$Species <- droplevels(df.red$Species)
manova <- MANOVA(df.red, "Species")
rstatix::box_m(df.red[,1:2], as.vector(as.matrix(df.red[,3])))
heplots::boxM(cbind(Angle, Width) ~ Species, data = df.red)
bartlett.test(Angle + Width ~ Species, data = df.red)
leveneTest(Angle + Width ~ Species, data = df.red)
manova
comparison <- compareMeanVectors(manova, 1, 2, confint_method = c("b", "h"))
comparison

#########################
# d. Confidence Ellipse #
#########################

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
plot <- ggplot(data = NULL) + xlab("Angle") + ylab("Width")
plot <- plotConfidenceEllipse(ellipse, plot)
plot <- plotConfidenceIntervals(comparison, plot)
plot
qplot(Angle, Width, data = df.red, colour = Species)

####################################
# e. Classification and Clustering #
####################################

calculateAPER <- function(confusion_matrix){
    correct <- sum(diag(confusion_matrix))
    N <- sum(confusion_matrix)
    return((N - correct)/N)
}


discriminantScore <- function(
  obs,
  p,
  S,
  xbar,
  costs,
  distribution = "mvnormal",
  method = "linear"
){
  if (method == "linear"){score <- t(xbar) %*% solve(S) %*% obs - 0.5 * t(xbar) %*% solve(S) %*% xbar + log(p)}
  if (method == "quadratic"){score <- -0.5 * log(det(S)) - 0.5 * t(obs) %*% solve(S) %*% obs + log(p)}
  return(score)
}


discriminantScores <- function(
  df,
  classes,
  costs,
  equal_priors = TRUE,
  equal_variance = TRUE,
  levels = NULL,
  method = "linear",
  priors = NULL
){
  if (classes %in% colnames(df)){
    .X <- as.matrix(df %>% dplyr::select(-classes))
    if (is.null(levels)){levels <- levels(df[[classes]])}
    if (is.null(priors)){priors <- generatePriors(df, classes, equal_priors = equal_priors, equal_variance = equal_variance)}
    scores <- tibble(classes = df[[classes]])
  } else{
    .X <- as.matrix(df)
    scores <- tibble(rep(NA, dim(.X)[1]))
  }
  n <- dim(.X)[1]
  k <- length(levels)
  p <- dim(.X)[2]
  
  for (i in 1:k){
    posteriors <- c()
    for (j in 1:n){
      delta <- .X[j,] - priors$xbars[[i]]
      if (method == "linear"){.score <- discriminantScore(t(t(.X[j,])), priors$p[[i]], priors$S, t(t(priors$xbars[[i]])), priors$costs, distribution, method)}
      if (method == "quadratic"){.score <- discriminantScore(t(t(delta)), priors$p[[i]], priors$S[[i]], t(t(priors$xbars[[i]])), priors$costs, distribution, method)}
      posteriors <- c(posteriors, .score)
    }
    scores <- scores %>% add_column(!! paste(levels[i], '_score', sep='') := posteriors)
  }
  new_classes <- c()
  for (i in 1:n){
    level <- order(scores[i, -1], decreasing = TRUE)[1]
    new_classes <- c(new_classes, levels[level])
  }

  scores$new_classes <- as.factor(new_classes)
  return(
    list(
      old_classes = df[[classes]],
      new_classes = scores$new_classes,
      df = df,
      scores = scores
    )
  )
}


expectationMaximization <- function(
  df,
  classes,
  random_seeds = FALSE,
  epsilon = 0.001,
  equal_priors = TRUE,
  equal_variance = FALSE,
  knn = FALSE,
  naive = FALSE,
  priors = NULL
){
  .X <- as.matrix(df %>% dplyr::select(-classes))
  levels <- levels(df[[classes]])
  if (is.null(priors)){priors <- generatePriors(df, classes, equal_priors = equal_priors, equal_variance = equal_variance)}
  n <- dim(.X)[1]
  k <- length(levels)
  p <- dim(.X)[2]
  if (naive && !knn) {
    for (i in 1:k){
      priors$S[[i]] <- diag(1, p)
    }
  }
  if (random_seeds) {
    for (i in 1:k){
      mu <- c()
      for (j in 1:p){
        mu <- c(mu, runif(1, range(df[,j])[1], range(df[,j])[2]))
      }
      priors$xbars[[i]] <- mu
    }
  }
  
  delta <- NULL
  iterations <- 0
  while (delta > epsilon || is.null(delta)){
    scores <- tibble(classes = df[[classes]])
    for (i in 1:k){
      posteriors <- c()
      if (knn) {
        .Xc <- t(apply(.X, 1, function(x) x - priors$xbars[[i]]))
        posteriors <- diag(.Xc %*% solve(priors$S[[i]]) %*% t(.Xc))
      } else{
        for (j in 1:n){posteriors <- c(posteriors, priors$p[[i]]*dmvnorm(.X[j,], priors$xbars[[i]], priors$S[[i]]))}
      }
      scores <- scores %>% add_column(!! paste(levels[i], '_score', sep='') := posteriors)
    }
    new_classes <- c()
    for (i in 1:n){
      row_scores <- scores[i, 2:(1+k)]/sum(scores[i, 2:(1+k)])
      new_classes <- c(new_classes, order(row_scores, decreasing = ifelse(knn, FALSE, TRUE))[1])
    }
    df$new_classes <- as.factor(levels[new_classes])
    delta <- 0
    old_priors <- priors
    priors <- generatePriors(df %>% dplyr::select(-classes), "new_classes")
    for (i in 1:k){
      delta <- delta + sum((old_priors$xbars[[i]] - priors$xbars[[i]])^2)
    }
    iterations <- iterations + 1
  }
  confusion_matrix <- generateConfusionMatrix(df[[classes]], levels[new_classes])$confusion_matrix
  remapping <- levels[unconfuseMatrix(confusion_matrix)]
  df$new_classes <- as.factor(remapping[new_classes])
  
  return(
    list(
      old_classes = df[[classes]],
      new_classes = df$new_classes,
      df = df,
      iterations = iterations,
      posteriors = priors
    )
  )
}


generateConfusionMatrix <- function(
  actual,
  class
){
  levels <- levels(df$Species)
  k <- length(levels)
  n <- length(actual)
  confusion <- matrix(0, ncol = k, nrow = k)
  cells <- cbind(
    sapply(actual, function(x) match(x, levels)),
    sapply(class, function(x) match(x, levels))
  )
  for (i in 1:n) {
    confusion[cells[i, 1], cells[i,2]] <- confusion[cells[i, 1], cells[i,2]] + 1
  }
  return(
    list(
      APER = calculateAPER(confusion),
      confusion_matrix = confusion
    )
  )
}


generateContourPoints <- function(
  var1, 
  var2,
  extra_margin = 0.0,
  length = 100
){
  
  generateSeq <- function(var){
    range <- range(var)
    scaled_margin <- extra_margin*(range[2] - range[1])
    return(seq(range[1] - scaled_margin, range[2] + scaled_margin, length.out = length))
  }
  
  range1 <- generateSeq(var1)
  range2 <- generateSeq(var2)
  return(as.tibble(expand.grid(range1, range2)))
}


generatePriors <- function(
  df,
  classes,
  equal_priors = FALSE,
  equal_variance = FALSE
){
  levels <- levels(df[[classes]])
  k <- length(levels)
  n <- list()
  p <- dim(df)[2] - 1
  P <- list()
  S <- list()
  xbars <- list()
  for (i in 1:k){
    df.tmp <- filterCategory(df, classes, levels[[i]])
    n[[i]] <- dim(df.tmp)[1]
    P[[i]] <- ifelse(equal_priors, 1/k, dim(df.tmp)[1]/dim(df)[1])
    xbars[[i]] <- colMeans(df.tmp)
    S[[i]] <- cov(df.tmp)
  }
  if (equal_variance){
    Spool <- matrix(data = 0, ncol = p, nrow = p) 
    for (i in 1:k){Spool <- Spool + (n[[i]] - 1) * S[[i]]}
    S <- Spool / (sum(unlist(n))-k)
  } 
  
  priors = list(p = P, xbars = xbars, S = S)
  return(priors)
}


unconfuseMatrix <- function(confusion_matrix){
  k <- dim(confusion_matrix)[2]
  remapping <- c()
  for (i in 1:k){
    order <- order(confusion_matrix[, i], decreasing = TRUE)
    remapping <- c(remapping, order[1])
  }
  return(remapping)
}


validationKFold <- function(
  df,
  classes,
  job = list(
    equal_priors = NULL,
    equal_variance = NULL,
    expression = NULL,
    k = NULL, 
    reps = NULL
  )
){
  aper <- 0
  g <- levels(df[[classes]])
  overall_confusion_matrix <- matrix(0, ncol = length(g), nrow = length(g))
  for (i in 1:job$reps){
    if (job$k == 0){
      folds <- list(1:dim(df)[1])
    } else{
      folds <- createFolds(df[[classes]], k = job$k, list = TRUE, returnTrain = FALSE)
    }
    for (fold in folds){
      if (job$k == 0){
        priors <- generatePriors(df, classes, equal_priors = job$equal_priors, equal_variance = job$equal_variance)
        df.tmp <- df
      } else{
        priors <- generatePriors(df[-fold, ], classes, equal_priors = job$equal_priors, equal_variance = job$equal_variance)
        df.tmp <- df[fold,]
      }
      .model <- base::eval(rlang::parse_expr(job$expression))
      .cm <- generateConfusionMatrix(.model$old_classes, .model$new_classes)
      aper <- aper + .cm$APER
      overall_confusion_matrix <- overall_confusion_matrix + .cm$confusion_matrix
    }
  }
  colnames(overall_confusion_matrix) <- g
  rownames(overall_confusion_matrix) <- g
  if(job$k == 0){
    iterations <- job$reps
  } else{
    iterations <- (job$k * job$reps)
  }
  
  return(
    list(
      job = job,
      overall_APER = aper / iterations,
      overall_confusion_matrix = overall_confusion_matrix / iterations
    )
  )
}


validation <- function(
  df,
  jobs = list(
    list(
      expression = NULL,
      attempts = 10,
      equal_priors = TRUE,
      equal_variance = TRUE,
      k = 5,
      reps = 1
    )
  )
){
  results <- list()
  for (job in jobs){
    print("Beginning validation job:")
    print(job$expression)
    attempt <- 0
    done <- FALSE
    while (attempt < job$attempts && !done){
      attempt <- attempt + 1
      print(paste("On attempt: ", as.character(attempt), sep=''))
      result <- NULL
      try(result <- validationKFold(df, "Species", job))
      if (!is.null(result)){
        results[[job$expression]] <- result
        done <- TRUE
      }
    }
  }
  
  return(results)
}


jobs.clustering.knn <- list(
  list(
    expression = 'expectationMaximization(df.tmp, "Species", equal_priors = TRUE, knn = TRUE, priors = priors)', 
    attempts = 20, 
    equal_priors = TRUE, 
    equal_variance = FALSE, 
    k = 5, 
    reps = 10
  ),
  list(
    expression = 'expectationMaximization(df.tmp, "Species", equal_priors = TRUE, knn = TRUE, priors = priors, random_seeds = TRUE)', 
    attempts = 1000, 
    equal_priors = TRUE, 
    equal_variance = FALSE, 
    k = 5, 
    reps = 1
  ),
  list(
    expression = 'expectationMaximization(df.tmp, "Species", equal_priors = FALSE, knn = TRUE, priors = priors)',
    attempts = 20,
    equal_priors = TRUE,
    equal_variance = FALSE,
    k = 5,
    reps = 10
  ),
  list(
    expression = 'expectationMaximization(df.tmp, "Species", equal_priors = FALSE, knn = TRUE, priors = priors, random_seeds = TRUE)',
    attempts = 1000,
    equal_priors = TRUE,
    equal_variance = FALSE,
    k = 5,
    reps = 1
  )
)
results.clustering.knn <- validation(df, jobs = jobs.clustering.knn)
results.clustering.knn

jobs.clustering.EM <- list(
  list(
    expression = 'expectationMaximization(df.tmp, "Species", equal_priors = TRUE, naive = FALSE, priors = priors)', 
    attempts = 100, 
    equal_priors = TRUE, 
    equal_variance = FALSE, 
    k = 5, 
    reps = 10
  ),
  list(
    expression = 'expectationMaximization(df.tmp, "Species", equal_priors = FALSE, naive = FALSE, priors = priors)',
    attempts = 100, 
    equal_priors = TRUE, 
    equal_variance = FALSE, 
    k = 5, 
    reps = 10
  ),
  list(
    expression = 'expectationMaximization(df.tmp, "Species", equal_priors = TRUE, naive = FALSE, priors = priors, random_seeds = TRUE)',
    attempts = 1000, 
    equal_priors = TRUE, 
    equal_variance = FALSE, 
    k = 5, 
    reps = 1
  )
)
results.clustering.EM <- validation(df, jobs = jobs.clustering.EM)
results.clustering.EM

jobs.clustering.EMnaive <- list(
  list(
    expression = 'expectationMaximization(df.tmp, "Species", equal_priors = TRUE, naive = TRUE, priors = priors)',
    attempts = 100, 
    equal_priors = TRUE, 
    equal_variance = FALSE, 
    k = 5, 
    reps = 10
  ),
  list(
    expression = 'expectationMaximization(df.tmp, "Species", equal_priors = FALSE, naive = TRUE, priors = priors)',
    attempts = 100, 
    equal_priors = TRUE, 
    equal_variance = FALSE, 
    k = 5, 
    reps = 10
  ),
  list(
    expression = 'expectationMaximization(df.tmp, "Species", equal_priors = TRUE, naive = TRUE, priors = priors, random_seeds = TRUE)',
    attempts = 1000, 
    equal_priors = TRUE, 
    equal_variance = FALSE, 
    k = 5, 
    reps = 1
  )
)
results.clustering.EMnaive <- validation(df, jobs = jobs.clustering.EMnaive)
results.clustering.EMnaive

jobs.classification.nofolds <- list(
  list(
    expression = 'discriminantScores(df.tmp, "Species", NULL, equal_priors = TRUE, method = "linear", priors = priors)',
    attempts = 1, 
    equal_priors = TRUE, 
    equal_variance = TRUE, 
    k = 0, 
    reps = 1
  ),
  list(
    expression = 'discriminantScores(df.tmp, "Species", NULL, equal_priors = TRUE, equal_variance = FALSE, method = "quadratic", priors = priors)',
    attempts = 1, 
    equal_priors = TRUE, 
    equal_variance = FALSE, 
    k = 0, 
    reps = 1
  )
)
results.classification.nofolds <- validation(df, jobs = jobs.classification.nofolds)
results.classification.nofolds

jobs.classification.OHAAT <- list(
  list(
    expression = 'discriminantScores(df.tmp, "Species", NULL, equal_priors = TRUE, method = "linear", priors = priors)',
    attempts = 1, 
    equal_priors = TRUE, 
    equal_variance = TRUE, 
    k = 74, 
    reps = 1
  ),
  list(
    expression = 'discriminantScores(df.tmp, "Species", NULL, equal_priors = TRUE, equal_variance = FALSE, method = "quadratic", priors = priors)',
    attempts = 1, 
    equal_priors = TRUE, 
    equal_variance = FALSE, 
    k = 74, 
    reps = 1
  )
)
results.classification.OHAAT <- validation(df, jobs = jobs.classification.OHAAT)
results.classification.OHAAT

jobs.classification.2folds <- list(
  list(
    expression = 'discriminantScores(df.tmp, "Species", NULL, equal_priors = TRUE, method = "linear", priors = priors)',
    attempts = 1, 
    equal_priors = TRUE, 
    equal_variance = TRUE, 
    k = 2, 
    reps = 1
  ),
  list(
    expression = 'discriminantScores(df.tmp, "Species", NULL, equal_priors = TRUE, equal_variance = FALSE, method = "quadratic", priors = priors)',
    attempts = 1, 
    equal_priors = TRUE, 
    equal_variance = FALSE, 
    k = 2, 
    reps = 1
  )
)
results.classification.2folds <- validation(df, jobs = jobs.classification.5folds)
results.classification.2folds

jobs.classification.5folds <- list(
  list(
    expression = 'discriminantScores(df.tmp, "Species", NULL, equal_priors = TRUE, method = "linear", priors = priors)',
    attempts = 1, 
    equal_priors = TRUE, 
    equal_variance = TRUE, 
    k = 5, 
    reps = 100
  ),
  list(
    expression = 'discriminantScores(df.tmp, "Species", NULL, equal_priors = TRUE, equal_variance = FALSE, method = "quadratic", priors = priors)',
    attempts = 1, 
    equal_priors = TRUE, 
    equal_variance = FALSE, 
    k = 5, 
    reps = 100
  )
)
results.classification.5folds <- validation(df, jobs = jobs.classification.5folds)
results.classification.5folds

jobs.classification.10folds <- list(
  list(
    expression = 'discriminantScores(df.tmp, "Species", NULL, equal_priors = TRUE, method = "linear", priors = priors)',
    attempts = 1, 
    equal_priors = TRUE, 
    equal_variance = TRUE, 
    k = 10, 
    reps = 100
  ),
  list(
    expression = 'discriminantScores(df.tmp, "Species", NULL, equal_priors = TRUE, equal_variance = FALSE, method = "quadratic", priors = priors)',
    attempts = 1, 
    equal_priors = TRUE, 
    equal_variance = FALSE, 
    k = 10, 
    reps = 100
  )
)
results.classification.10folds <- validation(df, jobs = jobs.classification.5folds)
results.classification.10folds

contour <- discriminantScores(contour_points, "Species", 0, levels = unique(df$Species), method = "linear", priors = generatePriors(df, "Species", equal_variance = TRUE))$new_classes
contour <- tibble(contour_points, z = as.numeric(contour))
ggplot(df, aes(x = Angle, y = Width)) +
  geom_point(aes(colour = Species)) + 
  geom_raster(data = contour, aes(Var1, Var2, fill = z), alpha = 0.3, show.legend = FALSE, interpolate = TRUE)

contour_points <- generateContourPoints(df$Angle, df$Width)
contour <- discriminantScores(contour_points, "Species", 0, levels = unique(df$Species), method = "quadratic", priors = generatePriors(df, "Species"))$new_classes
contour <- tibble(contour_points, z = as.numeric(contour))
ggplot(df, aes(x = Angle, y = Width)) +
  geom_point(aes(colour = Species)) + 
  geom_raster(data = contour, aes(Var1, Var2, fill = z), alpha = 0.3, show.legend = FALSE, interpolate = TRUE)
  
