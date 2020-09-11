library(tidyverse)

# 2.27 
A <- matrix(c(1, 2, 2, 1, -2, 2), nrow=3)
AAt <- A %*% t(A)02.R
eigen(AAt)

# 2.38
A <- matrix(c(13, -4, 2, -4, 13, -2, 2, -2, 10), nrow=3)
eigen(A)

# 1.2
x1 <- c(1, 2, 3, 3, 4, 5, 6, 8, 9, 11)
x2 <- c(18.95, 19.00, 17.95, 15.54, 14.00, 12.95, 8.94, 7.49, 6.00, 3.99)
qplot(x1, x2)
mean(x1)
mean(x2)
var(x1)
var(x2)
cov(x1, x2)
cor(x1, x2)

# 1.17
file <- file.choose()
DATA <- read_tsv(file, col_names = FALSE)
DATA %>% select(-X1) %>% summarise_all(funs(mean)) %>% as.matrix()
cov(DATA[,2:8])
cor(DATA[,2:8])

# 1.18
X2 <- c(100/DATA[,2])
X3 <- c(200/DATA[,3])
X4 <- c(400/DATA[,4])
X5 <- c(800/(DATA[,5]*60))
X6 <- c(1500/(DATA[,6]*60))
X7 <- c(3000/(DATA[,7]*60))
X8 <- c(42195/(DATA[,8]*60))
DATA.adj <- as.tibble(cbind(DATA[,1], X2, X3, X4, X5, X6, X7, X8))
DATA.adj %>% select(-X1) %>% summarise_all(funs(mean)) %>% as.matrix()
cov(DATA.adj[,2:8])
cor(DATA.adj[,2:8])

# 2.22
A = as.matrix(rbind(c(4, 8, 8), c(3, 6, -9)))
AAt <- eigen(A %*% t(A))
AtA <- eigen(t(A) %*% A)           
AAt$vectors[,1] %*% t(AtA$vectors[,1])
AAt$vectors[,2] %*% t(AtA$vectors[,2])
AAt
AtA
sqrt(AAt$values[1])*-1 * AAt$vectors[,1] %*% t(AtA$vectors[,1]) + sqrt(AAt$values[2]) * AAt$vectors[,2] %*% t(AtA$vectors[,2])
