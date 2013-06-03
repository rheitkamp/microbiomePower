#enter the parameters as reported in the paper
Proteobacteria <-c(75.39, 23.81)
Firmicutes <- c(17.22, 16.58)
Bacteroidetes <- c(3.39, 8.99)
Actinobacteria <- c(1.88, 3.36)
Fusobacteria <- c(1.75, 8.14)
Cyanobacteria <- c(0.27, 0.44)
OD1 <- c(0.067, 0.086)
TM7 <- c(0.045, 0.305)
DeinococcusThermus <- c(0.015, 0.030)
Nitrospira <- c(0.005, 0.012)
Planctomycetes <- c(0.001, 0.008)
Chloroflexi <- c(0.001, 0.005)
BRC1 <- c(0.000, 0.000)

Prot1<- as.numeric(Proteobacteria)

getBetaParams <- function(mean, sd) {
  m <- (1-mean)/mean
  n <- 1 + m
  alpha <- (1/n)(m/(sd^2*n^2)-1)
  beta <- m * alpha
  params <- c(alpha, beta)
  return(params)
}

getBetaParams(mean = 75.39, sd = 23.81)
