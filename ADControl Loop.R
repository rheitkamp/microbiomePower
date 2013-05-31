###This loop is not yet complete

library("PearsonDS")

#####################Backbone of AD loop#######################
x=10
y=5
ADControl <- matrix(nrow=x,ncol=y)
ADControl <- data.frame(ADControl)
colnames(ADControl) <- c("Firmicutes", "Actinobacteria", "Proteobacteria", "Bacteroidetes", "Others")
size <- dim(ADControl)[1]

for (i in 1:size){
  total <- 100
  ###Generate number for Firmicutes
  ADControl$Firmicutes[i] <- rpearsonI(n=1, a=0.9970149, b=2.0062602, location=10.6318818, scale=89.3067879)
  total <- total - ADControl$Firmicutes[i]
  
  ###Generate number for Actinobacteria
  rActino <- rpearsonII(n=1, a=0.840323, location=9.764036, scale=81.726210)
  ADControl$Actinobacteria[i] <- (total*rActino)/100
  total <- total - ADControl$Actinobacteria[i]
  
  ###Generate number for Proteobacteria
  rProteo <- rpearsonI(n=1, a=1.371166,  b=1.477854, location=45.429029, scale=41.340968)
  ADControl$Proteobacteria[i] <- (total*rProteo)/100
  total <- total - ADControl$Proteobacteria[i]
  
  ###Generate number for Bacteridetes
  rBacter <- rpearsonI(n=1, a=2.2098166, b=0.6798847, location=-12.4973479, scale=110.4140148)
  ADControl$Bacteroidetes[i] <- (total*rBacter)/100
  ADControl$Others[i] <- total - ADControl$Bacteroidetes[i]
  
}

ADControl <- data.matrix(ADControl)