###This loop is not yet complete

library("PearsonDS")

###AD Controls %remaining dataset
control <- read.table("ADcontrolsredist.txt")

###fitted distribution taxa generation 0-100
Firmicutes1 <- rpearsonI(n=1, a=0.9970149, b=2.0062602, location=10.6318818, scale=89.3067879)
Actinobacteria1 <- rpearsonII(n=1, a=0.840323, location=9.764036, scale=81.726210)
Proteobacteria1 <- rpearsonI(n=1, a=1.371166,  b=1.477854, location=45.429029, scale=41.340968)
Bacteroidetes1 <- rpearsonI(n=1, a=2.2098166, b=0.6798847, location=-12.4973479, scale=110.4140148)


#####################Backbone of AD loop#######################

total <- 100
###Generate number for Firmicutes
Firmicutes <- rpearsonI(n=1, a=0.9970149, b=2.0062602, location=10.6318818, scale=89.3067879)
total <- total - Firmicutes

###Generate number for Actinobacteria
rActino <- rpearsonII(n=1, a=0.840323, location=9.764036, scale=81.726210)
Actinobacteria <- (total*rActino)/100
total <- total - Actinobacteria

###Generate number for Proteobacteria
rProteo <- rpearsonI(n=1, a=1.371166,  b=1.477854, location=45.429029, scale=41.340968)
Proteobacteria <- (total*rProteo)/100
total <- total - Proteobacteria

###Generate number for Bacteridetes
rBacter <- rpearsonI(n=1, a=2.2098166, b=0.6798847, location=-12.4973479, scale=110.4140148)
Bacteroidetes <- (total*rBacter)/100
Others <- total - Bacteroidetes

