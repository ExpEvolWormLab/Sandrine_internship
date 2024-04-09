library(poolSeq)
library(ggplot2)
library(gridExtra)

###################################################################################################################
#######  Ne over generation for Androdioecy population ################################################
##################################################################################################################


#Generation

A_generation <-c(0, 100,36,50,68, 100,10,36,50,68, 100,10,36,50,68, 100,10,32,66, 100,10,32,66,10,32)
#Replicates
A_replicates<- c(0,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,5,5,5,5,6,6)

#Poolsize
Ne <- 10000

#Read sync file
CA.sync<-read.sync("CA_correct.sync",gen=A_generation,repl=A_replicates)

#extract alleles frequency
Allele.frequency= af(CA.sync,gen=A_generation,repl=A_replicates)

#extract coverage
Allele.coverage = coverage(CA.sync,gen=A_generation,repl=A_replicates)

#extract alleles frequency and coverage for generation 0
Allele.freq0= af(CA.sync,gen=0,repl=0)
Allele.coverage0 = coverage(CA.sync,gen=0,repl=0)

Ne_A <- numeric(length(A_generation) - 1) 

for(i in 2:length(A_generation)) {
  Ne_A[i-1] <- estimateNe(p0=Allele.freq0, pt=Allele.frequency[,i], cov0=Allele.coverage0, covt=Allele.coverage[,i], method="W.planI", t=A_generation[i], Ncensus=Ne)
}



generation=c("100"," 36"," 50"," 68","100"," 10"," 36"," 50"," 68","100"," 10"," 36"," 50"," 68","100"," 10"," 32"," 66","100"," 10"," 32"," 66"," 10"," 32")

replicate= c("A1","A1","A1","A1","A2","A2","A2","A2","A2","A3","A3","A3","A3","A3","A4","A4","A4","A4","A5","A5","A5","A5","A6","A6")

population=data.frame(replicate=replicate,generation=generation,Ne=Ne_A )

pdf(file = "CA_Ne.pdf")

ggplot(population, aes(x=generation, y=Ne,group = replicate)) + 
  geom_line(aes(colour = replicate)) + geom_point(aes(colour = replicate))

dev.off()


###################################################################################################################
#######  Ne over generation for Dioecy population ################################################
##################################################################################################################



#Generation
D_generation <-c(35,50,100,35,100,35,50,100,35,50,100,32,100,32,100,32,100,0)
#Replicates
D_replicates<- c(1,1,1,2,2,3,3,3,4,4,4,5,5,6,6,7,7,0)


#Read sync file

CD.sync<-read.sync("CD.sync",gen=D_generation,repl=D_replicates)

#extract alleles frequency
D_Allele.frequency= af(CD.sync,gen=D_generation,repl=D_replicates)

#extract coverage
D_Allele.coverage = coverage(CD.sync,gen=D_generation,repl=D_replicates)

#extract alleles frequency and coverage for generation 0
D_Allele.freq0= af(CD.sync,gen=0,repl=0)
D_Allele.coverage0 = coverage(CD.sync,gen=0,repl=0)

Ne_D <- numeric(length(D_generation) - 1) 

for(i in 2:length(D_generation)) {
  Ne_D[i-1] <- estimateNe(p0=D_Allele.freq0, pt=D_Allele.frequency[,i], cov0=D_Allele.coverage0, covt=D_Allele.coverage[,i], method="W.planI", t=D_generation[i], Ncensus=Ne)
}


Generation = c(" 35"," 50","100"," 35","100"," 35"," 50","100"," 35"," 50","100"," 32","100"," 32","100"," 32","100")
Replicate = c("D1","D1","D1","D2","D2","D3","D3","D3","D4","D4","D4","D5","D5","D6","D6","D7","D7")

D_population=data.frame(replicate=Replicate,generation=Generation,Ne=Ne_D )

pdf(file = "CD_Ne.pdf")

ggplot(D_population, aes(x=generation, y=Ne,group = replicate)) + 
  geom_line(aes(colour = replicate)) + geom_point(aes(colour = replicate))

dev.off()


####################################################################################################################
#######  Average Ne over generation for CD and CA ################################################
##################################################################################################################

# Function to calculate mean and confidence interval
calc_mean_ci <- function(data_vector) {
  mean_value <- mean(data_vector)
  ci <- t.test(data_vector)$conf.int
  return(c(mean_value, ci[1], ci[2]))
}


Dgeneration=c(" 32"," 35"," 50","100")

#combine Ne value per generation

D<- list(
  D32 = c(Ne_D[12],Ne_D[14],Ne_D[16]),
  D35 = c(Ne_D[1],Ne_D[4],Ne_D[6],Ne_D[9]),
  D50 = c(Ne_D[2],Ne_D[7],Ne_D[10]),
  D100 = c(Ne_D[3],Ne_D[5],Ne_D[8],Ne_D[11],Ne_D[13],Ne_D[15],Ne_D[17])
)

CD <- setNames(D,Dgeneration)

# Calculate mean and confidence interval for each vector
mean_ci_D  <- lapply(CD, calc_mean_ci)

# Create a data frame for plotting
plot_data_D <- data.frame(
  vector = rep(names(CD), each = 3),
  stat = rep(c("Mean", "CI_Lower", "CI_Upper"), length(CD)),
  value = unlist(mean_ci_D )
)


P_D <- ggplot(plot_data_D, aes(x = vector, y = value)) +
  geom_point(data = subset(plot_data_D, stat == "Mean"), size = 3) +
  geom_errorbar(aes(ymin = ifelse(stat == "CI_Lower", value, NA), 
                    ymax = ifelse(stat == "CI_Upper", value, NA)),
                width = 0.2, color = "black") +
  geom_line() +
  labs(x = "Generation", y = "Ne")


Ageneration=c(" 10"," 32"," 36"," 50"," 66"," 68","100")

#combine Ne value per generation

A<- list(
  A100 = c(Ne_A[1],Ne_A[5],Ne_A[10],Ne_A[15],Ne_A[19]),
  A10 = c(Ne_A[6],Ne_A[11],Ne_A[16],Ne_A[20],Ne_A[23]),
  A32 = c(Ne_A[17],Ne_A[21],Ne_A[24]),
  A36 = c(Ne_A[2],Ne_A[7],Ne_A[12]),
  A50 = c(Ne_A[3],Ne_A[8],Ne_A[13]),
  A66 = c(Ne_A[18],Ne_A[22]),
  A68 = c(Ne_A[4],Ne_A[9],Ne_A[14])
)

CA <- setNames(A, Ageneration)


# Calculate mean and confidence interval for each generation
mean_ci_A <- lapply(CA, calc_mean_ci)

# Create a data frame for plotting
plot_data_A <- data.frame(
  vector = rep(names(CA), each = 3),
  stat = rep(c("Mean", "CI_Lower", "CI_Upper"), length(CA)),
  value = unlist(mean_ci_A)
)


P_A <- ggplot(plot_data_A, aes(x = vector, y = value)) +
  geom_point(data = subset(plot_data_A, stat == "Mean"), size = 2, color = "red") +
  geom_errorbar(aes(ymin = ifelse(stat == "CI_Lower", value, NA), 
                    ymax = ifelse(stat == "CI_Upper", value, NA)),
                width = 0.2, color = "black") +
  geom_line() +
  labs(x = "Generation", y = "Ne")

pdf(file = "Ne_AvsD.pdf")
grid.arrange(P_A, P_D, ncol = 2)
dev.off()
