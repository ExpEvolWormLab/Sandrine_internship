library(poolSeq)
library(ggplot2)
library(gridExtra)

#Generation
D_generation <-c(35,50,100,35,100,35,50,100,35,50,100,32,100,32,100,32,100,0)
#Replicates
D_replicates<- c(1,1,1,2,2,3,3,3,4,4,4,5,5,6,6,7,7,0)

#Read sync file

CD.sync<-read.sync("CD.sync",gen=D_generation,repl=D_replicates)

#extract alleles frequency and coverage for generation 0
Allele.freq0= af(CD.sync,gen=0,repl=D_replicates)
Allele.coverage0 = coverage(CD.sync,gen=0,repl=D_replicates)

D0 <- Allele.freq0*Allele.coverage0
d0 <- Allele.coverage0 - D0


#extract alleles frequency and coverage for generation 32
Allele.freq32= af(CD.sync,gen=32,repl=D_replicates)

Allele.coverage32 = coverage(CD.sync,gen=32,repl=D_replicates)

D32<-Allele.freq32*Allele.coverage32
d32<-Allele.coverage32 - D32

Chr_32 <- gsub("\\..*", "", rownames(Allele.freq32))

# perform CMH-test for all empirical loci
p.values_cmh_32 <- cmh.test(A0=rbind(D0,D0,D0), a0=rbind(d0,d0,d0), At=t(D32), at=t(d32),min.cov=1, max.cov=1, min.cnt=1, log=TRUE)

cmh_32<-data.frame(value =p.values_cmh_32,position= c(1:length(p.values_cmh_32)),Chr = Chr_32 )

p_cmh_32<-ggplot(na.omit(cmh_32), aes(x=position, y=value,group=Chr))+
  geom_point() + geom_point(aes(colour = Chr))+labs(x = "Position ", y = "-log10(p) ",title = "CDG32")


#extract alleles frequency and coverage for generation 35
Allele.freq35= af(CD.sync,gen=35,repl=D_replicates)

Allele.coverage35 = coverage(CD.sync,gen=35,repl=D_replicates)

D35<-Allele.freq35*Allele.coverage35
d35<-Allele.coverage35 - D35

Chr_35 <- gsub("\\..*", "", rownames(Allele.freq35))

# perform CMH-test for all empirical loci
p.values_cmh_35 <- cmh.test(A0=rbind(D0,D0,D0,D0), a0=rbind(d0,d0,d0,d0), At=t(D35), at=t(d35),min.cov=1, max.cov=1, min.cnt=1, log=TRUE)

cmh_35<-data.frame(value =p.values_cmh_35,position= c(1:length(p.values_cmh_35)),Chr = Chr_35)

p_cmh_35<-ggplot(na.omit(cmh_35), aes(x=position/1e4, y=value,group=Chr))+
  geom_point() + geom_point(aes(colour = Chr))+labs(x = "Position ", y = "-log10(p) ",title = "CDG35")


#extract alleles frequency and coverage for generation 50
Allele.freq50= af(CD.sync,gen=50,repl=D_replicates)

Allele.coverage50 = coverage(CD.sync,gen=50,repl=D_replicates)

D50<-Allele.freq50*Allele.coverage50
d50<-Allele.coverage50 - D50

Chr_50 <- gsub("\\..*", "", rownames(Allele.freq50))

# perform CMH-test for all empirical loci
p.values_cmh_50 <- cmh.test(A0=rbind(D0,D0,D0), a0=rbind(d0,d0,d0), At=t(D50), at=t(d50),min.cov=1, max.cov=1, min.cnt=1, log=TRUE)

cmh_50<-data.frame(value =p.values_cmh_50,position= c(1:length(p.values_cmh_50)),Chr = Chr_50 )

p_cmh_50<-ggplot(na.omit(cmh_50), aes(x=position, y=value,group=Chr))+
  geom_point() + geom_point(aes(colour = Chr))+labs(x = "Position ", y = "-log10(p) ",title = "CDG50")


#extract alleles frequency and coverage for generation 100
Allele.freq100= af(CD.sync,gen=100,repl=D_replicates)

Allele.coverage100 = coverage(CD.sync,gen=100,repl=D_replicates)

D100<-Allele.freq100*Allele.coverage100
d100<-Allele.coverage100 - D100

Chr_100 <- gsub("\\..*", "", rownames(Allele.freq100))

# perform CMH-test for all empirical loci
p.values_cmh_100 <- cmh.test(A0=rbind(D0,D0,D0,D0,D0,D0,D0), a0=rbind(d0,d0,d0,d0,d0,d0,d0), At=t(D100), at=t(d100),min.cov=1, max.cov=1, min.cnt=1, log=TRUE)

cmh_100<-data.frame(value =p.values_cmh_100,position= c(1:length(p.values_cmh_100)),Chr = Chr_100)

p_cmh_100_D<-ggplot(na.omit(cmh_100), aes(x=position/1e4, y=value,group=Chr))+
  geom_point() + geom_point(aes(colour = Chr))+labs(x = "Position ", y = "-log10(p) ",title = "CDG100")




N= max(p.values_cmh_100, na.rm = TRUE)

grid.arrange(p_cmh_32+coord_cartesian(ylim = c(0, N)),p_cmh_35+coord_cartesian(ylim = c(0, N)), ncol = 2)
grid.arrange(p_cmh_50+coord_cartesian(ylim = c(0, N)),p_cmh_100+coord_cartesian(ylim = c(0, N)), ncol = 2) 



######################################################################################
############# alllele frequency##########################################################
###################################################################################

cmh_100<-data.frame(value =p.values_cmh_100,position= c(1:length(p.values_cmh_100)),Chr = Chr_100)

D_Allele.freq= af(CD.sync,gen=D_generation,repl=D_replicates)
D_Allele.coverage = coverage(CD.sync,gen=D_generation,repl=D_replicates)



D_Chr <- gsub("\\..*", "", rownames(D_Allele.freq))
D_position<- gsub("^.*?\\.(.*?)$", "\\1", rownames(D_Allele.freq))

rownames(cmh_100)<-rownames(D_Allele.freq)

#extract data set with the chromosome V
D_CHRIV <- subset(cmh_100, gsub("\\..*", "", rownames(cmh_100)) == "IV")

#high p.values_cmh_100
indices <- which(D_CHRIV$value >=75 )

data_filtered1 <-numeric(0)
for (i in indices){
  data_filtered1 <- rbind(data_filtered1,D_CHRIV[ i , ])
}

# Find index of the row with correlated allele 
index_max <- which.max(data_filtered1$value)

# define indiex of rows defore and after
start_index <- max(1, index_max -5)
end_index <- min(nrow(data_filtered1), index_max + (nrow(data_filtered1)-index_max))

# Filter the data.frame of correlated allele
data_filtered <- data_filtered1[start_index:end_index, ]
#Create frequency data
data_D<-as.data.frame(matrix(NA,ncol(D_Allele.freq), nrow(data_filtered )))

for(i in 1:nrow(data_filtered )) {
  data_D[,i] <- cbind(t(subset(D_Allele.freq, rownames(D_Allele.freq) == rownames(data_filtered)[i])))
  
}

#extract allele frequency for replicate 1
frequence_D <- numeric(0)

D_data<-cbind(t(data_D)[,1],t(data_D)[,2:4])
for (i in 1:nrow(D_data)) {
  frequence_D <- c(frequence_D, D_data[i, ])
}

CD1<- data.frame(
  generation = rep(c(0,35,50,100),nrow(data_filtered)),Chr = rep(c(1:nrow(data_filtered )), each=4),
  frequency = frequence_D)

#Median frequency replicate 1
median_freq1 <- CD1 %>%
  group_by(generation) %>%
  summarize(median_frequency = median(frequency))

median_freq1<- cbind(median_freq1,Chr=rep("a",4))

#Ploting allele frequncy change
D1<-ggplot(CD1, aes(x = generation, y = frequency,group = as.character(Chr))) +
  geom_line()+
  geom_line(data = median_freq1, aes(x = generation, y = median_frequency), color = "red", linewidth = 1)  +
  coord_cartesian(ylim = c(0, 0.60)) + 
  ggtitle("1")+theme_classic()+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))


###############samething is perform for other replicate #############

D_data2<-cbind(t(data_D)[,1],t(data_D)[,5:6])

frequence_D2 <- numeric(0) 

for (i in 1:nrow(D_data2)) {
  frequence_D2 <- c(frequence_D2, D_data2[i, ])
}

CD2<- data.frame(
  generation = rep(c(0,35,100),nrow(data_filtered)),Chr = rep(c(1:nrow(data_filtered )), each=3),
  frequency = frequence_D2)

median_freq2 <- CD2 %>%
  group_by(generation) %>%
  summarize(median_frequency = median(frequency))

median_freq2<- cbind(median_freq2,Chr=rep("a",3))


D2<-ggplot(CD2, aes(x = generation, y = frequency,group = as.character(Chr))) +
  geom_line()+
  geom_line(data = median_freq2, aes(x = generation, y = median_frequency), color = "red", linewidth = 1)  +
  coord_cartesian(ylim = c(0, 0.60)) + 
  ggtitle("2")+theme_classic()+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))





D_data3<-cbind(t(data_D)[,1],t(data_D)[,7:9])

frequence_D3 <- numeric(0) 

for (i in 1:nrow(D_data3)) {
  frequence_D3 <- c(frequence_D3, D_data3[i, ])
}

CD3<- data.frame(
  generation = rep(c(0,35,50,100),nrow(data_filtered)),Chr = rep(c(1:nrow(data_filtered )), each=4),
  frequency = frequence_D3)

median_freq3 <- CD3 %>%
  group_by(generation) %>%
  summarize(median_frequency = median(frequency))

median_freq3<- cbind(median_freq3,Chr=rep("a",4))



D3<-ggplot(CD3, aes(x = generation, y = frequency,group = as.character(Chr))) +
  geom_line()+
  geom_line(data = median_freq3, aes(x = generation, y = median_frequency), color = "red", linewidth = 1)  +
  coord_cartesian(ylim = c(0, 0.60)) + 
  ggtitle("3")+theme_classic()+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))


D_data4<-cbind(t(data_D)[,1],t(data_D)[,10:12])

frequence_D4 <- numeric(0) 

for (i in 1:nrow(D_data4)) {
  frequence_D4 <- c(frequence_D4, D_data4[i, ])
}

CD4<- data.frame(
  generation = rep(c(0,35,50,100),nrow(data_filtered)),Chr = rep(c(1:nrow(data_filtered )), each=4),
  frequency = frequence_D4)

median_freq4 <- CD4 %>%
  group_by(generation) %>%
  summarize(median_frequency = median(frequency))

median_freq4<- cbind(median_freq4,Chr=rep("a",4))


D4<-ggplot(CD4, aes(x = generation, y = frequency,group = as.character(Chr))) +
  geom_line() +
  geom_line(data = median_freq4, aes(x = generation, y = median_frequency), color = "red", linewidth = 1) +
  coord_cartesian(ylim = c(0, 0.60)) + 
  ggtitle("4")+theme_classic()+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))


D_data5<-cbind(t(data_D)[,1],t(data_D)[,13:14])

frequence_D5 <- numeric(0) 

for (i in 1:nrow(D_data5)) {
  frequence_D5 <- c(frequence_D5, D_data5[i, ])
}

CD5<- data.frame(
  generation = rep(c(0,32,100),nrow(data_filtered)),Chr = rep(c(1:nrow(data_filtered )), each=3),
  frequency = frequence_D5)

median_freq5 <- CD5 %>%
  group_by(generation) %>%
  summarize(median_frequency = median(frequency))

median_freq5<- cbind(median_freq5,Chr=rep("a",3))


D5<-ggplot(CD5, aes(x = generation, y = frequency,group = as.character(Chr))) +
  geom_line()+
  geom_line(data = median_freq5, aes(x = generation, y = median_frequency), color = "red", linewidth = 1)  +
  coord_cartesian(ylim = c(0, 0.60)) + 
  ggtitle("5")+theme_classic()+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

D_data6<-cbind(t(data_D)[,1],t(data_D)[,15:16])

frequence_D6 <- numeric(0) 

for (i in 1:nrow(D_data6)) {
  frequence_D6 <- c(frequence_D6, D_data6[i, ])
}

CD6<- data.frame(
  generation = rep(c(0,32,100),nrow(data_filtered)),Chr = rep(c(1:nrow(data_filtered )), each=3),
  frequency = frequence_D6)

median_freq6 <- CD6 %>%
  group_by(generation) %>%
  summarize(median_frequency = median(frequency))

median_freq6<- cbind(median_freq6,Chr=rep("a",3))


D6<-ggplot(CD6, aes(x = generation, y = frequency,group = as.character(Chr))) +
  geom_line()+
  geom_line(data = median_freq6, aes(x = generation, y = median_frequency), color = "red", linewidth = 1)  +
  coord_cartesian(ylim = c(0, 0.60)) + 
  ggtitle("6")+theme_classic()+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))


D_data7<-cbind(t(data_D)[,1],t(data_D)[,17:18])

frequence_D7 <- numeric(0) 

for (i in 1:nrow(D_data7)) {
  frequence_D7 <- c(frequence_D7, D_data6[i, ])
}

CD7<- data.frame(
  generation = rep(c(0,32,100),nrow(data_filtered)),Chr = rep(c(1:nrow(data_filtered )), each=3),
  frequency = frequence_D7)

median_freq7 <- CD7 %>%
  group_by(generation) %>%
  summarize(median_frequency = median(frequency))

median_freq7<- cbind(median_freq7,Chr=rep("a",3))


D7<-ggplot(CD7, aes(x = generation, y = frequency,group = as.character(Chr))) +
  geom_line() +
  geom_line(data = median_freq7, aes(x = generation, y = median_frequency), color = "red", linewidth = 1)  +
  coord_cartesian(ylim = c(0, 0.60)) + 
  ggtitle("7")+theme_classic()+ 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))


grid.arrange(D1 , D2, D3, D4 , D5 , D6, D7, ncol = 7)

####################################################################################
################################################################################

