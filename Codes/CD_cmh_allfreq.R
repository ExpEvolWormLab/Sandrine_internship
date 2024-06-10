library(poolSeq)
library(ggplot2)
library(gridExtra)

#Generation
D_generation <-c(35,50,100,35,100,35,50,100,35,50,100,32,100,32,100,32,100,0)
#Replicates
D_replicates<- c(1:18)

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
  geom_point() + geom_point(aes(colour = Chr))+labs(x = "Position ", y = "-log10(p) ",title = "DG32")

# perform SHI-test for all empirical loci
for(i in 1:ncol(D32)){
  p.values_chisq_32 <- cbind(chi.sq.test(A0=D0, a0=d0, At=D32[,i] ,at=d32[,i],min.cov=20, max.cov=0.99, min.cnt=1, log=TRUE))
}

chisq_32<-data.frame(value =p.values_chisq_32,position= c(1:length(p.values_chisq_32)),Chr = Chr_32 )

p_chisq_32<-ggplot(na.omit(chisq_32), aes(x=position, y=value,group=Chr))+
  geom_point() + geom_point(aes(colour = Chr))+labs(x = "Position ", y = "-log10(p) ",title = "DG32")


#extract alleles frequency and coverage for generation 35
Allele.freq35= af(CD.sync,gen=35,repl=D_replicates)

Allele.coverage35 = coverage(CD.sync,gen=35,repl=D_replicates)

D35<-Allele.freq35*Allele.coverage35
d35<-Allele.coverage35 - D35

Chr_35 <- gsub("\\..*", "", rownames(Allele.freq35))

# perform CMH-test for all empirical loci
p.values_cmh_35 <- cmh.test(A0=rbind(D0,D0,D0,D0), a0=rbind(d0,d0,d0,d0), At=t(D35), at=t(d35),min.cov=1, max.cov=1, min.cnt=1, log=TRUE)

cmh_35<-data.frame(value =p.values_cmh_35,position= c(1:length(p.values_cmh_35)),Chr = Chr_35)

p_cmh_35<-ggplot(na.omit(cmh_35), aes(x=position, y=value,group=Chr))+
  geom_point() + geom_point(aes(colour = Chr))+labs(x = "Position ", y = "-log10(p) ",title = "DG35")

# perform CHI-test for all empirical loci
for(i in 1:ncol(D35)){
  p.values_chisq_35 <- cbind(chi.sq.test(A0=D0, a0=d0, At=D35[,i] ,at=d35[,i],min.cov=20, max.cov=0.99, min.cnt=1, log=TRUE))
}

chisq_35<-data.frame(value =p.values_chisq_35,position= c(1:length(p.values_chisq_35)),Chr = Chr_35 )

p_chisq_35<-ggplot(na.omit(chisq_35), aes(x=position, y=value,group=Chr))+
  geom_point() + geom_point(aes(colour = Chr))+labs(x = "Position ", y = "-log10(p) ",title = "DG35")

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
  geom_point() + geom_point(aes(colour = Chr))+labs(x = "Position ", y = "-log10(p) ",title = "DG50")

# perform CHI-test for all empirical loci
for(i in 1:ncol(D50)){
  p.values_chisq_50 <- cbind(chi.sq.test(A0=D0, a0=d0, At=D50[,i] ,at=d50[,i],min.cov=20, max.cov=0.99, min.cnt=1, log=TRUE))
}

chisq_50<-data.frame(value =p.values_chisq_50,position= c(1:length(p.values_chisq_50)),Chr = Chr_50 )

p_chisq_50<-ggplot(na.omit(chisq_50), aes(x=position, y=value,group=Chr))+
  geom_point() + geom_point(aes(colour = Chr))+labs(x = "Position ", y = "-log10(p) ",title = "DG50")



#extract alleles frequency and coverage for generation 100
Allele.freq100= af(CD.sync,gen=100,repl=D_replicates)

Allele.coverage100 = coverage(CD.sync,gen=100,repl=D_replicates)

D100<-Allele.freq100*Allele.coverage100
d100<-Allele.coverage100 - D100

Chr_100 <- gsub("\\..*", "", rownames(Allele.freq100))

# perform CMH-test for all empirical loci
p.values_cmh_100 <- cmh.test(A0=rbind(D0,D0,D0,D0,D0,D0,D0), a0=rbind(d0,d0,d0,d0,d0,d0,d0), At=t(D100), at=t(d100),min.cov=1, max.cov=1, min.cnt=1, log=TRUE)

cmh_100<-data.frame(value =p.values_cmh_100,position= c(1:length(p.values_cmh_100)),Chr = Chr_100)

p_cmh_100<-ggplot(na.omit(cmh_100), aes(x=position, y=value,group=Chr))+
  geom_point() + geom_point(aes(colour = Chr))+labs(x = "Position ", y = "-log10(p) ",title = "DG100")

# perform CHI-test for all empirical loci
for(i in 1:ncol(D100)){
  p.values_chisq_100 <- cbind(chi.sq.test(A0=D0, a0=d0, At=D100[,7] ,at=d100[,7],min.cov=20, max.cov=0.99, min.cnt=1, log=TRUE))
}

chisq_100<-data.frame(value =p.values_chisq_100,position= c(1:length(p.values_chisq_100)),Chr = Chr_100)

p_chisq_100<-ggplot(na.omit(chisq_100), aes(x=position, y=value,group=Chr))+
  geom_point() + geom_point(aes(colour = Chr))+labs(x = "Position ", y = "-log10(p) ",title = "DG100")

N= max(p.values_cmh_100, na.rm = TRUE)

grid.arrange(p_cmh_32+coord_cartesian(ylim = c(0, N)),p_cmh_35+coord_cartesian(ylim = c(0, N)), ncol = 2)
grid.arrange(p_cmh_50+coord_cartesian(ylim = c(0, N)),p_cmh_100+coord_cartesian(ylim = c(0, N)), ncol = 2) 
#grid.arrange(p_cmh_10,p_cmh_36, p_cmh_50,p_cmh_66,p_cmh_68, p_cmh_100, ncol = 2) 

M= max(p.values_chisq_100, na.rm = TRUE)

grid.arrange(p_chisq_32+coord_cartesian(ylim = c(0, M)),p_chisq_35+coord_cartesian(ylim = c(0, M)),p_chisq_50+coord_cartesian(ylim = c(0, M)),p_chisq_100+coord_cartesian(ylim = c(0, M)), ncol = 2)
grid.arrange(p_chisq_50+coord_cartesian(ylim = c(0, M)),p_chisq_100+coord_cartesian(ylim = c(0, M)), ncol = 2) 
