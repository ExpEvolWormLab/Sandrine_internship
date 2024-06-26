library(poolSeq)
library(ggplot2)
library(gridExtra)


#Generation
A_generation <-c(0,100,36,50,68,100,10,36,50,68,100,10,36,50,68,100,10,32,66,100,10,32,66,10,32)
#Replicates
A_replicates<- c(0,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,5,5,5,5,6,6)

#Read sync file

CA.sync<-read.sync("CA.sync",gen=A_generation,repl=A_replicates)

#extract alleles frequency and coverage for generation 0
Allele.freq0= af(CA.sync,gen=0,repl=A_replicates)
Allele.coverage0 = coverage(CA.sync,gen=0,repl=A_replicates)

# alternative allele count for generation 0
A0 <- Allele.freq0*Allele.coverage0

# reference allele count for generation 0
a0 <- Allele.coverage0 - A0

#extract alleles frequency and coverage for generation 10
Allele.freq10= af(CA.sync,gen=10,repl=A_replicates)

Allele.coverage10 = coverage(CA.sync,gen=10,repl=A_replicates)

# alternative allele count for generation 10
A10<-Allele.freq10*Allele.coverage10

# reference allele count for generation 10
a10<-Allele.coverage10 - A10

Chr_10 <- gsub("\\..*", "", rownames(Allele.freq10))

# perform CMH-test for all loci
p.values_cmh_10 <- cmh.test(A0=rbind(A0,A0,A0,A0,A0), a0=rbind(a0,a0,a0,a0,a0), At=t(A10), at=t(a10),min.cov=1, max.cov=1, min.cnt=1, log=TRUE)

cmh_10<-data.frame(value =p.values_cmh_10,position= c(1:length(p.values_cmh_10)),Chr = Chr_10 )

p_cmh_10<-ggplot(na.omit(cmh_10), aes(x=position, y=value,group=Chr))+
  geom_point() + geom_point(aes(colour = Chr))+labs(x = "Position ", y = "-log10(p) ",title = "CAG10")


#extract alleles frequency and coverage for generation 32
Allele.freq32= af(CA.sync,gen=32,repl=A_replicates)

Allele.coverage32 = coverage(CA.sync,gen=32,repl=A_replicates)

# alternative and reference allele count for generation 32
A32<-Allele.freq32*Allele.coverage32
a32<-Allele.coverage32 - A32

#extract number for generation 32
Chr_32 <- gsub("\\..*", "", rownames(Allele.freq32))
# perform CMH-test for all empirical loci
p.values_cmh_32 <- cmh.test(A0=rbind(A0,A0,A0), a0=rbind(a0,a0,a0), At=t(A32), at=t(a32),min.cov=1, max.cov=1, min.cnt=1, log=TRUE)

cmh_32<-data.frame(value =p.values_cmh_32,position= c(1:length(p.values_cmh_32)),Chr = Chr_32 )

p_cmh_32<-ggplot(na.omit(cmh_32), aes(x=position, y=value,group=Chr))+
  geom_point() + geom_point(aes(colour = Chr))+labs(x = "Position ", y = "-log10(p) ",title = "CAG32")



#extract alleles frequency and coverage for generation 36
Allele.freq36= af(CA.sync,gen=36,repl=A_replicates)

Allele.coverage36 = coverage(CA.sync,gen=36,repl=A_replicates)

#alternative end reference allele count
A36<-Allele.freq36*Allele.coverage36
a36<-Allele.coverage36 - A36

#extract chromosome number for all loci
Chr_36 <- gsub("\\..*", "", rownames(Allele.freq36))

# perform CMH-test for all empirical loci
p.values_cmh_36 <- cmh.test(A0=rbind(A0,A0,A0), a0=rbind(a0,a0,a0), At=t(A36), at=t(a36),min.cov=1, max.cov=1, min.cnt=1, log=TRUE)

cmh_36<-data.frame(value =p.values_cmh_36,position= c(1:length(p.values_cmh_36)),Chr = Chr_36 )

p_cmh_36<-ggplot(na.omit(cmh_36), aes(x=position/1e4, y=value,group=Chr))+
  geom_point() + geom_point(aes(colour = Chr))+labs(x = "Position ", y = "-log10(p) ",title = "CAG36")


#extract alleles frequency and coverage for generation 50
Allele.freq50= af(CA.sync,gen=50,repl=A_replicates)

Allele.coverage50 = coverage(CA.sync,gen=50,repl=A_replicates)

#alternative end reference allele count
A50<-Allele.freq50*Allele.coverage50
a50<-Allele.coverage50 - A50

Chr_50 <- gsub("\\..*", "", rownames(Allele.freq50))

# perform CMH-test for all empirical loci
p.values_cmh_50 <- cmh.test(A0=rbind(A0,A0,A0), a0=rbind(a0,a0,a0), At=t(A50), at=t(a50),min.cov=1, max.cov=1, min.cnt=1, log=TRUE)

cmh_50<-data.frame(value =p.values_cmh_50,position= c(1:length(p.values_cmh_50)),Chr = Chr_50 )

p_cmh_50<-ggplot(na.omit(cmh_50), aes(x=position, y=value,group=Chr))+
  geom_point() + geom_point(aes(colour = Chr))+labs(x = "Position ", y = "-log10(p) ",title = "CAG50")



#extract alleles frequency and coverage for generation 66
Allele.freq66= af(CA.sync,gen=66,repl=A_replicates)

Allele.coverage66 = coverage(CA.sync,gen=66,repl=A_replicates)

#alternative end reference allele count
A66<-Allele.freq66*Allele.coverage66
a66<-Allele.coverage66 - A66

Chr_66 <- gsub("\\..*", "", rownames(Allele.freq66))

# perform CMH-test for all empirical loci
p.values_cmh_66 <- cmh.test(A0=rbind(A0,A0), a0=rbind(a0,a0), At=t(A66), at=t(a66),min.cov=1, max.cov=1, min.cnt=1, log=TRUE)

cmh_66<-data.frame(value =p.values_cmh_66,position= c(1:length(p.values_cmh_66)),Chr = Chr_66 )

p_cmh_66<-ggplot(cmh_66, aes(x=position, y=value,group=Chr))+
  geom_point() + geom_point(aes(colour = Chr))+labs(x = "Position ", y = "-log10(p) ",title = "CAG66")


#extract alleles frequency and coverage for generation 68
Allele.freq68= af(CA.sync,gen=68,repl=A_replicates)

Allele.coverage68 = coverage(CA.sync,gen=68,repl=A_replicates)

#alternative end reference allele count
A68<-Allele.freq68*Allele.coverage68
a68<-Allele.coverage68 - A68

Chr_68 <- gsub("\\..*", "", rownames(Allele.freq68))

# perform CMH-test for all empirical loci
p.values_cmh_68 <- cmh.test(A0=rbind(A0,A0,A0), a0=rbind(a0,a0,a0), At=t(A68), at=t(a68),min.cov=1, max.cov=1, min.cnt=1, log=TRUE)

cmh_68<-data.frame(value =p.values_cmh_68,position= c(1:length(p.values_cmh_68)),Chr = Chr_68 )

p_cmh_68<-ggplot(cmh_68, aes(x=position, y=value,group=Chr))+
  geom_point() + geom_point(aes(colour = Chr))+labs(x = "Position ", y = "-log10(p) ",title = "CAG68")



#extract alleles frequency and coverage for generation 100
Allele.freq100= af(CA.sync,gen=100,repl=A_replicates)

Allele.coverage100 = coverage(CA.sync,gen=100,repl=A_replicates)


A100<-Allele.freq100*Allele.coverage100
a100<-Allele.coverage100 - A100

Chr_100 <- gsub("\\..*", "", rownames(Allele.freq100))

# perform CMH-test for all empirical loci
p.values_cmh_100 <- cmh.test(A0=rbind(A0,A0,A0,A0,A0), a0=rbind(a0,a0,a0,a0,a0), At=t(A100), at=t(a100),min.cov=1, max.cov=1, min.cnt=1, log=TRUE)

cmh_100<-data.frame(value =p.values_cmh_100,position= c(1:length(p.values_cmh_100)),Chr = Chr_100)

p_cmh_100<-ggplot(na.omit(cmh_100), aes(x=position/1e4, y=value,group=Chr))+
  geom_point() + geom_point(aes(colour = Chr))+labs(x = "Position ", y = "-log10(p) ",title = "CAG100")



N= max(p.values_cmh_10, na.rm = TRUE)

grid.arrange(p_cmh_10+coord_cartesian(ylim = c(0, N)),p_cmh_50+coord_cartesian(ylim = c(0, N)), ncol = 2)
grid.arrange(p_cmh_66+coord_cartesian(ylim = c(0, N)),p_cmh_100+coord_cartesian(ylim = c(0, N)), ncol = 2) 


