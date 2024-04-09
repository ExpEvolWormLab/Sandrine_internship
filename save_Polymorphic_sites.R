#!/usr/bin/env Rscript

#import library
library(data.table)
library(poolSeq)


#import the corresponding sync file containing all bam file of different population and population
args = commandArgs(trailingOnly=TRUE)

#call population name
Sample_names=t(read.table(args[1]))
number_samples=length(Sample_names)

#Androdioecy Generation
Androdioecy_generation <- c(0,0,100,10,36,50,68,100,10,36,50,68,100,10,36,50,68,100,10,32,66,100,10,32,66,10,32)

# Androdioecy Replicates 
Androdioecy_replicates<- c(0,0,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,5,5,5,5,6,6)

#Dioecy Generation
Dioecy_generation <-c(35,50,100,35,100,35,50,100,35,50,100,32,100,32,100,32,100,0)

# Dioecy Replicates
Dioecy_replicates<- c(1,1,1,2,2,3,3,3,4,4,4,5,5,6,6,7,7,0)

#Androdioecy $ Dioecy Generation

All_Generation <- c(Androdioecy_generation,Dioecy_generation)

#Androdioecy $ Dioecy Generation

All_Replicates <-c(Androdioecy_replicates,Dioecy_replicates)


#Read Poolseq.sync file
Poolseq.sync=read.sync(args[2],gen=All_Generation,repl=All_Replicates )

#extract alleles frequency
all.alleles_frequency = af(Poolseq.sync,gen=All_Generation,repl=All_Replicates )

#extract alleles coverage
all.alleles_coverage = coverage(Poolseq.sync,gen=All_Generation,repl=All_Replicates )
rm(Poolseq.sync);gc()

## Only keep lines for which we have at least 10 non-NA values
na_values <- rowSums(is.na(all.alleles_frequency))
keep = (na_values<=(number_samples-10))

## Subset for mean coverage < 125 & >15
mean_cov <- rowMeans(all.alleles_coverage)
keep = keep & mean_cov<=125 & mean_cov>=15

# Extract polymorphic sites - at least 10 counts
alternative_counts <- all.alleles_frequency*all.alleles_coverage
summed_alternative_alleles <- rowSums(round(alternative_counts),na.rm=TRUE)
keep = keep & (summed_alternative_alleles> 10)

#create table with position and chromosome number 
temp=row.names(all.alleles_frequency)[keep]
temp.CHR=tstrsplit(temp,"[.]")[[1]]
temp.POS=tstrsplit(temp,"[.]")[[2]]

kept_pos <- data.frame(CHR=temp.CHR,POS=temp.POS)
write.table(kept_pos,file=sprintf("%s.mean",args[3]),quote=FALSE,row.names = FALSE,col.names=FALSE)


