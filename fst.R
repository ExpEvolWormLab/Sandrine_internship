library(poolfstat)
library(ggplot2)


####################################################################################################################
#######  Fst over generation for Androdioecy population ################################################
##################################################################################################################


#Convert Popoolation Sync files into a pooldata object
CA1G100 <- popsync2pooldata("CA1G100.sync", poolsizes = c(1000,1000))

#compute Fst considering window
fstA1G100<-computeFST(CA1G100, sliding.window.size=50)


CA1G10 <- popsync2pooldata("CA1G10.sync", poolsizes = c(1000,1000))
fstA1G10<-computeFST(CA1G10, sliding.window.size=50)

CA1G36 <- popsync2pooldata("CA1G36.sync", poolsizes = c(1000,1000))
fstA1G36<-computeFST(CA1G36,sliding.window.size=50)



CA1G50 <- popsync2pooldata("CA1G50.sync", poolsizes = c(1000,1000))
fstA1G50<-computeFST(CA1G50, sliding.window.size=50)

CA1G68 <- popsync2pooldata("CA1G68.sync", poolsizes = c(1000,1000))
fstA1G68<-computeFST(CA1G68, sliding.window.size=50)

CA2G100 <- popsync2pooldata("CA2G100.sync", poolsizes = c(1000,1000))
fstA2G100<-computeFST(CA2G100, sliding.window.size=50)


CA2G10 <- popsync2pooldata("CA2G10.sync", poolsizes = c(1000,1000))
fstA2G10<-computeFST(CA2G10, sliding.window.size=50)


CA2G36 <- popsync2pooldata("CA2G36.sync", poolsizes = c(1000,1000))
fstA2G36<-computeFST(CA2G36, sliding.window.size=50)


CA2G50 <- popsync2pooldata("CA2G50.sync", poolsizes = c(1000,1000))
fstA2G50<-computeFST(CA2G50, sliding.window.size=50)

CA2G68 <- popsync2pooldata("CA2G68.sync", poolsizes = c(1000,1000))
fstA2G68<-computeFST(CA2G68, sliding.window.size=50)

CA3G100 <- popsync2pooldata("CA3G100.sync", poolsizes = c(1000,1000))
fstA3G100<-computeFST(CA3G100, sliding.window.size=50)

CA3G10 <- popsync2pooldata("CA3G10.sync", poolsizes = c(1000,1000))
fstA3G10<-computeFST(CA3G10, sliding.window.size=50)

CA3G36 <- popsync2pooldata("CA3G36.sync", poolsizes = c(1000,1000))
fstA3G36<-computeFST(CA3G36, sliding.window.size=50)


CA3G50 <- popsync2pooldata("CA3G50.sync", poolsizes = c(1000,1000))
fstA3G50<-computeFST(CA3G50, sliding.window.size=50)


CA3G68 <- popsync2pooldata("CA3G68.sync", poolsizes = c(1000,1000))
fstA3G68<-computeFST(CA3G68, sliding.window.size=50)


CA4G100 <- popsync2pooldata("CA4G100.sync", poolsizes = c(1000,1000))
fstA4G100<-computeFST(CA4G100, sliding.window.size=50)

CA4G10 <- popsync2pooldata("CA4G10.sync", poolsizes = c(1000,1000))
fstA4G10<-computeFST(CA4G10, sliding.window.size=50)

CA4G32 <- popsync2pooldata("CA4G32.sync", poolsizes = c(1000,1000))
fstA4G32<-computeFST(CA4G32, sliding.window.size=50)


CA4G66 <- popsync2pooldata("CA4G66.sync", poolsizes = c(1000,1000))
fstA4G66<-computeFST(CA4G66, sliding.window.size=50)

CA5G100 <- popsync2pooldata("CA5G100.sync", poolsizes = c(1000,1000))
fstA5G100<-computeFST(CA5G100, sliding.window.size=50)

CA5G10 <- popsync2pooldata("CA5G10.sync", poolsizes = c(1000,1000))
fstA5G10<-computeFST(CA5G10, sliding.window.size=50)

CA5G32 <- popsync2pooldata("CA5G32.sync", poolsizes = c(1000,1000))
fstA5G32<-computeFST(CA5G32, sliding.window.size=50)
CA5G66 <- popsync2pooldata("CA5G66.sync", poolsizes = c(1000,1000))
fstA5G66<-computeFST(CA5G66, sliding.window.size=50)


CA6G10 <- popsync2pooldata("CA6G10.sync", poolsizes = c(1000,1000))
fstA6G10<-computeFST(CA6G10, sliding.window.size=50)

CA6G32 <- popsync2pooldata("CA6G32.sync", poolsizes = c(1000,1000))
fstA6G32<-computeFST(CA6G32, sliding.window.size=50)

FSTA1=c(fstA1G36$FST,fstA1G50$FST,fstA1G68$FST,fstA1G100$FST)
FSTA2=c(fstA2G10$FST,fstA2G36$FST,fstA2G50$FST,fstA2G68$FST,fstA2G100$FST)
FSTA3=c(fstA3G10$FST,fstA3G36$FST,fstA3G50$FST,fstA3G68$FST,fstA3G100$FST)
FSTA4=c(fstA4G10$FST,fstA4G32$FST,fstA4G66$FST,fstA4G100$FST)
FSTA5=c(fstA5G10$FST,fstA5G32$FST,fstA5G66$FST,fstA5G100$FST)
FSTA6=c(fstA6G10$FST,fstA6G32$FST)

######   Plot every replicate fst in the same plot ############################

A_fst= c(FSTA1,FSTA2,FSTA3,FSTA4,FSTA5,FSTA6)
A_generations=c(" 36"," 50"," 68","100"," 10"," 36"," 50"," 68","100"," 10"," 36"," 50"," 68","100"," 10"," 32"," 66","100"," 10"," 32"," 66","100"," 10"," 32")

population_A=data.frame(replicate=c(rep(c("CA1","CA2","CA3","CA4","CA5","CA6"),c(4,5,5,4,4,2))),generation=A_generations,Fst=A_fst)

p_A <- ggplot(population_A, aes(x=generation, y=Fst,group = replicate)) + 
  geom_line(aes(colour = replicate)) + geom_point(aes(colour = replicate))+coord_cartesian(ylim = c(0, 1)) 


###### Save image plot in  PDF file ############

pdf(file = "A_Fst.pdf")

ggplot(population_A, aes(x=generation, y=Fst,group = replicate)) + 
  geom_line(aes(colour = replicate)) + geom_point(aes(colour = replicate))+coord_cartesian(ylim = c(0, 0.03))

dev.off()


####################################################################################################################
#######  Fst over generation for Dioecy population ################################################
##################################################################################################################

CD1G35 <- popsync2pooldata("CD1G35.sync", poolsizes = c(1000,1000))
fstD1G35<-computeFST(CD1G35, sliding.window.size=50)

CD1G50 <- popsync2pooldata("CD1G50.sync", poolsizes = c(1000,1000))
fstD1G50<-computeFST(CD1G50, sliding.window.size=50)

CD1G100 <- popsync2pooldata("CD1G100.sync", poolsizes = c(1000,1000))
fstD1G100<-computeFST(CD1G100, sliding.window.size=50)

CD2G35 <- popsync2pooldata("CD2G35.sync", poolsizes = c(1000,1000))
fstD2G35<-computeFST(CD2G35, sliding.window.size=50)

CD2G100 <- popsync2pooldata("CD2G100.sync", poolsizes = c(1000,1000))
fstD2G100<-computeFST(CD2G100, sliding.window.size=50)

CD3G35 <- popsync2pooldata("CD3G35.sync", poolsizes = c(1000,1000))
fstD3G35<-computeFST(CD3G35, sliding.window.size=50)

CD3G50 <- popsync2pooldata("CD3G50.sync", poolsizes = c(1000,1000))
fstD3G50<-computeFST(CD3G50, sliding.window.size=50)


CD3G100 <- popsync2pooldata("CD3G100.sync", poolsizes = c(1000,1000))
fstD3G100<-computeFST(CD3G100, sliding.window.size=50)

CD4G35 <- popsync2pooldata("CD4G35.sync", poolsizes = c(1000,1000))
fstD4G35<-computeFST(CD4G35, sliding.window.size=50)

CD4G50 <- popsync2pooldata("CD4G50.sync", poolsizes = c(1000,1000))
fstD4G50<-computeFST(CD4G50, sliding.window.size=50)

CD4G100 <- popsync2pooldata("CD4G100.sync", poolsizes = c(1000,1000))
fstD4G100<-computeFST(CD4G100, sliding.window.size=50)

CD5G32 <- popsync2pooldata("CD5G32.sync", poolsizes = c(1000,1000))
fstD5G32<-computeFST(CD5G32, sliding.window.size=50)

CD5G100 <- popsync2pooldata("CD5G100.sync", poolsizes = c(1000,1000))
fstD5G100<-computeFST(CD5G100, sliding.window.size=50)

CD6G32 <- popsync2pooldata("CD6G32.sync", poolsizes = c(1000,1000))
fstD6G32<-computeFST(CD6G32, sliding.window.size=50)

CD6G100 <- popsync2pooldata("CD6G100.sync", poolsizes = c(1000,1000))
fstD6G100<-computeFST(CD6G100, sliding.window.size=50)

CD7G32 <- popsync2pooldata("CD6G32.sync", poolsizes = c(1000,1000))
fstD7G32<-computeFST(CD7G32, sliding.window.size=50)

CD7G100 <- popsync2pooldata("CD7G100.sync", poolsizes = c(1000,1000))
fstD7G100<-computeFST(CD7G100, sliding.window.size=50)

FSTD1=c(fstD1G35$FST,fstD1G50$FST,fstD1G100$FST)
FSTD2=c(fstD2G35$FST,fstD2G100$FST)
FSTD3=c(fstD3G35$FST,fstD3G50$FST,fstD3G100$FST)
FSTD4=c(fstD4G35$FST,fstD4G50$FST,fstD4G100$FST)
FSTD5=c(fstD5G32$FST,fstD5G100$FST)
FSTD6=c(fstD6G32$FST,fstD6G100$FST)
FSTD7=c(fstD7G32$FST,fstD7G100$FST)

######   Plot every replicate fst in the same plot ############################
  

D_fst = c(FSTD1,FSTD2,FSTD3,FSTD4,FSTD5,FSTD6,FSTD7)
D_generations = c(" 35"," 50","100"," 35","100"," 35"," 50","100"," 35"," 50","100"," 32","100"," 32","100"," 32","100")

population_D=data.frame(replicate=c(rep(c("CD1","CD2","CD3","CD4","CD5","CD6","CD7"),c(3,2,3,3,2,2,2))),generation=D_generations,Fst=D_fst)

###### Save image plot in  PDF file ############
pdf(file = "D_Fst.pdf")

ggplot(population_D, aes(x=generation, y=Fst,group = replicate)) + 
  geom_line(aes(colour = replicate)) + geom_point(aes(colour = replicate))+coord_cartesian(ylim = c(0, 0.06))

dev.off()


####################################################################################################################
#######  Average Fst over generation for CD and CA ################################################
##################################################################################################################

# Function to calculate mean and confidence interval
calc_mean_ci <- function(data_vector) {
  mean_value <- mean(data_vector)
  ci <- t.test(data_vector)$conf.int
  return(c(mean_value, ci[1], ci[2]))
}


Dgeneration=c(" 32"," 35"," 50","100")

#combine fst value per generation

D<- list(DG32 = c(fstD5G32$FST, fstD6G32$FST, fstD7G32$FST),
         DG35 = c(fstD1G35$FST, fstD2G35$FST, fstD3G35$FST, fstD4G35$FST),
         DG50 = c(fstD1G50$FST, fstD3G50$FST, fstD4G50$FST),
         DG100 = c(fstD1G100$FST, fstD2G100$FST, fstD3G100$FST, 
                   fstD4G100$FST, fstD5G100$FST, fstD6G100$FST, fstD7G100$FST)
)

CD <- setNames(D,Dgeneration)

# Calculate mean and confidence interval for each generation
mean_ci_D  <- lapply(CD, calc_mean_ci)

# Create a data frame for plotting
plot_data_D <- data.frame(
  vector = rep(names(CD), each = 3),
  stat = rep(c("Mean", "CI_Lower", "CI_Upper"), length(CD)),
  value = unlist(stats)
)


Ageneration=c(" 10"," 32"," 36"," 50"," 66"," 68","100")

#combine fst value per generation

A<- list(
        AG10= c(fstA2G10$FST,fstA3G10$FST,fstA4G10$FST,fstA5G10$FST,fstA6G10$FST),
        AG32= c(fstA4G32$FST,fstA5G32$FST,fstA6G32$FST),
        AG36= c(fstA1G36$FST,fstA2G36$FST,fstA3G36$FST),
        AG50= c(fstA1G50$FST,fstA2G50$FST,fstA3G50$FST),
        AG66= c(fstA4G66$FST,fstA5G66$FST),
        AG68= c(fstA1G68$FST,fstA2G68$FST,fstA3G68$FST),
        AG100= c(fstA1G100$FST,fstA2G100$FST,fstA3G100$FST,fstA4G100$FST,fstA5G100$FST)
)


CA <- setNames(A,Ageneration)

# Calculate mean and confidence interval for each generation
mean_ci_A <- lapply(CA, calc_mean_ci)

# Create a data frame for plotting
plot_data_A <- data.frame(
  vector = rep(names(CA), each = 3),
  stat = rep(c("Mean", "CI_Lower", "CI_Upper"), length(CA)),
  value = unlist(mean_ci_A)
)


#Add a column to identify each group
plot_data_A$group <- "CA"
plot_data_D$group <- "CD"

# Combine both data frames
combined_plot_data <- rbind(plot_data_A, plot_data_D)

#plotting
p <- ggplot(combined_plot_data, aes(x = vector, y = value, color = group , linetype =group)) +
  geom_point(data = subset(combined_plot_data, stat == "Mean"), size = 3) +
  geom_errorbar(aes(ymin = ifelse(stat == "CI_Lower", value, NA), 
                    ymax = ifelse(stat == "CI_Upper", value, NA)),
                width = 0.2, color = "black") +
  geom_line(data = subset(combined_plot_data, stat == "Mean")) +
  labs(x = "Generation", y = "Fst")
 
##### Save image plot in  PDF file ############
pdf(file = "DvsA_Fst.pdf")
ggplot(combined_plot_data, aes(x = vector, y = value, color = group , linetype =group)) +
  geom_point(data = subset(combined_plot_data, stat == "Mean"), size = 3) +
  geom_errorbar(aes(ymin = ifelse(stat == "CI_Lower", value, NA), 
                    ymax = ifelse(stat == "CI_Upper", value, NA)),
                width = 0.2, color = "black") +
  geom_line(data = subset(combined_plot_data, stat == "Mean")) +
  labs(x = "Generation", y = "Fst")
dev.off()



#############################################################################
############ plot fst along the genome for Androdioecy population######################################
#############################################################################

p1_100<-ggplot(fstA1G100$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.40))+labs(x = " ",
                                                                        y = " ",title = "CA1G100")

p2_100<-ggplot(fstA2G100$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.40))+labs(x = " ",
                                                                        y = " ",title = "CA2G100")

p3_100<-ggplot(fstA3G100$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.40))+labs(x = " ",
                                                                        y = " ",title = "CA3G100")

p4_100<-ggplot(fstA4G100$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.40))+labs(x = "",
                                                                        y = " ",title = "CA4G100")

p5_100<-ggplot(fstA5G100$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.40))+labs(x = "Position ",
                                                                        y = " ",title = "CA5G100")

p1_10<-ggplot(fstA2G10$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.40))+labs(x = " ",
                                                                        y = "FST",title = "CA2G10")

p2_10<-ggplot(fstA3G10$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.40))+labs(x = " ",
                                                                        y = "FST",title = "CA3G10")

p3_10<-ggplot(fstA4G10$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.40))+labs(x = " ",
                                                                        y = "FST",title = "CA4G10")

p4_10<-ggplot(fstA5G10$sliding.windows.fst, aes(x=CumulatedPosition, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.40))+labs(x = " Position",
                                                                        y = "FST",title = "CA5G10")

p5_10<-ggplot(fstA6G10$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.40))+labs(x = " Position",
                                                                        y = "FST",title = "CA6G10")

p1_68<-ggplot(fstA1G68$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.40))+labs(x = " ",
                                                                        y = "",title = "CA1G68")

p2_68<-ggplot(fstA2G68$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.40))+labs(x = " ",
                                                                        y = "",title = "CA2G68")

p3_68<-ggplot(fstA3G68$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.35))+labs(x = " ",
                                                                        y = "",title = "CA3G68")

p4_66<-ggplot(fstA4G66$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.40))+labs(x = " ",
                                                                        y = " ",title = "CA4G66")

p5_66<-ggplot(fstA5G66$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.40))+labs(x = " Position",
                                                                        y = "",title = "CA5G66")

p1_50<-ggplot(fstA1G50$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.40))+labs(x = " ",
                                                                        y = "",title = "CA1G50")

p2_50<-ggplot(fstA2G50$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.40))+labs(x = " ",
                                                                        y = "",title = "CA2G50")

p3_50<-ggplot(fstA3G50$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.35))+labs(x = " ",
                                                                        y = "",title = "CA3G50")

###### Save image plot in  PDF file ############
pdf(file = "CA_Fst_position.pdf")

grid.arrange(p1_10,p2_50, p2_100,p2_10,p3_50, p3_100,p3_10,p4_66,p4_100, p4_10,p5_66,p5_100, ncol = 3) # Affiche p1 et p2 côte à côte

dev.off()
############




#############################################################################
############ plot fst along the genome for Dioecy population######################################
#############################################################################


p1_100<-ggplot(fstD1G100$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.70))+labs(x = " ",
                                                                        y = " ",title = "CD1G100")

p2_100<-ggplot(fstD2G100$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.70))+labs(x = " ",
                                                                        y = " ",title = "CD2G100")

p3_100<-ggplot(fstD3G100$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.70))+labs(x = " ",
                                                                        y = " ",title = "CD3G100")

p4_100<-ggplot(fstD4G100$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.70))+labs(x = "",
                                                                        y = " ",title = "CD4G100")

p5_100<-ggplot(fstD5G100$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.70))+labs(x = "Position ",
                                                                        y = " ",title = "CD5G100")

p6_100<-ggplot(fstD6G100$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.70))+labs(x = " ",
                                                                        y = "",title = "CD6G100")

p7_100<-ggplot(fstD7G100$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.70))+labs(x = "CumulatedPosition ",
                                                                        y = "",title = "CD7G100")



p1_35<-ggplot(fstD1G35$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.70))+labs(x = " ",
                                                                        y = " FST",title = "CD1G35")

p2_35<-ggplot(fstD2G35$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.70))+labs(x = " ",
                                                                        y = "FST ",title = "CD2G35")

p3_35<-ggplot(fstD3G35$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.70))+labs(x = " ",
                                                                        y = " FST",title = "CD3G35")

p4_35<-ggplot(fstD4G35$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.70))+labs(x = "",
                                                                        y = "FST ",title = "CD4G35")

p5_32<-ggplot(fstD5G32$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.70))+labs(x = "Position ",
                                                                        y = "FST ",title = "CD5G32")

p6_32<-ggplot(fstD6G32$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.70))+labs(x = " ",
                                                                        y = "FST",title = "CD6G32")

p7_32<-ggplot(fstD7G32$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.70))+labs(x = " CumulatedPosition",
                                                                        y = "FST",title = "CD7G32")

###### Save image plot in  PDF file ############
pdf(file = "CD_Fst_position.pdf")

grid.arrange(p1_35,p1_100,p2_35,p2_100, p3_35,p3_100, p4_35,p4_100, p5_32,p5_100, ncol = 2)

dev.off()
############

#############################################################################
############ plot fst along the genome for Androdioecy and Dioecy ancestor agains A6140######################################
#############################################################################

A6140_A00 <- popsync2pooldata("A6140vsA00.sync", poolsizes = c(1000,1000))

A6140_A00.fst<-computeFST(A6140_A00, sliding.window.size=50)

A6140_D00 <- popsync2pooldata("A6140vsD00.sync", poolsizes = c(1000,1000))

A6140_D00.fst<-computeFST(A6140_D00, sliding.window.size=50)

## plotting
p_A00<-ggplot(A6140_A00.fst$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.70))+
  labs(x = "Position",y = "FST",title = "A6140 vs A00")+
  geom_hline(yintercept = A6140_A00.fst$FST, linetype = "dashed", color = "black", size=1)  

p_D00<-ggplot(A6140_D00.fst$sliding.windows.fst, aes(x=CumulatedPosition/1e6, y=MultiLocusFst,group = Chr)) + 
  geom_point(aes(colour = Chr))+coord_cartesian(ylim = c(0, 0.70))+
  labs(x = " Position",y = "",title = "A6140 vs D00")+
  geom_hline(yintercept = A6140_D00.fst$FST, linetype = "dashed", color = "black", size=1) 

grid.arrange( p_A00,p_D00, ncol = 2)


