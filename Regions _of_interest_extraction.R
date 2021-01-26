##CHR 17 SNPS
setwd("~/Documents/Summer-project")
library(ggplot2)
library(sprof)
library(plyr)

chr17SNPS <- read.csv("CHR17SNPs.csv")
View(chr17SNPS)
##Plot the SNPS using SNP vs PHASTCONS
ggplot(data= chr17SNPS, mapping= aes(x=PHASTCONS, y=SNPS)) + geom_point() + labs(title= "SNP count across highly conserved region of chr 17")
#look at most conserved with >75 SNPs
chr17_0.95 <- chr17SNPS[which(chr17SNPS$PHASTCONS >=0.95),] 
chr17_0.95 <- chr17_0.95[which((chr17_0.95$SNPS >=75)),]
View(chr17_0.95)
#convert from list to df to do extraction
chr17_0.95<- list.as.matrix(chr17_0.95, byrow=F)
chr17_0.95<- data.frame(chr17_0.95)
#extract the first region of interest
chr17_0.95_1 <- chr17_0.95[1:3,1:3]
View(chr17_0.95_1)
#find average of the scores/counts over that region
mean(chr17_0.95_1$SNPS)
mean((chr17_0.95_1$PHASTCONS))
#extract the 2nd region of interest
chr17_0.95_2 <- chr17_0.95[6:7,1:3]
View(chr17_0.95_2)
mean(chr17_0.95_2$SNPS)
mean((chr17_0.95_2$PHASTCONS))
#extract the 3rdr egion of interest
chr17_0.95_3 <- chr17_0.95[8:11,1:3]
View(chr17_0.95_3)
mean(chr17_0.95_3$SNPS)
mean((chr17_0.95_3$PHASTCONS))
#extract problem_4 region of interest
chr17_0.95_4 <- chr17_0.95[15:16,1:3]
View(chr17_0.95_4)
mean(chr17_0.95_4$SNPS)
mean((chr17_0.95_4$PHASTCONS))

##extract chr17 0.90-0.95 regions of interest
chr17_0.90 <- chr17SNPS[which(chr17SNPS$PHASTCONS >=0.90),]
chr17_0.90 <-chr17_0.90[which(chr17_0.90$PHASTCONS <0.95),]
ggplot(data=chr17_0.90, mapping = aes(x= PHASTCONS, y=SNPS)) + geom_point()
View(chr17_0.90)
chr17_0.90 <-chr17_0.90[which(chr17_0.90$SNPS >=80),]
chr17_0.90 <- list.as.matrix(chr17_0.90)
chr17_0.90 <- data.frame(chr17_0.90)
#extract the 1st region of interest
chr17_0.90_1 <- chr17_0.90[3:4,1:3]
View(chr17_0.90_1)
mean(chr17_0.90_1$SNPS)
mean((chr17_0.90_1$PHASTCONS))
#extract the 2nd region of interest
chr17_0.90_2 <- chr17_0.90[5:6,1:3]
View(chr17_0.90_2)
mean(chr17_0.90_2$SNPS)
mean((chr17_0.90_2$PHASTCONS))

##Chr17 with outlier SNPs coloured
ggplot(data= chr17SNPS, mapping= aes(x=PHASTCONS, y=SNPS)) + geom_point() + labs(title= "SNP count across highly conserved region of chr 17")
chr17SNPS
plot17<-ggplot(data=chr17SNPS, aes(x=PHASTCONS, y=SNPS)) + geom_point()+geom_point(data= sub1795, colour="red")
plot17
sub1795<- subset(chr17SNPS, Start %in% c(117379, 117399,117419,117839,118139,7484359,7484379,7571819,7571839,7571859,7571879,43100579,46771679,48593799,78760439,78760459))
sub1795
sub1790<- subset(chr17SNPS, Start %in% c(117339,117359,	117799,117819,118099,118119,7484339,8186879,38867639,39153739,43100599,57682279,59255779,	71610319,78760419))
sub1790
plot17<-ggplot(chr17SNPS, aes(x=PHASTCONS, y=SNPS)) + geom_point(alpha=0.3) +geom_point(data=sub1795, colour="red") +geom_point(data=sub1790, colour="blue")+ labs (title = "SNPs count across conserved regions of chromosome 17")
plot17

#### Chr21SNPS
chr21SNPs <- read.csv("OutputCHR21SNPs.csv", header = TRUE)
#colnames(chr21SNPs) <- c("Start", "SNPs", "PHASTCONS")
View(chr21SNPs)
ggplot(data= chr21SNPs, mapping= aes(x=PHASTCONS, y=SNPS)) +geom_point() + labs (title = "SNPs across conserved regions on Chr21")
ggplot(data = chr21SNPs, mapping = aes(x=SNPs,)) +geom_bar(colour="black", fill="white") +labs(title = "SNP density across Chr21")
Prob_21_95 <- chr21SNPs[which(chr21SNPs$PHASTCONS >=0.95),]
View(Prob_21_95)
Prob_21_95 <- Prob_21_95[which(Prob_21_95$SNPS >= 75),]
Prob_21_95 <- list.as.matrix(Prob_21_95)
Prob_21_95 <- data.frame(Prob_21_95)
#extract first region of interesr
p21_95_1 <- Prob_21_95[1:2, 1:3]
View(p21_95_1)
mean (p21_95_1$PHASTCONS)
mean(p21_95_1$SNPS)
#extract 2nd region of intereest
p21_95_2 <- Prob_21_95[5:16, 1:3]
View(p21_95_2)
mean (p21_95_2$PHASTCONS)
mean(p21_95_2$SNPS)
#extract 3rd region of interest
p21_95_3 <- Prob_21_95[17:18, 1:3]
View(p21_95_3)
mean (p21_95_3$PHASTCONS)
mean(p21_95_3$SNPS)
#extract 4th region of interest
p21_95_4 <- Prob_21_95[19:26, 1:3]
View(p21_95_4)
mean (p21_95_4$PHASTCONS)
mean(p21_95_4$SNPS)
#extract 5th region of interest
p21_95_5 <- Prob_21_95[35:49, 1:3]
View(p21_95_5)
mean (p21_95_5$PHASTCONS)
mean(p21_95_5$SNPS)
#extarct 6th region of interest
p21_95_6 <- Prob_21_95[50:51, 1:3]
View(p21_95_6)
mean (p21_95_6$PHASTCONS)
mean(p21_95_6$SNPS)
#extract 7th region of interest
p21_95_7 <- Prob_21_95[56:60, 1:3]
View(p21_95_7)
mean (p21_95_7$PHASTCONS)
mean(p21_95_7$SNPS)
#extract 8th region of interest
p21_95_8 <- Prob_21_95[61:73, 1:3]
View(p21_95_8)
mean (p21_95_8$PHASTCONS)
mean(p21_95_8$SNPS)
#extract 9th region of interest
p21_95_9 <- Prob_21_95[76:77, 1:3]
View(p21_95_9)
mean (p21_95_9$PHASTCONS)
mean(p21_95_9$SNPS)
#extract 10th region of interest
p21_95_10 <- Prob_21_95[78:81, 1:3]
View(p21_95_10)
mean (p21_95_10$PHASTCONS)
mean(p21_95_10$SNPS)

##extract 0.90-0.95 region for chr 21
Prob_21_90 <- chr21SNPs[which(chr21SNPs$PHASTCONS >=0.90),]
Prob_21_90 <- Prob_21_90[which(Prob_21_90$PHASTCONS <0.95),]
View(Prob_21_90)
Prob_21_90 <- Prob_21_90[which(Prob_21_90$SNPS >= 100),]
Prob_21_90 <- list.as.matrix(Prob_21_90)
Prob_21_90 <- data.frame(Prob_21_90)
#extract 1st region of interest
p21_90_1 <- Prob_21_90[1:4,1:3]
View(p21_90_1)
mean(p21_90_1$PHASTCONS)
mean(p21_90_1$SNPS)
#extract 12nd region of interest
p21_90_2 <- Prob_21_90[8:10,1:3]
View(p21_90_2)
mean(p21_90_2$PHASTCONS)
mean(p21_90_2$SNPS)

##chromosome 21 coloured graph
ggplot(data= chr21SNPs, mapping= aes(x=PHASTCONS, y=SNPS)) + geom_point() + labs(title= "SNP count across highly conserved region of chr21")
chr17SNPS
plot21<-ggplot(data=chr21SNPs, aes(x=PHASTCONS, y=SNPS)) + geom_point()+geom_point(data= sub2195, colour="red")
plot21
sub2195 <- subset(chr21SNPs, Start %in% c(8777600,8777620,8990520,8990540,8990640,8990660,8990680,8990700,	8990720,8990740,	8990760,8990780,8990800,	8990820,	8990840,8990860,9091160,9091180,9092780,9092800,9092820,9092840,9092860,9092880,9092900,9093040,9093280,9096240,9908100,9908680,10330740,10330760,10330820,10330840,	10330980,10331000,10331020,10331040,10331220,10331240,10331260,10331280,10331300,10331320,10331340,	10331360,10331380,10331420,10331440,10364220,10364240,	10379400,	10735260,10735280,10735300,10735680,10735700,	10735720,10735740,10735760,10736600,10736620,10736640,10736660,10737000,10737020,10737040,10737060,10737080,10737120,10737140,10737160,10737180,10776620,10776640,36004940,	36004960,	41475120,41475140,41475160,	41475180))
sub2195
sub2190<- subset(chr21SNPs, Start %in% c(8990560,8990580,	8990600,8990620,8990880,9092760,41475100,41475260,41475280,41475300))
sub2190
plot21<-ggplot(chr21SNPs, aes(x=PHASTCONS, y=SNPS)) + geom_point(alpha=0.3) +geom_point(data=sub2195, colour="red") +geom_point(data=sub2190, colour="blue")+ labs (title = "SNPs across conserved regions of chromosome 21")
plot21

##CHR6
CHR6SNPS <- read.csv("OutCHR6SNPs.csv")
View(CHR6SNPS)
ggplot(data = CHR6SNPS, mapping= aes(x=PHASTCONS, y=SNPS)) + geom_point() + labs (title = "SNP count across highly conserved regions of Chromosome 6")
#extract 0.95-0.1 regions of interest for chr 6
Prob_6_95 <- CHR6SNPS[which(CHR6SNPS$PHASTCONS >=0.95),]
View(Prob_6_95)
Prob_6_95 <- Prob_6_95[which (Prob_6_95$SNPS >=70),]
Prob_6_95 <- list.as.matrix(Prob_6_95)
Prob_6_95<- data.frame(Prob_6_95)
#extract 1st region of interest
p6_95_1 <- Prob_6_95[2:3,1:3]
View(p6_95_1)
mean(p6_95_1$SNPS)
mean(p6_95_1$PHASTCONS)
#extract 2nd region of interest 
p6_95_2 <- Prob_6_95[4:5,1:3]
View(p6_95_2)
mean(p6_95_2$SNPS)
mean(p6_95_2$PHASTCONS)
#extract 3rd region of interest
p6_95_3 <- Prob_6_95[6:10,1:3]
View(p6_95_3)
mean(p6_95_3$SNPS)
mean(p6_95_3$PHASTCONS)
#extract 4th region of interest
p6_95_4 <- Prob_6_95[12:13,1:3]
View(p6_95_4)
mean(p6_95_4$SNPS)
mean(p6_95_4$PHASTCONS)
#extract 5th region of interest
p6_95_5<- Prob_6_95[16:18,1:3]
View(p6_95_5)
mean(p6_95_5$SNPS)
mean(p6_95_5$PHASTCONS)
#ectract 6th region of interest
p6_95_6 <- Prob_6_95[19:24,1:3]
View(p6_95_6)
mean(p6_95_6$SNPS)
mean(p6_95_6$PHASTCONS)

#chr6 0.90-0.95 regions of interest
Prob_6_90 <- CHR6SNPS[which(CHR6SNPS$PHASTCONS >=0.90),]
Prob_6_90 <- Prob_6_90[which(Prob_6_90$PHASTCONS <=0.95),]
View(Prob_6_90)
Prob_6_90 <- Prob_6_90[which(Prob_6_90$SNPS >=85),]
Prob_6_90 <- list.as.matrix(Prob_6_90)
Prob_6_90<- data.frame(Prob_6_90)
#extract 1st region of interest
p6_90_1 <- Prob_6_90[7:8, 1:3]
View(p6_90_1)
mean(p6_90_1$SNPS)
mean(p6_90_1$PHASTCONS)
#extract 2nd region of interst
p6_90_2 <- Prob_6_90[9:10, 1:3]
View(p6_90_2)
mean(p6_90_2$SNPS)
mean(p6_90_2$PHASTCONS)

##CHr6 with outlier SNPs coloured
ggplot(data= CHR6SNPS, mapping= aes(x=PHASTCONS, y=SNPS)) + geom_point() + labs(title= "SNP count across conserved regions of Chromosome 6")
plot6<-ggplot(data=CHR6SNPS, aes(x=PHASTCONS, y=SNPS)) + geom_point()+geom_point(data= sub695, colour="red")
plot6
sub695<- subset(CHR6SNPS, Start %in% c(6006945,42782340,42782355,52995705,52995720,54770745,54770760,54770775,54770790,54770805,85290255,90296130,90296145,106083195,139374585,156662490,156662505,156662520,168615060,168615075,168615090,168615105,168615120,168615135))
sub695
sub690<- subset(CHR6SNPS, Start %in% c(6006930,26328135,27177210,27777885,28944570,42782325,89339535,89339550,90296115,90296160))
sub690
plot6<-ggplot(CHR6SNPS, aes(x=PHASTCONS, y=SNPS)) + geom_point(alpha=0.3) +geom_point(data=sub695, colour="red") +geom_point(data=sub690, colour="blue")+ labs (title = "SNPs count across conserved regions of chromosome 6")
plot6

