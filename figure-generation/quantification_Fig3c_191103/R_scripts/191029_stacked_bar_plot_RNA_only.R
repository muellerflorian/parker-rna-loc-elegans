#####
#This script plots P granule and RNA clustering data in a series of stacked bar plots
#####


#####
#Open libraries
#####

library(plyr)
library(dplyr)
library(ggplot2)


#####
#Open P granule and RNA cluster data, append P granule Falses to RNA.clusters data frame, and normalize data set #1
#####
getwd()
chs1.RNA.clusters1 <- read.csv(file="../chs-1/Image_01_RNA_clusters.csv", header = TRUE, sep=";")
as.data.frame(chs1.RNA.clusters1)
#Opens RNA cluster data file and reads as data frame

chs1.RNA.clusters1$cell <- rep("Image_01",nrow(chs1.RNA.clusters1))
#Appends P granule falses to RNA cluster data and adds a column named after image number.

chs1.values1 <- table(chs1.RNA.clusters1$coloc)
#Gets number of True, false, and P granule false.

chs1.false1 <- chs1.values1["False"]
chs1.true1 <- chs1.values1["True"]
#Separates by colocalization type

chs1.total1 <- sum(chs1.false1, chs1.true1)
#Sums all these values

chs1.normalized.false1 <- chs1.false1/chs1.total1
chs1.normalized.true1 <- chs1.true1/chs1.total1
#Normalize

chs1.normalized.plot1 <- data.frame("normalized_values" = c(chs1.normalized.false1, chs1.normalized.true1), "cell" = as.character(rep("Image_01")), "coloc" = c("False","True"))
#makes data frame for normalized data

#####
#Repeat first block with new data
#####

chs1.RNA.clusters2 <- read.csv(file="../chs-1/Image_05_RNA_clusters.csv", header = TRUE, sep=";")
as.data.frame(chs1.RNA.clusters2)

chs1.RNA.clusters2$cell <- rep("Image_02",nrow(chs1.RNA.clusters2))

chs1.values2 <- table(chs1.RNA.clusters2$coloc)
chs1.false2 <- chs1.values2["False"]
chs1.true2 <- chs1.values2["True"]
chs1.total2 <- sum(chs1.false2, chs1.true2)
chs1.normalized.false2 <- chs1.false2/chs1.total2
chs1.normalized.true2 <- chs1.true2/chs1.total2

chs1.normalized.plot2 <- data.frame("normalized_values" = c(chs1.normalized.false2, chs1.normalized.true2), "cell" = as.character(rep("Image_02")), "coloc" = c("False","True"))

#####
#Repeat first block with new data
#####

chs1.RNA.clusters3 <- read.csv(file="../chs-1/Image_07_RNA_clusters.csv", header = TRUE, sep=";")
as.data.frame(chs1.RNA.clusters3)

chs1.RNA.clusters3$cell <- rep("Image_03",nrow(chs1.RNA.clusters3))

chs1.values3 <- table(chs1.RNA.clusters3$coloc)
chs1.false3 <- chs1.values3["False"]
chs1.true3 <- chs1.values3["True"]
chs1.total3 <- sum(chs1.false3, chs1.true3)
chs1.normalized.false3 <- chs1.false3/chs1.total3
chs1.normalized.true3 <- chs1.true3/chs1.total3

chs1.normalized.plot3 <- data.frame("normalized_values" = c(chs1.normalized.false3, chs1.normalized.true3), "cell" = as.character(rep("Image_03")), "coloc" = c("False","True"))

#####
#Repeat first block with new data
#####

chs1.RNA.clusters4 <- read.csv(file="../chs-1/Image_18_RNA_clusters.csv", header = TRUE, sep=";")
as.data.frame(chs1.RNA.clusters4)

chs1.RNA.clusters4$cell <- rep("Image_04",nrow(chs1.RNA.clusters4))

chs1.values4 <- table(chs1.RNA.clusters4$coloc)
chs1.false4 <- chs1.values4["False"]
chs1.true4 <- chs1.values4["True"]
chs1.total4 <- sum(chs1.false4, chs1.true4)
chs1.normalized.false4 <- chs1.false4/chs1.total4
chs1.normalized.true4 <- chs1.true4/chs1.total4

chs1.normalized.plot4 <- data.frame("normalized_values" = c(chs1.normalized.false4, chs1.normalized.true4), "cell" = as.character(rep("Image_04")), "coloc" = c("False","True"))


#####
#Synthesize data frames
#####

chs1.allclusters <- rbind_list(chs1.RNA.clusters1, chs1.RNA.clusters2, chs1.RNA.clusters3, chs1.RNA.clusters4)
chs1.normalized.plot <- rbind_list(chs1.normalized.plot1, chs1.normalized.plot2, chs1.normalized.plot3, chs1.normalized.plot4)
#Creates large data frames from all the different sets for both raw and normalized data

#####
#Make normalized averages
#####

chs1.normalized.total.true <- chs1.normalized.plot %>% filter(coloc == "True")
chs1.normalized.average.true <- sum(chs1.normalized.total.true$normalized_values)/4
chs1.normalized.total.false <- chs1.normalized.plot %>% filter(coloc == "False")
chs1.normalized.average.false <- sum(chs1.normalized.total.false$normalized_values)/4
chs1.normalized.average <- data.frame("normalized_values" = c(chs1.normalized.average.false, chs1.normalized.average.true), "cell" = as.character(rep("Norm_Avg")), "coloc" = c("False","True"))


#####
#Make a combined raw data data frame
#####

chs1.observations.true <- chs1.allclusters %>% filter(coloc == "True")
chs1.observations.false <- chs1.allclusters %>% filter(coloc == "False")
chs1.alldata <- rbind(chs1.observations.true, chs1.observations.false)
chs1.alldata$avg <- rep("Avg",nrow(chs1.alldata))

as.data.frame(chs1.alldata)[1:100,]
head(chs1.alldata)
str(chs1.alldata)
dim(chs1.alldata)

chs.df <- chs1.alldata %>% 
  mutate(transcript = "chs-1")
#####
#Repeat everything for nos-2 RNA
#####

nos2.RNA.clusters1 <- read.csv(file="../nos-2/Image_01_RNA_clusters.csv", header = TRUE, sep=";")
as.data.frame(nos2.RNA.clusters1)

nos2.RNA.clusters1$cell <- rep("Image_01",nrow(nos2.RNA.clusters1))

nos2.values1 <- table(nos2.RNA.clusters1$coloc)
nos2.false1 <- nos2.values1["False"]
nos2.true1 <- nos2.values1["True"]
nos2.total1 <- sum(nos2.false1, nos2.true1)
nos2.normalized.false1 <- nos2.false1/nos2.total1
nos2.normalized.true1 <- nos2.true1/nos2.total1

nos2.normalized.plot1 <- data.frame("normalized_values" = c(nos2.normalized.false1, nos2.normalized.true1), "cell" = as.character(rep("Image_01")), "coloc" = c("False","True"))

#####
#Repeat first block with new data
#####

nos2.RNA.clusters2 <- read.csv(file="../nos-2/Image_05_RNA_clusters.csv", header = TRUE, sep=";")
as.data.frame(nos2.RNA.clusters2)

nos2.RNA.clusters2$cell <- rep("Image_02",nrow(nos2.RNA.clusters2))

nos2.values2 <- table(nos2.RNA.clusters2$coloc)
nos2.false2 <- nos2.values2["False"]
nos2.true2 <- nos2.values2["True"]
nos2.total2 <- sum(nos2.false2, nos2.true2)
nos2.normalized.false2 <- nos2.false2/nos2.total2
nos2.normalized.true2 <- nos2.true2/nos2.total2

nos2.normalized.plot2 <- data.frame("normalized_values" = c(nos2.normalized.false2, nos2.normalized.true2), "cell" = as.character(rep("Image_02")), "coloc" = c("False","True"))

#####
#Repeat first block with new data
#####

nos2.RNA.clusters3 <- read.csv(file="../nos-2/Image_07_RNA_clusters.csv", header = TRUE, sep=";")
as.data.frame(nos2.RNA.clusters3)

nos2.RNA.clusters3$cell <- rep("Image_03",nrow(nos2.RNA.clusters3))

nos2.values3 <- table(nos2.RNA.clusters3$coloc)
nos2.false3 <- nos2.values3["False"]
nos2.true3 <- nos2.values3["True"]
nos2.total3 <- sum(nos2.false3, nos2.true3)
nos2.normalized.false3 <- nos2.false3/nos2.total3
nos2.normalized.true3 <- nos2.true3/nos2.total3

nos2.normalized.plot3 <- data.frame("normalized_values" = c(nos2.normalized.false3, nos2.normalized.true3), "cell" = as.character(rep("Image_03")), "coloc" = c("False","True"))

#####
#Repeat first block with new data
#####

nos2.RNA.clusters4 <- read.csv(file="../nos-2/Image_18_RNA_clusters.csv", header = TRUE, sep=";")
as.data.frame(nos2.RNA.clusters4)

nos2.RNA.clusters4$cell <- rep("Image_04",nrow(nos2.RNA.clusters4))

nos2.values4 <- table(nos2.RNA.clusters4$coloc)
nos2.false4 <- nos2.values4["False"]
nos2.true4 <- nos2.values4["True"]
nos2.total4 <- sum(nos2.false4, nos2.true4)
nos2.normalized.false4 <- nos2.false4/nos2.total4
nos2.normalized.true4 <- nos2.true4/nos2.total4

nos2.normalized.plot4 <- data.frame("normalized_values" = c(nos2.normalized.false4, nos2.normalized.true4), "cell" = as.character(rep("Image_04")), "coloc" = c("False","True"))


#####
#Synthesize data frames
#####

nos2.allclusters <- rbind_list(nos2.RNA.clusters1, nos2.RNA.clusters2, nos2.RNA.clusters3, nos2.RNA.clusters4)
nos2.normalized.plot <- rbind_list(nos2.normalized.plot1, nos2.normalized.plot2, nos2.normalized.plot3, nos2.normalized.plot4)
#Creates large data frames from all the different sets for both raw and normalized data

#####
#Make normalized averages
#####

nos2.normalized.total.true <- nos2.normalized.plot %>% filter(coloc == "True")
nos2.normalized.average.true <- sum(nos2.normalized.total.true$normalized_values)/4
nos2.normalized.total.false <- nos2.normalized.plot %>% filter(coloc == "False")
nos2.normalized.average.false <- sum(nos2.normalized.total.false$normalized_values)/4
nos2.normalized.average <- data.frame("normalized_values" = c(nos2.normalized.average.false, nos2.normalized.average.true), "cell" = as.character(rep("Norm_Avg")), "coloc" = c("False","True"))


#####
#Make a combined raw data data frame
#####

nos2.observations.true <- nos2.allclusters %>% filter(coloc == "True")
nos2.observations.false <- nos2.allclusters %>% filter(coloc == "False")
nos2.alldata <- rbind(nos2.observations.true, nos2.observations.false)
nos2.alldata$avg <- rep("Avg",nrow(nos2.alldata))


nos.df <- nos2.alldata %>% 
  mutate(transcript = "nos-2")
nos.df

#####
#Repeat everything for clu-1 RNA
#####

clu1.RNA.clusters1 <- read.csv(file="../clu-1/Image_02_RNA_clusters.csv", header = TRUE, sep=";")
as.data.frame(clu1.RNA.clusters1)

clu1.RNA.clusters1$cell <- rep("Image_02",nrow(clu1.RNA.clusters1))

clu1.values1 <- table(clu1.RNA.clusters1$coloc)
clu1.false1 <- clu1.values1["False"]
clu1.true1 <- clu1.values1["True"]
clu1.total1 <- sum(clu1.false1, clu1.true1)
clu1.normalized.false1 <- clu1.false1/clu1.total1
clu1.normalized.true1 <- clu1.true1/clu1.total1

clu1.normalized.plot1 <- data.frame("normalized_values" = c(clu1.normalized.false1, clu1.normalized.true1), "cell" = as.character(rep("Image_01")), "coloc" = c("False","True"))

#####
#Repeat first block with new data
#####

clu1.RNA.clusters2 <- read.csv(file="../clu-1/Image_06_RNA_clusters.csv", header = TRUE, sep=";")
as.data.frame(clu1.RNA.clusters2)

clu1.RNA.clusters2$cell <- rep("Image_06",nrow(clu1.RNA.clusters2))

clu1.values2 <- table(clu1.RNA.clusters2$coloc)
clu1.false2 <- clu1.values2["False"]
clu1.true2 <- clu1.values2["True"]
clu1.total2 <- sum(clu1.false2, clu1.true2)
clu1.normalized.false2 <- clu1.false2/clu1.total2
clu1.normalized.true2 <- clu1.true2/clu1.total2

clu1.normalized.plot2 <- data.frame("normalized_values" = c(clu1.normalized.false2, clu1.normalized.true2), "cell" = as.character(rep("Image_02")), "coloc" = c("False","True"))

#####
#Repeat first block with new data
#####

clu1.RNA.clusters3 <- read.csv(file="../clu-1/Image_06-2_RNA_clusters.csv", header = TRUE, sep=";")
as.data.frame(clu1.RNA.clusters3)

clu1.RNA.clusters3$cell <- rep("Image_06-2",nrow(clu1.RNA.clusters3))

clu1.values3 <- table(clu1.RNA.clusters3$coloc)
clu1.false3 <- clu1.values3["False"]
clu1.true3 <- clu1.values3["True"]
clu1.total3 <- sum(clu1.false3, clu1.true3)
clu1.normalized.false3 <- clu1.false3/clu1.total3
clu1.normalized.true3 <- clu1.true3/clu1.total3

clu1.normalized.plot3 <- data.frame("normalized_values" = c(clu1.normalized.false3, clu1.normalized.true3), "cell" = as.character(rep("Image_03")), "coloc" = c("False","True"))

#####
#Repeat first block with new data
#####

clu1.RNA.clusters4 <- read.csv(file="../clu-1/Image_10_RNA_clusters.csv", header = TRUE, sep=";")
as.data.frame(clu1.RNA.clusters4)

clu1.RNA.clusters4$cell <- rep("Image_10",nrow(clu1.RNA.clusters4))

clu1.values4 <- table(clu1.RNA.clusters4$coloc)
clu1.false4 <- clu1.values4["False"]
clu1.true4 <- clu1.values4["True"]
clu1.total4 <- sum(clu1.false4, clu1.true4)
clu1.normalized.false4 <- clu1.false4/clu1.total4
clu1.normalized.true4 <- clu1.true4/clu1.total4

clu1.normalized.plot4 <- data.frame("normalized_values" = c(clu1.normalized.false4, clu1.normalized.true4), "cell" = as.character(rep("Image_04")), "coloc" = c("False","True"))


#####
#Synthesize data frames
#####

clu1.allclusters <- rbind_list(clu1.RNA.clusters1, clu1.RNA.clusters2, clu1.RNA.clusters3, clu1.RNA.clusters4)
clu1.normalized.plot <- rbind_list(clu1.normalized.plot1, clu1.normalized.plot2, clu1.normalized.plot3, clu1.normalized.plot4)
#Creates large data frames from all the different sets for both raw and normalized data

#####
#Make normalized averages
#####

clu1.normalized.total.true <- clu1.normalized.plot %>% filter(coloc == "True")
clu1.normalized.average.true <- sum(clu1.normalized.total.true$normalized_values)/4
clu1.normalized.total.false <- clu1.normalized.plot %>% filter(coloc == "False")
clu1.normalized.average.false <- sum(clu1.normalized.total.false$normalized_values)/4
clu1.normalized.average <- data.frame("normalized_values" = c(clu1.normalized.average.false, clu1.normalized.average.true), "cell" = as.character(rep("Norm_Avg")), "coloc" = c("False","True"))


#####
#Make a combined raw data data frame
#####

clu1.observations.true <- clu1.allclusters %>% filter(coloc == "True")
clu1.observations.false <- clu1.allclusters %>% filter(coloc == "False")
clu1.alldata <- rbind(clu1.observations.true, clu1.observations.false)
clu1.alldata$avg <- rep("Avg",nrow(clu1.alldata))


clu.df <- clu1.alldata %>% 
  mutate(transcript = "clu-1")
clu.df

#####
#Repeat everything for cpg-2 RNA
#####

cpg2.RNA.clusters1 <- read.csv(file="../cpg-2/Image_02_RNA_clusters.csv", header = TRUE, sep=";")
as.data.frame(cpg2.RNA.clusters1)

cpg2.RNA.clusters1$cell <- rep("Image_02",nrow(cpg2.RNA.clusters1))

cpg2.values1 <- table(cpg2.RNA.clusters1$coloc)
cpg2.false1 <- cpg2.values1["False"]
cpg2.true1 <- cpg2.values1["True"]
cpg2.total1 <- sum(cpg2.false1, cpg2.true1)
cpg2.normalized.false1 <- cpg2.false1/cpg2.total1
cpg2.normalized.true1 <- cpg2.true1/cpg2.total1

cpg2.normalized.plot1 <- data.frame("normalized_values" = c(cpg2.normalized.false1, cpg2.normalized.true1), "cell" = as.character(rep("Image_01")), "coloc" = c("False","True"))

#####
#Repeat first block with new data
#####

cpg2.RNA.clusters2 <- read.csv(file="../cpg-2/Image_06_RNA_clusters.csv", header = TRUE, sep=";")
as.data.frame(cpg2.RNA.clusters2)

cpg2.RNA.clusters2$cell <- rep("Image_06",nrow(cpg2.RNA.clusters2))

cpg2.values2 <- table(cpg2.RNA.clusters2$coloc)
cpg2.false2 <- cpg2.values2["False"]
cpg2.true2 <- cpg2.values2["True"]
cpg2.total2 <- sum(cpg2.false2, cpg2.true2)
cpg2.normalized.false2 <- cpg2.false2/cpg2.total2
cpg2.normalized.true2 <- cpg2.true2/cpg2.total2

cpg2.normalized.plot2 <- data.frame("normalized_values" = c(cpg2.normalized.false2, cpg2.normalized.true2), "cell" = as.character(rep("Image_02")), "coloc" = c("False","True"))

#####
#Repeat first block with new data
#####

cpg2.RNA.clusters3 <- read.csv(file="../cpg-2/Image_06-2_RNA_clusters.csv", header = TRUE, sep=";")
as.data.frame(cpg2.RNA.clusters3)

cpg2.RNA.clusters3$cell <- rep("Image_06-2",nrow(cpg2.RNA.clusters3))

cpg2.values3 <- table(cpg2.RNA.clusters3$coloc)
cpg2.false3 <- cpg2.values3["False"]
cpg2.true3 <- cpg2.values3["True"]
cpg2.total3 <- sum(cpg2.false3, cpg2.true3)
cpg2.normalized.false3 <- cpg2.false3/cpg2.total3
cpg2.normalized.true3 <- cpg2.true3/cpg2.total3

cpg2.normalized.plot3 <- data.frame("normalized_values" = c(cpg2.normalized.false3, cpg2.normalized.true3), "cell" = as.character(rep("Image_03")), "coloc" = c("False","True"))

#####
#Repeat first block with new data
#####

cpg2.RNA.clusters4 <- read.csv(file="../cpg-2/Image_10_RNA_clusters.csv", header = TRUE, sep=";")
as.data.frame(cpg2.RNA.clusters4)

cpg2.RNA.clusters4$cell <- rep("Image_10",nrow(cpg2.RNA.clusters4))

cpg2.values4 <- table(cpg2.RNA.clusters4$coloc)
cpg2.false4 <- cpg2.values4["False"]
cpg2.true4 <- cpg2.values4["True"]
cpg2.total4 <- sum(cpg2.false4, cpg2.true4)
cpg2.normalized.false4 <- cpg2.false4/cpg2.total4
cpg2.normalized.true4 <- cpg2.true4/cpg2.total4

cpg2.normalized.plot4 <- data.frame("normalized_values" = c(cpg2.normalized.false4, cpg2.normalized.true4), "cell" = as.character(rep("Image_04")), "coloc" = c("False","True"))


#####
#Synthesize data frames
#####

cpg2.allclusters <- rbind_list(cpg2.RNA.clusters1, cpg2.RNA.clusters2, cpg2.RNA.clusters3, cpg2.RNA.clusters4)
cpg2.normalized.plot <- rbind_list(cpg2.normalized.plot1, cpg2.normalized.plot2, cpg2.normalized.plot3, cpg2.normalized.plot4)
#Creates large data frames from all the different sets for both raw and normalized data

#####
#Make normalized averages
#####

cpg2.normalized.total.true <- cpg2.normalized.plot %>% filter(coloc == "True")
cpg2.normalized.average.true <- sum(cpg2.normalized.total.true$normalized_values)/4
cpg2.normalized.total.false <- cpg2.normalized.plot %>% filter(coloc == "False")
cpg2.normalized.average.false <- sum(cpg2.normalized.total.false$normalized_values)/4
cpg2.normalized.average <- data.frame("normalized_values" = c(cpg2.normalized.average.false, cpg2.normalized.average.true), "cell" = as.character(rep("Norm_Avg")), "coloc" = c("False","True"))


#####
#Make a combined raw data data frame
#####

cpg2.observations.true <- cpg2.allclusters %>% filter(coloc == "True")
cpg2.observations.false <- cpg2.allclusters %>% filter(coloc == "False")
cpg2.alldata <- rbind(cpg2.observations.true, cpg2.observations.false)
cpg2.alldata$avg <- rep("Avg",nrow(cpg2.alldata))


cpg.df <- cpg2.alldata %>% 
  mutate(transcript = "cpg-2")
cpg.df


### Merge dataframes together and export as source data:
fulldata.df <- rbind(nos.df, chs.df, cpg.df, clu.df)
fulldata.df <- as.data.frame(fulldata.df)

write.table(fulldata.df, file = "../figure3_sourcedata1_191103.txt", sep = "\t", quote = FALSE)


#####
#Plot normalized averages
#SHOWN IN FIGURE 3C!!!
#####

chs1.normal.average <- ggplot(chs1.normalized.average, aes(x=cell))
chs1.normal.average + geom_bar(stat="identity", aes(fill = coloc, y = normalized_values))

nos2.normal.average <- ggplot(nos2.normalized.average, aes(x=cell))
nos2.normal.average + geom_bar(stat="identity", aes(fill = coloc, y = normalized_values))

clu1.normal.average <- ggplot(clu1.normalized.average, aes(x=cell))
clu1.normal.average + geom_bar(stat="identity", aes(fill = coloc, y = normalized_values))

cpg2.normal.average <- ggplot(cpg2.normalized.average, aes(x=cell))
cpg2.normal.average + geom_bar(stat="identity", aes(fill = coloc, y = normalized_values))

#####
#Plot total combined observations
#DATA NOT SHOWN
#####

chs1.total.observations <- ggplot(chs1.alldata, aes(x=avg))
chs1.total.observations + geom_bar( aes(fill = coloc))

nos2.total.observations <- ggplot(nos2.alldata, aes(x=avg))
nos2.total.observations + geom_bar( aes(fill = coloc))

clu1.total.observations <- ggplot(clu1.alldata, aes(x=avg))
clu1.total.observations + geom_bar( aes(fill = coloc))

cpg2.total.observations <- ggplot(cpg2.alldata, aes(x=avg))
cpg2.total.observations + geom_bar( aes(fill = coloc))

#####
#Plot raw data by cell
#DATA NOT SHOWN
#####

chs1.raw.data <- ggplot(chs1.allclusters, aes(x=cell))
chs1.raw.data + geom_bar(aes(fill = coloc))

nos2.raw.data <- ggplot(nos2.allclusters, aes(x=cell))
nos2.raw.data + geom_bar(aes(fill = coloc))

clu1.raw.data <- ggplot(clu1.allclusters, aes(x=cell))
clu1.raw.data + geom_bar(aes(fill = coloc))

cpg2.raw.data <- ggplot(cpg2.allclusters, aes(x=cell))
cpg2.raw.data + geom_bar(aes(fill = coloc))

#####
#Plot normalized data by cell
#DATA NOT SHOWN
#####

chs1.embryo.normalized.data <- ggplot(chs1.normalized.plot, aes(x=cell))
chs1.embryo.normalized.data + geom_bar(stat="identity", aes(fill = coloc, y = normalized_values))

nos2.embryo.normalized.data <- ggplot(nos2.normalized.plot, aes(x=cell))
nos2.embryo.normalized.data + geom_bar(stat="identity", aes(fill = coloc, y = normalized_values))

clu1.embryo.normalized.data <- ggplot(clu1.normalized.plot, aes(x=cell))
clu1.embryo.normalized.data + geom_bar(stat="identity", aes(fill = coloc, y = normalized_values))

cpg2.embryo.normalized.data <- ggplot(cpg2.normalized.plot, aes(x=cell))
cpg2.embryo.normalized.data + geom_bar(stat="identity", aes(fill = coloc, y = normalized_values))


