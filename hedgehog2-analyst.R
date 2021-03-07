#Mackenzie Knox
#BF528 Spring 2021
#Hedgehog
#Project 2 - analyst

library(ggplot2)

df <- read.csv("/projectnb/bf528/users/hedgehog/project_2/programmer/run_cufflinks/cuffdiff/cuffdiff_out/gene_exp.diff", 
                 sep='\t')

#6.1
data <- df[order(df$q_value),]
sub <- data[1:10,c("gene_id", "gene", "value_1", "value_2", "log2.fold_change.", "p_value", "q_value")]


#6.2
hist(data$log2.fold_change., breaks=50, col="purple")
ggplot(data, aes(x=log2.fold_change.)) + 
  geom_histogram(binwidth=.25, fill="plum4") +
  labs(title="Log2 Fold Change (all data)", x="Log2 Fold Change", y="Count")

#6.3
sigdata <- data[which(data$significant=="yes"),]

#6.4
hist(sigdata$log2.fold_change., breaks=50, col="darkblue")
ggplot(sigdata, aes(x=log2.fold_change.)) + 
  geom_histogram(binwidth=.25, fill="lightblue4") +
  labs(title="Log2 Fold Change (significant data)", x="Log2 Fold Change", y="Count")
#bimodal distribution rather than the previous normal distribution
#the counts that were near 0 for log2.fold_change. were non-significant, as the 
#highest count now is <120 rather than 15,000

#6.5
updata <- sigdata[which(sigdata$log2.fold_change. > 0),]
updata <- updata[["gene"]]
updata <- gsub(',', '\n', updata)
downdata <- sigdata[which(sigdata$log2.fold_change. < 0),]
downdata <- downdata[["gene"]]
downdata <- gsub(',', '\n', downdata)
length(updata) #642
length(downdata) #995

#6.6

write.csv(updata, 
          file="/projectnb/bf528/users/hedgehog/project_2/analyst/hedgehog2-upreg.csv", 
          row.names=F, quote=F)
write.csv(downdata, 
          file="/projectnb/bf528/users/hedgehog/project_2/analyst/hedgehog2-downreg.csv", 
          row.names=F, quote=F)
#manually remove the x from the beginning of the csv files
#I can't make it go away
