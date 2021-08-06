
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")
BiocManager::install("GenomicAlignments")
library("Biostrings")
##IMPORTANT##########################
###Make sure to switch the lables of the model if you swapped the labels during the trainig!
###########################################################
data_p = read.csv("~/myfutdata_NA.csv",sep= ",",header=T,stringsAsFactors=FALSE)
data_p = data_p[,2:ncol(data_p)]
length(which(data_p$outcome==1))/nrow(data_p)
fut_tip = data_p$tip
data_p = data_p[,-80]
data_p = scaleData(data_p)
data_predict = data_p
#the best model
SVM_pred=predict(svm.fit,data_predict)
successful_tips_ind  = which(SVM_pred == 1)
length(successful_tips_ind )
successful_tips = fut_tip[successful_tips_ind]
#read the epitope distance of the tips
epi_Dis = read.csv("~/epitopeDis.csv",sep= ",",header=T,stringsAsFactors=FALSE)
epi_Dis = epi_Dis[,2:ncol(epi_Dis)]
tip_ind = match(successful_tips,epi_Dis[,3])
length(which(is.na(tip_ind)))
epi_Dis_tips = epi_Dis[tip_ind ,]
plot(density(epi_Dis_tips[,2]))
density(epi_Dis_tips[,2])
median(epi_Dis_tips[,2])
#===================================================================================
#===================================================================================
#we choose our vaccine among the set of sequences after 2019/1
L=sapply(epi_Dis_tips$tip_name , function(x) strsplit(x,"--"))
head(L)
Date = sapply(L, function(x) x[length(x)])
dd = as.Date(Date,"%m/%d/%Y")
ind = which(dd>="2019-1-1")

fut_table = epi_Dis_tips[ind,]
plot(density(fut_table[,2]))
density(fut_table[,2])
median(fut_table[,2])

vaccine_tip = which(fut_table[,2] >= 6.598476 -0.01 & fut_table[,2] <= 6.598476 +0.01)
vaccine_tip_name = fut_table[vaccine_tip,3]
vaccine_tip_name
#===================================================================================
source("~/VaccineDist.R")
dis=numeric()
for(mytip in vaccine_tip_name){
  tipinf= getVaccineDist_new(mytip, Aux_data,hdata, Pdata,Pdata_2020, year = 2020)
  dis = c(dis,tipinf$tipinfo)
}
dis
boxplot(dis)
