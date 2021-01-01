## prepareJonesData_edition_003.R
## subtracts each s's own min RT to eliminated effect of T0
## this is the one to use
## changes the mapping of correct , incorrect for going into Noah fit code
## corrYN2 is 1 if right and 2 if wrong and -1 if missing

## prepareJonesData_edition_002.R
## changes coding for rep5new
## this is the one to use

## Analysis of speed accuracy data
## From E.J. Wagenmakers and Chris Donkin
## To me from Matthew Jones website
## downloaded August 9, 2020
## See Wagenmakers, EJ, Ratcliff, R, Gomez, P, and McKoon, G. (2008)
## A diffusion model account of criterion shifts in a lexical decision task
## Journal of Memory and Language
## Volume 58
## 140-159
#rm(list=ls())
#setwd("/Users/OUConline/Desktop/cppdissertationfiles")

setwd("C:/Users/Robert Gore/Desktop/dissertation/20200930")
dataset = read.table ("SpeedAccData.txt",header=FALSE,col.names=c("pnum","blocknum","practiceYN","speedaccuracyYN","stimulusID","freqCondition","Ch","RTsec","censorYN"))
attach(dataset)
library(tidyverse)

RHO = .50
ALPHA = numeric()
BETA = numeric()

ALPHA[1]=1
BETA[1]=1

for (i in 2:length(RTsec))
{
  if (pnum[i] == pnum[i-1])
  {
    if (blocknum[i] == blocknum[i-1])
    {
      if (freqCondition[i] < 4)
      {
        ALPHA[i] = RHO * (ALPHA[i-1]-1) + 1 + 1  
        BETA[i]  = RHO * (BETA[i-1]-1) + 1
      }
      else
      {
        ALPHA[i] = RHO * (ALPHA[i-1]-1) + 1   
        BETA[i]  = RHO * (BETA[i-1]-1) + 1 + 1
      }
    }
    else
    {
      ALPHA[i]=BETA[i]=1
    }
  }
  else
  {
    ALPHA[i]=BETA[i]=1
  }
}

MUz2word = ALPHA/(ALPHA+BETA)
MUz2nonword = BETA/(ALPHA+BETA)
MUZ2THIS = numeric()
for (i in 1:length(RTsec))
{
  if (freqCondition[i]<4)
  {
    MUZ2THIS[i] = MUz2word[i]
  }
  else # if (freqCondition[i]<9)
  {
    MUZ2THIS[i] = MUz2nonword[i]
  }
}

dataset$MUZ2THIS = MUZ2THIS

# blockPlot = function(PNUM, BLOCKNUM)
# {
#    	png(filename=paste0("blockplot_p",PNUM,"_b",BLOCKNUM,".png"))
#     plot(RTsec[pnum==PNUM & blocknum==BLOCKNUM],ylim=c(0,max(RTsec[pnum==PNUM & blocknum==BLOCKNUM])),type="l",main=paste0("Person: ",PNUM,", Block: ",BLOCKNUM,", Response Times, All Trials in Block"),ylab="Response Time (Sec)",xlab="Trial Number (Within Block)")
# 	iqr = quantile(RTsec[pnum==PNUM & blocknum==BLOCKNUM],.75)-quantile(RTsec[pnum==PNUM & blocknum==BLOCKNUM],.25)
# 	lowlim  = median(RTsec[pnum==PNUM & blocknum==BLOCKNUM]) - 1.5*iqr
# 	highlim = median(RTsec[pnum==PNUM & blocknum==BLOCKNUM]) + 1.5*iqr
# 	n = length(RTsec[pnum==PNUM & blocknum==BLOCKNUM])
# 	arrows(0,lowlim,n,lowlim,len=0)
# 	arrows(0,highlim,n,highlim,len=0)
# 	RTvector = RTsec[pnum==PNUM & blocknum == BLOCKNUM]
# 	for (i in 1:length(RTvector))
# 	{
# 		if (Ch[i]==-1)
# 		{
# 			arrows(i,0,i,max(RTvector),len=0,col="black")
# 		}
# 	}
# 	dev.off()
# }

# for (p in 1:max(pnum))
# {
#   for (b in 1:max(blocknum[pnum==p]))	
#    {
#      blockPlot(p,b)
#    }
# }

trialnum=numeric()
trialnum[1]=1
for (i in 2:length(RTsec))
{
	if (pnum[i] == pnum[i-1])
	{
		trialnum[i] = trialnum[i-1] + 1
	}
	else
	{
		trialnum[i] = 1
	}
}
dataset$trialnum=trialnum

trialNumWithinBlock=numeric()
trialNumWithinBlock[1]=1
for (i in 2:length(RTsec))
{
	if (pnum[i] == pnum[i-1] & blocknum[i] == blocknum[i-1])
	{
		trialnum[i] = trialnum[i-1] + 1
	}
	else
	{
		trialnum[i] = 1
	}
}
table(trialnum)
dataset$trialNumWithinBlock=trialnum

wordYN = rep(NA,length(RTsec))
wordYN[freqCondition < 4] = 1
wordYN[freqCondition > 3] = 0
dataset$wordYN=wordYN

stimHist5 = rep(NA,length(RTsec))
for (i in 5:length(RTsec))
{
	if (pnum[i] == pnum[i-4] & blocknum[i] == blocknum[i-4])
	{
		stimHist5[i] = 1e4*wordYN[i-4]+1e3*wordYN[i-3]+1e2*wordYN[i-2]+1e1*wordYN[i-1]+wordYN[i]
	}
}
dataset$stimHist5=stimHist5
#table(stimHist5)

corrYN = rep(NA, length(RTsec))
for (i in 1:length(RTsec))
{
    if (wordYN[i] == Ch[i] & Ch[i] > -1)
        {
 	        corrYN[i] = 1
        }
    else if (abs(wordYN[i]-Ch[i])==1 & Ch[i] > -1)
        {
        	corrYN[i] = 0
    }
  else if (Ch[i] == -1)
  {
    corrYN[i] = -1
  }
}
dataset$corrYN=corrYN
#table(corrYN)

#data2aggregate = as.data.frame(cbind(dataset$RTsec,dataset$MUZ2THIS,dataset$speedaccuracyYN,dataset$corrYN))
dataset$MUZ2THIS = round(dataset$MUZ2THIS,2)
data2aggregate = subset(dataset, censorYN != 1 & pnum != 2, select=c(RTsec, MUZ2THIS, speedaccuracyYN, corrYN,pnum))
aggregatedData = aggregate(data2aggregate, by=list(data2aggregate$MUZ2THIS,data2aggregate$speedaccuracyYN,data2aggregate$corrYN,data2aggregate$pnum),FUN="mean")

aggregData1 = subset(aggregatedData, Group.2==0 & Group.3==0, select=c(RTsec, Group.1,Group.4))
aggregData2 = subset(aggregatedData, Group.2==0 & Group.3==1, select=c(RTsec, Group.1,Group.4))
aggregData3 = subset(aggregatedData, Group.2==1 & Group.3==0, select=c(RTsec, Group.1,Group.4))
aggregData4 = subset(aggregatedData, Group.2==1 & Group.3==1, select=c(RTsec, Group.1,Group.4))

aggregData1$meanRTwrongAcc   = aggregData1$RTsec
aggregData2$meanRTcorrAcc    = aggregData2$RTsec
aggregData3$meanRTwrongSpeed = aggregData3$RTsec
aggregData4$meanRTcorrSpeed  = aggregData4$RTsec


mergedAggA = merge(aggregData1,aggregData2,by=c("Group.4","Group.1"))
mergedAggB = merge(aggregData3,aggregData4,by=c("Group.4","Group.1"))
mergedAggC = merge(mergedAggA, mergedAggB, by=c("Group.4","Group.1"))

mergedAggC$zeta = (mergedAggC$meanRTcorrSpeed - mergedAggC$meanRTwrongSpeed)*(mergedAggC$meanRTcorrAcc+mergedAggC$meanRTwrongAcc)/((mergedAggC$meanRTwrongAcc - mergedAggC$meanRTcorrAcc)*(mergedAggC$meanRTcorrSpeed+mergedAggC$meanRTwrongSpeed))
#plot(mergedAggC$zeta~mergedAggC$Group.1)
#arrows(0,0,1,0,len=0,lty=3)
#abline(lm(mergedAggC$zeta~mergedAggC$Group.1))

fourmeans = function(stimhist)
{
	mean1 = mean(RTsec[stimHist5==stimhist & speedaccuracyYN == 1 & corrYN == 1 & censorYN == 0],na.rm=TRUE)
	mean2 = mean(RTsec[stimHist5==stimhist & speedaccuracyYN == 1 & corrYN == 0 & censorYN == 0],na.rm=TRUE)
	mean3 = mean(RTsec[stimHist5==stimhist & speedaccuracyYN == 0 & corrYN == 1 & censorYN == 0],na.rm=TRUE)
	mean4 = mean(RTsec[stimHist5==stimhist & speedaccuracyYN == 0 & corrYN == 0 & censorYN == 0],na.rm=TRUE)
	return(list(speededCorrects = mean1, speededErrors = mean2, accuracyCorrects = mean3, accuracyErrors = mean4))
}

fourmeans(11111)

# rep 5 is number of past 5 stimuli that were repetitions - range 0 to 4

rep5 = rep(NA,length(RTsec))
rep5new = rep(NA,length(RTsec))

rep5new[stimHist5==00000]=4
rep5new[stimHist5==11111]=4

rep5new[stimHist5==00010]=3
rep5new[stimHist5==00100]=3
rep5new[stimHist5==01000]=3
rep5new[stimHist5==10000]=3
rep5new[stimHist5==01111]=3
rep5new[stimHist5==10111]=3
rep5new[stimHist5==11011]=3
rep5new[stimHist5==11101]=3

rep5new[stimHist5==00110]=2
rep5new[stimHist5==00111]=2
rep5new[stimHist5==01010]=2
rep5new[stimHist5==01011]=2
rep5new[stimHist5==01100]=2
rep5new[stimHist5==01101]=2
rep5new[stimHist5==10010]=2
rep5new[stimHist5==10011]=2
rep5new[stimHist5==10100]=2
rep5new[stimHist5==10101]=2
rep5new[stimHist5==11000]=2
rep5new[stimHist5==11001]=2

rep5new[stimHist5==00011]=1
rep5new[stimHist5==00101]=1
rep5new[stimHist5==01001]=1
rep5new[stimHist5==10001]=1
rep5new[stimHist5==01110]=1
rep5new[stimHist5==10110]=1
rep5new[stimHist5==11010]=1
rep5new[stimHist5==11100]=1

rep5new[stimHist5==00001]=0
rep5new[stimHist5==11110]=0

#plot(rep5new ~ MUZ2THIS)
#plot(rep5new[corrYN==1] ~ MUZ2THIS[corrYN==1])
dataset$rep5new=rep5new
# code below had a problem
#rep5[stimHist5==11111 | stimHist5==0] = 4
#rep5[stimHist5== 1111 | stimHist5==10000] = 3
#rep5[stimHist5==  111 | stimHist5==10111 | stimHist5== 1000 | stimHist5==11000] = 2
#rep5[stimHist5==   11 | stimHist5== 1011 | stimHist5==10011 | stimHist5==11011] = 1
#rep5[stimHist5==  100 | stimHist5== 1100 | stimHist5==10100 | stimHist5==11100] = 1
#rep5[stimHist5==    1 | stimHist5==  101 | stimHist5== 1001 | stimHist5== 1101] = 0
#rep5[stimHist5==10001 | stimHist5==10101 | stimHist5==11001 | stimHist5==11101] = 0
#rep5[stimHist5==   10 | stimHist5==  110 | stimHist5== 1010 | stimHist5== 1110] = 0
#rep5[stimHist5==10010 | stimHist5==10110 | stimHist5==11010 | stimHist5==11110] = 0
#dataset$rep5=rep5
dataset$rep5=rep5new

completeCaseVector = complete.cases(dataset$pnum, dataset$RTsec, dataset$corrYN, dataset$speedaccuracyYN, dataset$rep5)
#table(completeCaseVector)
dataset$completeCaseVector=completeCaseVector

#plotData1 <- data.frame (
#  pnum=dataset$pnum[dataset$completeCaseVector & dataset$censorYN != 1],
#  RTsec=dataset$RTsec[dataset$completeCaseVector & dataset$censorYN != 1], 
#  corrYN=dataset$corrYN[dataset$completeCaseVector & dataset$censorYN != 1],
#  speedaccuracyYN=dataset$speedaccuracyYN[dataset$completeCaseVector & dataset$censorYN != 1],
#  rep5=dataset$rep5[dataset$completeCaseVector & dataset$censorYN != 1]
#  )
#summary(dataset)
write.csv (dataset,file="speedAccData4Stata.csv")

dataset$RTsubjectMIN  = rep(NA,length(pnum))
dataset$RTsecMinusMin = rep(NA,length(pnum))
for (iii in 1:length(dataset$pnum))
{
  dataset$RTsubjectMIN[iii] = min(RTsec[pnum==pnum[iii]])
  dataset$RTsecMinusMin[iii] = RTsec[iii] - 0.9*min(RTsec[pnum==pnum[iii]])
}
write.csv (dataset,file="speedAccData4Stata.csv")

dataArray = rep(NA,(16*20*96*2))
dim(dataArray) = c(16,20,96,2)
dimnames (dataArray) = list(rep(NA,16),rep(NA,20),rep(NA,96),c("t_list","resp_list"))
dataArray

dataArrayExtras = rep(NA,(16*20*96*4))
dim(dataArrayExtras) = c(16,20,96,4)
dimnames(dataArrayExtras) = list(rep(NA,16),rep(NA,20),rep(NA,96),c("wordYN","speedaccuracyYN","freqCondition","censorYN"))
dataArrayExtras

dataset$corrYN2=rep(NA,length(dataset$corrYN))
for (i in 1:length(dataset$corrYN))
{
  if (dataset$corrYN[i]==1) {dataset$corrYN2[i]=1}
  else if (dataset$corrYN[i]==0) {dataset$corrYN2[i]=2}
  else if (dataset$corrYN[i]==-1) {dataset$corrYN2[i]=-1}
}
#table(dataset$corrYN,dataset$corrYN2,useNA="ifany")

for (i in 1:length(dataset$pnum))
{
  if (dataset$blocknum[i] > 0 & pnum[i]==1)
  {
    dataArray[dataset$pnum[i],dataset$blocknum[i],dataset$trialNumWithinBlock[i],1]=dataset$RTsecMinusMin[i]
    dataArray[dataset$pnum[i],dataset$blocknum[i],dataset$trialNumWithinBlock[i],2]=dataset$corrYN2[i]
    dataArrayExtras[dataset$pnum[i],dataset$blocknum[i],dataset$trialNumWithinBlock[i],1]=dataset$wordYN[i]
    dataArrayExtras[dataset$pnum[i],dataset$blocknum[i],dataset$trialNumWithinBlock[i],2]=dataset$speedaccuracyYN[i]
    dataArrayExtras[dataset$pnum[i],dataset$blocknum[i],dataset$trialNumWithinBlock[i],3]=dataset$freqCondition[i]
    dataArrayExtras[dataset$pnum[i],dataset$blocknum[i],dataset$trialNumWithinBlock[i],4]=dataset$censorYN[i]    
  }
  else if (dataset$blocknum[i] > 0 & pnum[i]>2)
  {
    dataArray[dataset$pnum[i]-1,dataset$blocknum[i],dataset$trialNumWithinBlock[i],1]=dataset$RTsecMinusMin[i]
    dataArray[dataset$pnum[i]-1,dataset$blocknum[i],dataset$trialNumWithinBlock[i],2]=dataset$corrYN2[i]
    dataArrayExtras[dataset$pnum[i]-1,dataset$blocknum[i],dataset$trialNumWithinBlock[i],1]=dataset$wordYN[i]
    dataArrayExtras[dataset$pnum[i]-1,dataset$blocknum[i],dataset$trialNumWithinBlock[i],2]=dataset$speedaccuracyYN[i]
    dataArrayExtras[dataset$pnum[i]-1,dataset$blocknum[i],dataset$trialNumWithinBlock[i],3]=dataset$freqCondition[i]
    dataArrayExtras[dataset$pnum[i]-1,dataset$blocknum[i],dataset$trialNumWithinBlock[i],4]=dataset$censorYN[i]    
  }
  
}


JonesData = array_tree (dataArray,margin=c(1,4,2))
#str(JonesData)

JonesDataExtras = array_tree(dataArrayExtras,margin=c(1,4,2))
#str(JonesDataExtras)

#detach(dataset)
# 
# attach(plotData1)
# aggData1 <- aggregate(plotData1, by=list(plotData1$corrYN,plotData1$speedaccuracyYN,plotData1$rep5), FUN=mean)  
# detach(plotData1)
# 
# attach(plotData1)
# aggData2 <- aggregate(plotData1, by=list(plotData1$corrYN,plotData1$speedaccuracyYN,plotData1$rep5,plotData1$pnum), FUN=mean)  
# detach(plotData1)
# 
# aggData1$speedaccuracyYN=as.factor(aggData1$speedaccuracyYN)
# levels(aggData1$speedaccuracyYN)[levels(aggData1$speedaccuracyYN)==1] = "Speed"
# levels(aggData1$speedaccuracyYN)[levels(aggData1$speedaccuracyYN)==0] = "Accuracy"
# 
# aggData2$speedaccuracyYN=as.factor(aggData2$speedaccuracyYN)
# levels(aggData2$speedaccuracyYN)[levels(aggData2$speedaccuracyYN)==1] = "Speed"
# levels(aggData2$speedaccuracyYN)[levels(aggData2$speedaccuracyYN)==0] = "Accuracy"
# 
# aggData1$rep5=as.factor(aggData1$rep5)
# levels(aggData1$rep5)[levels(aggData1$rep5)==0] = "0 Previous\nRepetitions" 
# levels(aggData1$rep5)[levels(aggData1$rep5)==1] = "1 Previous\nRepetition"
# levels(aggData1$rep5)[levels(aggData1$rep5)==2] = "2 Previous\nRepetitions"
# levels(aggData1$rep5)[levels(aggData1$rep5)==3] = "3 Previous\nRepetitions"
# levels(aggData1$rep5)[levels(aggData1$rep5)==4] = "4 Previous\nRepetitions"
# 
# aggData2$rep5=as.factor(aggData2$rep5)
# levels(aggData2$rep5)[levels(aggData2$rep5)==0] = "0 Previous\nRepetitions" 
# levels(aggData2$rep5)[levels(aggData2$rep5)==1] = "1 Previous\nRepetition"
# levels(aggData2$rep5)[levels(aggData2$rep5)==2] = "2 Previous\nRepetitions"
# levels(aggData2$rep5)[levels(aggData2$rep5)==3] = "3 Previous\nRepetitions"
# levels(aggData2$rep5)[levels(aggData2$rep5)==4] = "4 Previous\nRepetitions"
# 
# aggData2_1=aggData2[aggData2$pnum==1,]
# aggData2_2=aggData2[aggData2$pnum==2,]
# aggData2_3=aggData2[aggData2$pnum==3,]
# aggData2_4=aggData2[aggData2$pnum==4,]
# aggData2_5=aggData2[aggData2$pnum==5,]
# aggData2_6=aggData2[aggData2$pnum==6,]
# aggData2_7=aggData2[aggData2$pnum==7,]
# aggData2_8=aggData2[aggData2$pnum==8,]
# aggData2_9=aggData2[aggData2$pnum==9,]
# aggData2_10=aggData2[aggData2$pnum==10,]
# aggData2_11=aggData2[aggData2$pnum==11,]
# aggData2_12=aggData2[aggData2$pnum==12,]
# aggData2_13=aggData2[aggData2$pnum==13,]
# aggData2_14=aggData2[aggData2$pnum==14,]
# aggData2_15=aggData2[aggData2$pnum==15,]
# aggData2_16=aggData2[aggData2$pnum==16,]
# aggData2_17=aggData2[aggData2$pnum==17,]
# 
# #ggplot(aggData1, aes(x=as.factor(corrYN), y=RTsec, fill=as.factor(corrYN)))+geom_bar(stat="identity")+facet_grid(speedaccuracyYN~rep5)+scale_fill_discrete(name="Correct\nResponse")
# ggplot(aggData1, aes(x=as.factor(corrYN), y=RTsec, fill=as.factor(corrYN)))+geom_bar(stat="identity")+facet_grid(speedaccuracyYN~rep5)+scale_fill_discrete(name="Response",breaks=c(0,1),labels=c("Wrong","Correct"))+xlab("Response Correctness")+ylab("Response Time (S)")+ggtitle("Overall Means Across Participants")+theme_minimal() + ggsave("RR_fig1.png", width=6, height=4)
# 
# ggplot(aggData2_1, aes(x=as.factor(corrYN), y=RTsec, fill=as.factor(corrYN)))+geom_bar(stat="identity")+facet_grid(speedaccuracyYN~rep5)+scale_fill_discrete(name="Response",breaks=c(0,1),labels=c("Wrong","Correct"))+xlab("Response Correctness")+ylab("Response Time (S)")+ggtitle("Participant 1")+theme_minimal() + ggsave("RR_fig2_1.png", width=6, height=4)
# 
# ggplot(aggData2_2, aes(x=as.factor(corrYN), y=RTsec, fill=as.factor(corrYN)))+geom_bar(stat="identity")+facet_grid(speedaccuracyYN~rep5)+scale_fill_discrete(name="Response",breaks=c(0,1),labels=c("Wrong","Correct"))+xlab("Response Correctness")+ylab("Response Time (S)")+ggtitle("Participant 2")+theme_minimal() + ggsave("RR_fig2_2.png", width=6, height=4)
# 
# ggplot(aggData2_3, aes(x=as.factor(corrYN), y=RTsec, fill=as.factor(corrYN)))+geom_bar(stat="identity")+facet_grid(speedaccuracyYN~rep5)+scale_fill_discrete(name="Response",breaks=c(0,1),labels=c("Wrong","Correct"))+xlab("Response Correctness")+ylab("Response Time (S)")+ggtitle("Participant 3")+theme_minimal() + ggsave("RR_fig2_3.png", width=6, height=4)
# 
# ggplot(aggData2_4, aes(x=as.factor(corrYN), y=RTsec, fill=as.factor(corrYN)))+geom_bar(stat="identity")+facet_grid(speedaccuracyYN~rep5)+scale_fill_discrete(name="Response",breaks=c(0,1),labels=c("Wrong","Correct"))+xlab("Response Correctness")+ylab("Response Time (S)")+ggtitle("Participant 4")+theme_minimal() + ggsave("RR_fig2_4.png", width=6, height=4)
# 
# ggplot(aggData2_5, aes(x=as.factor(corrYN), y=RTsec, fill=as.factor(corrYN)))+geom_bar(stat="identity")+facet_grid(speedaccuracyYN~rep5)+scale_fill_discrete(name="Response",breaks=c(0,1),labels=c("Wrong","Correct"))+xlab("Response Correctness")+ylab("Response Time (S)")+ggtitle("Participant 5")+theme_minimal() + ggsave("RR_fig2_5.png", width=6, height=4)
# 
# ggplot(aggData2_6, aes(x=as.factor(corrYN), y=RTsec, fill=as.factor(corrYN)))+geom_bar(stat="identity")+facet_grid(speedaccuracyYN~rep5)+scale_fill_discrete(name="Response",breaks=c(0,1),labels=c("Wrong","Correct"))+xlab("Response Correctness")+ylab("Response Time (S)")+ggtitle("Participant 6")+theme_minimal() + ggsave("RR_fig2_6.png", width=6, height=4)
# 
# ggplot(aggData2_7, aes(x=as.factor(corrYN), y=RTsec, fill=as.factor(corrYN)))+geom_bar(stat="identity")+facet_grid(speedaccuracyYN~rep5)+scale_fill_discrete(name="Response",breaks=c(0,1),labels=c("Wrong","Correct"))+xlab("Response Correctness")+ylab("Response Time (S)")+ggtitle("Participant 7")+theme_minimal() + ggsave("RR_fig2_7.png", width=6, height=4)
# 
# ggplot(aggData2_8, aes(x=as.factor(corrYN), y=RTsec, fill=as.factor(corrYN)))+geom_bar(stat="identity")+facet_grid(speedaccuracyYN~rep5)+scale_fill_discrete(name="Response",breaks=c(0,1),labels=c("Wrong","Correct"))+xlab("Response Correctness")+ylab("Response Time (S)")+ggtitle("Participant 8")+theme_minimal() + ggsave("RR_fig2_8.png", width=6, height=4)
# 
# ggplot(aggData2_9, aes(x=as.factor(corrYN), y=RTsec, fill=as.factor(corrYN)))+geom_bar(stat="identity")+facet_grid(speedaccuracyYN~rep5)+scale_fill_discrete(name="Response",breaks=c(0,1),labels=c("Wrong","Correct"))+xlab("Response Correctness")+ylab("Response Time (S)")+ggtitle("Participant 9")+theme_minimal() + ggsave("RR_fig2_9.png", width=6, height=4)
# 
# ggplot(aggData2_10, aes(x=as.factor(corrYN), y=RTsec, fill=as.factor(corrYN)))+geom_bar(stat="identity")+facet_grid(speedaccuracyYN~rep5)+scale_fill_discrete(name="Response",breaks=c(0,1),labels=c("Wrong","Correct"))+xlab("Response Correctness")+ylab("Response Time (S)")+ggtitle("Participant 10")+theme_minimal() + ggsave("RR_fig2_10.png", width=6, height=4)
# 
# ggplot(aggData2_11, aes(x=as.factor(corrYN), y=RTsec, fill=as.factor(corrYN)))+geom_bar(stat="identity")+facet_grid(speedaccuracyYN~rep5)+scale_fill_discrete(name="Response",breaks=c(0,1),labels=c("Wrong","Correct"))+xlab("Response Correctness")+ylab("Response Time (S)")+ggtitle("Participant 11")+theme_minimal() + ggsave("RR_fig2_11.png", width=6, height=4)
# 
# ggplot(aggData2_12, aes(x=as.factor(corrYN), y=RTsec, fill=as.factor(corrYN)))+geom_bar(stat="identity")+facet_grid(speedaccuracyYN~rep5)+scale_fill_discrete(name="Response",breaks=c(0,1),labels=c("Wrong","Correct"))+xlab("Response Correctness")+ylab("Response Time (S)")+ggtitle("Participant 12")+theme_minimal() + ggsave("RR_fig2_12.png", width=6, height=4)
# 
# ggplot(aggData2_13, aes(x=as.factor(corrYN), y=RTsec, fill=as.factor(corrYN)))+geom_bar(stat="identity")+facet_grid(speedaccuracyYN~rep5)+scale_fill_discrete(name="Response",breaks=c(0,1),labels=c("Wrong","Correct"))+xlab("Response Correctness")+ylab("Response Time (S)")+ggtitle("Participant 13")+theme_minimal() + ggsave("RR_fig2_13.png", width=6, height=4)
# 
# ggplot(aggData2_14, aes(x=as.factor(corrYN), y=RTsec, fill=as.factor(corrYN)))+geom_bar(stat="identity")+facet_grid(speedaccuracyYN~rep5)+scale_fill_discrete(name="Response",breaks=c(0,1),labels=c("Wrong","Correct"))+xlab("Response Correctness")+ylab("Response Time (S)")+ggtitle("Participant 14")+theme_minimal() + ggsave("RR_fig2_14.png", width=6, height=4)
# 
# ggplot(aggData2_15, aes(x=as.factor(corrYN), y=RTsec, fill=as.factor(corrYN)))+geom_bar(stat="identity")+facet_grid(speedaccuracyYN~rep5)+scale_fill_discrete(name="Response",breaks=c(0,1),labels=c("Wrong","Correct"))+xlab("Response Correctness")+ylab("Response Time (S)")+ggtitle("Participant 15")+theme_minimal() + ggsave("RR_fig2_15.png", width=6, height=4)
# 
# ggplot(aggData2_16, aes(x=as.factor(corrYN), y=RTsec, fill=as.factor(corrYN)))+geom_bar(stat="identity")+facet_grid(speedaccuracyYN~rep5)+scale_fill_discrete(name="Response",breaks=c(0,1),labels=c("Wrong","Correct"))+xlab("Response Correctness")+ylab("Response Time (S)")+ggtitle("Participant 16")+theme_minimal() + ggsave("RR_fig2_16.png", width=6, height=4)
# 
# ggplot(aggData2_17, aes(x=as.factor(corrYN), y=RTsec, fill=as.factor(corrYN)))+geom_bar(stat="identity")+facet_grid(speedaccuracyYN~rep5)+scale_fill_discrete(name="Response",breaks=c(0,1),labels=c("Wrong","Correct"))+xlab("Response Correctness")+ylab("Response Time (S)")+ggtitle("Participant 17")+theme_minimal() + ggsave("RR_fig2_17.png", width=6, height=4)
# 
# ## getting moments and quantiles to use in developing data fit parameters for BAGL ##
# RTms = plotData1$RTsec[plotData1$pnum==11]*1000
# mean(RTms)
