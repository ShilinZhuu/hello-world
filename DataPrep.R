library(lubridate)
library(ggplot2)
library(reshape)
library(bsts)
library(stringr)
library(DescTools)
library(urca)f
library(tseries)

ExtractDatabyDateRange<-function(StartDate, EndDate, InputDF, Freq="Monthly")
# Purpose: extract the variables data by start and end date, by given frequency;
#          It has to be extracted from more granular frequency to less granular frequency, i.e., Monthly --> quarterly
#          It's better that column "Date" is at first column; "Date" is at character format; the date has to be consecutive;
  # StartDate & EndDate - have to be in character format at month end, quarter end, semi-annual end or year end, i.e., "3/31/2018"
  # InputDF - has to be a data frame with column "Date"; column "Date" is at character format
  # Freq -currently only support "Monthly","Quarterly","Semi-annually","Annually"
# the "Date" column in output data frame is already at Date object format
  {
  SD<-as.Date(StartDate,format="%m/%d/%Y")
  ED<-as.Date(EndDate,format="%m/%d/%Y")
  
  if (Freq=="Monthly") 
    {Intervals<-MonthDistance(ED,SD)+1} else 
      if (Freq=="Quarterly") 
      {Intervals<-(as.yearqtr(ED)-as.yearqtr(SD))*4+1} else if (Freq=="Semi-annually")
      {Intervals<-(year(ED)-year(SD))*2+1} else if (Freq=="Annually")
      {Intervals<-(year(ED)-year(SD))+1} else {}
  
  DateSeries<-matrix(1,Intervals,2)
  DateSeries<-as.data.frame(DateSeries)
  colnames(DateSeries)<-c("Date","dummy")
  
  DateSeries$Date<-as.Date(DateSeries$Date,format="%m/%d/%Y")
  DateSeries$Date[1]<-SD
  
  for (l in 2:nrow(DateSeries))
  {DateSeries$Date[l]<-LastDayOfMonth(
                                      AddMonths(DateSeries$Date[l-1],
                                                ifelse(Freq=="Monthly",1,
                                                       ifelse(Freq=="Quarterly",3,
                                                              ifelse(Freq=="Semi-annually",6,
                                                                     ifelse(Freq=="Annually",12,0)
                                                                     )
                                                              )
                                                       )
                                                )
                                      )
  }
  
  
  InputDF$Date <- as.Date(InputDF$Date,format="%m/%d/%Y")
  VarData<-merge(DateSeries,InputDF,by="Date",all.x=TRUE)
  VarData<-VarData[,-2]
  
  return(VarData)}

TransformVarData<-function(VarData,MaxLag=2,OrigDataOut=FALSE)
# Purpose: transform  the data into yoy, qoq, yd,qd, log, log_yd, log_qd, and then lag it with 1 to MaxLag
  # VaRData - expected to have "Date" column, the "Date" column needs to be in Date object format
  # MaxLag - the maximum lag number to append to the data set
  # OrigDataOut - whether to output the original data 
{

  SignMatrix<-matrix(0,1,ncol(VarData))
  for (i in 1:ncol(VarData))
  {
    if (sum(na.omit(VarData[,i]<=0))>0)
      SignMatrix[i]=1 
  }
  
  NoNegCol<-colnames(VarData) %in% colnames(VarData)[which(SignMatrix==1)]
  n<-nrow(VarData)
  
  if (MonthDistance(VarData$Date[2],VarData$Date[1])==12) {IntM=12} else if 
  (MonthDistance(VarData$Date[2],VarData$Date[1])==6) {IntM=6} else if 
  (MonthDistance(VarData$Date[2],VarData$Date[1])==3) {IntM=3} else {IntM=1}
  
  VarData_noDate<-VarData[,!(colnames(VarData) %in% c("Date"))]
  
  ##yoy
  yoy <- (VarData_noDate[(12/IntM+1):n,] / VarData_noDate[1:(n-12/IntM),])-1
  yoy_temp <- matrix(NA, 12/IntM, ncol(yoy))
  colnames(yoy_temp) <- colnames(yoy)
  yoy <- rbind.data.frame(yoy_temp, yoy)
  colnames(yoy) <- paste0(colnames(VarData_noDate), '_yoy', sep='')
  
  ##qoq
  qoq <- (VarData_noDate[(3/IntM+1):n,] /VarData_noDate[1:(n-3/IntM),])-1
  qoq_temp <- matrix(NA, 3/IntM, ncol(qoq))
  colnames(qoq_temp) <- colnames(qoq)
  qoq <- rbind.data.frame(qoq_temp, qoq)
  colnames(qoq) <- paste0(colnames(VarData_noDate), '_qoq', sep='')
  
  ##log
  logsub <- VarData_noDate[,!ifelse(sum(NoNegCol)==0,0,NoNegCol)]
  log_temp <- log(logsub)
  colnames(log_temp) <- paste0(colnames(logsub), '_log', sep='')
  
  #combine original and log data
  VarData_noDate<-cbind.data.frame(VarData_noDate,log_temp)
  
  ##yd
  yd <- VarData_noDate[(12/IntM+1):n,] - VarData_noDate[1:(n-12/IntM),]
  yd_temp <- matrix(NA, 12/IntM, ncol(yd))
  colnames(yd_temp) <- colnames(yd)
  yd <- rbind.data.frame(yd_temp, yd)
  colnames(yd) <- paste0(colnames(VarData_noDate), '_yd', sep='')
  
  ##qd
  qd <- VarData_noDate[(3/IntM+1):n,] - VarData_noDate[1:(n-3/IntM),]
  qd_temp <- matrix(NA, 3/IntM, ncol(qd))
  colnames(qd_temp) <- colnames(qd)
  qd <- rbind.data.frame(qd_temp, qd)
  colnames(qd) <- paste0(colnames(VarData_noDate), '_qd', sep='')
  
  VarData_noDate <- cbind.data.frame(VarData_noDate, yd, qd, yoy, qoq)
  
  ## Add Lags
  
  VarData_noDate_orig<-VarData_noDate
  
  for (i in 1:MaxLag)
    
  {l_i <- VarData_noDate_orig[1:(n-i),]
  t_i <- matrix(NA, i, dim(VarData_noDate_orig)[2])
  colnames(t_i) <- colnames(l_i)
  l_i <- rbind.data.frame(t_i, l_i)
  colnames(l_i) <- paste0(colnames(l_i), '_lag',i, sep='')
  
  if (OrigDataOut==TRUE) 
    {VarData_noDate <- cbind.data.frame(VarData_noDate, l_i)} else
    { if (i==1) {VarData_noDate<-l_i} else
        {VarData_noDate <- cbind.data.frame(VarData_noDate, l_i)}
      }
  }
  
  VarData_output<-cbind(VarData$Date,VarData_noDate)
  return(VarData_output)
}

CoIntTest<-function(TargData,PredData,SignLvl=1)
  # TargData and PredData should be without Date
  # SignLvl (1,2,3) --> (10%,5%,1%)
{
  VarData<-cbind(TargData,PredData)
  #combn_var_n<-combn(c(1:(ncol(VarData))),2)
  CoInt<-cor(as.matrix(VarData))
  CoInt<-CoInt[,1:ncol(TargData)]
  for (j in 1:ncol(TargData))
  {
    for (i in 1:nrow(VarData))
      { if (inherits(try(ca.jo(na.omit(data.frame(VarData[,i],VarData[,j])),type="trace", K=2, ecdet="none", spec="longrun")),"try-error")) 
        {coint_ts_i<-10000} else 
        {coint_i<-ca.jo(na.omit(data.frame(VarData[,i],VarData[,j])),type="trace", K=2, ecdet="none", spec="longrun")
        coint_ts_i<-coint_i@teststat[2]-coint_i@cval[2,SignLvl]
        }
      CoInt[i,j]<-coint_ts_i
    
      }
  }
  return(CoInt)  
}

CorrelTest<-function(TargData,PredData,corr_use="pairwise.complete.obs")
  # TargData and PredData should be without Date
  # corr_use will feed into the use parameters for cor function: "everything", "all.obs","complete.obs","na.or.complete","pairwise.complete.obs"
{
  VarData<-cbind(TargData,PredData)
  #combn_var_n<-combn(c(1:(ncol(VarData))),2)
  Correl<-cor(as.matrix(VarData),use=corr_use)
  Correl<-Correl[,1:ncol(TargData)]

  return(Correl)  
}

StationTest<-function(VarData)
  # ADF stationarity test for variables, return the p-values, if <0.1 or 0.05, it's stationary.
{
  StationarityTest<-matrix(1,1,ncol(VarData))
  for (i in 2:ncol(StationarityTest))
  {
    x_temp <- VarData[,i]
    x_temp<-na.omit(as.numeric(x_temp))
    if (sum(!is.finite(x_temp))>0 )
        {x_temp<-x_temp[which(is.finite(x_temp)==TRUE)]
    }
    if (inherits(try(adf.test(x_temp)),"try-error")) {StationarityTest[1,i]<-100}
    else {StationarityTest[1,i] <-adf.test(x_temp)$p.value}
    #adf.test(as.numeric(x_temp))$p.value
    #adfTest(na.omit(as.numeric(x_temp)), lags=0,type="ct")@test$p.value
    
  }
  colnames(StationarityTest)<-colnames(VarData)
  return(StationarityTest)
  }

setwd("C:\\Works\\CMBS EW")

Data1 <- read.csv('BBG_Variables.csv')
Data2 <- read.csv('Other_Variables.csv')
Data3 <- read.csv('Target_Variable.csv')

VarData1<-ExtractDatabyDateRange('6/30/2000','9/30/2018',Data1,"Quarterly")
VarData2<-ExtractDatabyDateRange('6/30/2000','9/30/2018',Data2,"Quarterly")
TargetVarData<-ExtractDatabyDateRange('9/30/2000','9/30/2018',Data3,"Quarterly")

VarData<-merge(VarData1,VarData2,by="Date",all.x=TRUE)
#VarData_q<-ExtractDatabyDateRange('9/30/2000','9/30/2018',VarData,"Quarterly")

VarData_Trans<-TransformVarData(VarData,2,FALSE)
colnames(VarData_Trans)[1]<-"Date"

TargetVarData_new<-cbind(TargetVarData,TargetVarData)
TargetVarData_new<-TargetVarData_new[,-3]
colnames(TargetVarData_new)[3]<-"Dummy"
TargetVarData_Trans<-TransformVarData(TargetVarData_new,2,TRUE)
colnames(TargetVarData_Trans)[1]<-"Date"
TargetVarData_Trans<-TargetVarData_Trans[,-which(substring(colnames(TargetVarData_Trans),1,5)=="Dummy")]
TargetVarData_Trans<-TargetVarData_Trans[,-which(str_detect(colnames(TargetVarData_Trans), "_lag")==TRUE)]

write.csv(VarData_Trans,"VarData_Trans.csv",row.names = FALSE)
write.csv(TargetVarData_Trans,"TargetVarData_Trans.csv",row.names = FALSE)

CoInt<-CoIntTest(TargetVarData_Trans[,-1],VarData_Trans[-1,-1],SignLvl=1)
StaTest<-StationTest(cbind(TargetVarData_Trans[,-1],VarData_Trans[-1,-1]))

Correl<-CorrelTest(TargetVarData_Trans[,-1],VarData_Trans[-1,-1],corr_use="pairwise.complete.obs")
Correl<-Correl[-(1:8),]
Correl_order<-Correl[order(Correl[,1]),]
Correl_order<-na.omit(Correl_order)
View(Correl_order)
write.csv(Correl_order,"Correl_order.csv",row.names = TRUE)
