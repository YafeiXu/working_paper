##############################################################################################################################################################################################################################################################################################################################
##############################################################################################################################################################################################################################################################################################################################
############### multiple vine cpm detection ---->>> 20170116 测试 20170121 测试 ###############################################################################################################################################################################################################################################################################################################
##############################################################################################################################################################################################################################################################################################################################
##############################################################################################################################################################################################################################################################################################################################
# 20170125 增添防止2020问题，即剩余数据集不足以做hn, i.e. warmup is less than 20
# 20170201 adding new quarantine with warnNr=110
# vine多点识别机unique change point detection
maxEneryCPM=function(indx){
tryCatch({
library(MASS)
library(NPMVCP)
simNr=300
Sigma2 <- .2*diag(5, 10, 10)
Mean2=rep(0,10)
sim2=(mvrnorm(n = simNr, Mean2, Sigma2))

Sigma3 <- diag(5, 10, 10)
Mean3=rep(0,10)
sim3=(mvrnorm(n = simNr, Mean3, Sigma3))

# mu <- rep(0,5); sigma <-  matrix(c(1,0,0,0,0,  0,1,0,0,0,   0,0,1,0,0,  0,0,0,1,0,  0,0,0,0,1),5,5); nu <- 5
# n <- simNr # Number of draws
# sim5 <- t(t(mvrnorm(n, rep(0, length(mu)), sigma) / sqrt(nu / rchisq(n, nu))) + mu)

# mu <- rep(0,5); sigma <-  matrix(c(1,0,0,0,0,  0,1,0,0,0,   0,0,1,0,0,  0,0,0,1,0,  0,0,0,0,1),5,5); nu <- 5
# n <- simNr # Number of draws
# sim4 <- t(t(mvrnorm(n, rep(0, length(mu)), sigma) / sqrt(nu / rchisq(n, nu))) + mu)

sim6=rbind(sim2[1:32,],(sim3)[1:20,])

# ptm <- proc.time()
# maxEneryCPM(data1=sim6,wNr=30,permNr=200,alpha=1/200)
# ((proc.time() - ptm)/60)[3]
data1=sim6
wNr=32
permNr=200
alpha=1/200

 library(copula)
library(VineCopula)
library(NPMVCP)
############################################################################################################################################################################################# vine单点识别机unique change point detection
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################

uniqueVineCPMenergy1=function(dataInput,confLev,simNr,warmNr,stepLength,runIndex){
#incomeDataUpdate is 20pts+2,i.e 22 points, in which 20pt are for training(simulation or warmup points)
#stepLength从t1到往后需要模拟的步数
#dataInput=read.csv("C:/Rcode-P2/sp500AfterPQRres200-tau03-Mood-20170105-2.csv")[,-1][1:30,];simNr=10;stepLength=11;confLev=.99;library(copula);dataInput=pobs(dataInput);warmNr=20;confLev=.99
#dataInput
#dataInput=read.csv("C:/Rcode-P2/sp500AfterPQRres200-tau03-Mood-20170105-2.csv")[,-1];simNr=100;stepLength=40;confLev=.99;dataInput=pobs(dataInput);warmNr=20;confLev=.99
##dataInput=setComb[1:40,];simNr=10;stepLength=50;confLev=.99;warmNr=20;confLev=.99;N=2
library(VineCopula)
library(copula)


# rnorm.cop01 =matrix(rnorm(450, mean = 100, sd = 1),ncol=3)
# rnorm.cop09 =matrix(runif(450, min = 0, max = 1),ncol=3)
# rnorm.cop02 =matrix(rnorm(450, mean = 10, sd = 1),ncol=3)


# setComb90=rbind(rnorm.cop01,rnorm.cop09,rnorm.cop02)
#setComb90=rbind(rnorm.cop01,rnorm.cop09,rnorm.cop02)
# setComb90=rbind(rnorm.cop01,rnorm.cop09,rnorm.cop02)
# dataInput=setComb90[1:300,];simNr=10;stepLength=50;confLev=.99;warmNr=110;confLev=.99;N=2
# runIndex=1





##################################################################################################### max energy test function
##################################################################################################### max energy test function
##################################################################################################### max energy test function
maxEnergyVector=function(data2,permSampleSize){# give the cuttoffs for max energy test statistic

        maxeg1=function(data1){# return only one max energy statistic, i.e. a number
        #data1=permData[[1]]
        library(energy)
        dataLength=dim(data1)[1]
            if(dataLength<wNr){
                #print("ptNotEnough")
                return(999999)
            }else{
                runs=dataLength-wNr+1
                Dkn=numeric()
                for(i in 1:runs){
                #i=401
                    Dkn[i]=eqdist.e(data1, c((wNr/2+i-1),(dataLength-(wNr/2+i-1))))         
                }
                return(max(Dkn))
            }

        }
# rnorm.cop01 =matrix(rnorm(450, mean = 0, sd = 1),ncol=3)
# rnorm.cop09 =matrix(runif(450, min = 0, max = 1),ncol=3)
# rnorm.cop02 =matrix(rnorm(450, mean = 10, sd = 1),ncol=3)


# setComb90=rbind(rnorm.cop01,rnorm.cop09,rnorm.cop02)
#setComb90=rbind(rnorm.cop01,rnorm.cop09,rnorm.cop02)
#data2=setComb90=rbind(rnorm.cop01,rnorm.cop09,rnorm.cop02);permSampleSize=10
        #data2=setComb90;permSampleSize=10

        dataLength=dim(data2)[1]
        permMatrix=list()
        #dataLength=10;permSampleSize=5
        permMatrix=t(as.matrix(sample(1:dataLength, dataLength, replace = FALSE, prob = NULL)))
        for(j in 1:permSampleSize){
        #permMatrix=rbind(permMatrix, sample(1:dataLength, dataLength, replace = FALSE, prob = NULL))
        sameVec=TRUE
            while(sameVec==TRUE){
            newSample=sample(1:dataLength, dataLength, replace = FALSE, prob = NULL)
            compareVector=apply(permMatrix, 1, function(x, want) isTRUE(all.equal(x, want)), newSample)
                if(length(which(compareVector==TRUE))>=1){
                sameVec=TRUE
                }else{
                permMatrix=rbind(newSample,permMatrix )# permMatrix is permSampleSize*dataLength
                sameVec=FALSE
                }
            }

        }
        #library(gtools)
        #permMatrix=permutations(n = dataLength,r=dataLength)
        #permSampleIndex=sample(1:dim(permMatrix)[1], permSampleSize, replace = FALSE, prob = NULL)
        permData=list()
        for(i in 1:permSampleSize){
        permData[[i]]=data2[permMatrix[i,],]


        }
        #maxeg1(permData[[1]])
                library(parallel)#做并行计算obtain a vecor of max energy
                cl <- makeCluster(detectCores())  
                one=parSapply(cl, X=permData, FUN=maxeg1)# 并行计算函数 return a vector of max energy test stat.               
                stopCluster(cl)
                return(one)

}
##################################################################################################### max energy test function
##################################################################################################### max energy test function
##################################################################################################### max energy test function

#############################################################
dataLength=dim(as.matrix(dataInput))[1]
d=dataInput; detectionTime=numeric(); estimatedCP=numeric()
#dataSetUpdate=d
confLev=confLev; simNr=simNr; warmNr=warmNr; dIndex=0
#??????????????????????????????????????? how to deal with this N
#!!!!!!!!!!!!!!!!!! this N is important !!!!!!!!!!!!!!!!
incomeDataLength=warmNr+1#-----------------------------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!此处需要进入循环来更新incomeDataUpdate此处N1是一个循环变量类似于i
#stepLengthUpdate=stepLength# 此处有疑问函数底部缺少该变量的更新变量，解决
source("C:/Rcode-P2/hn20170201-hnqEnergy-test.r") # 导入目标函数
        # tryCatch({     
if(dim(d)[1]>=(wNr+1)){#qrt
#cutoff=hnqEnergy1(d[1:110,], simNr, stepLengthUpdate, confLev=.99,warmNr)[[2]]#求出模拟出来的control limit????????????????????????????可能不是d
#cutoff=maxEnergyVector(d[1:110,],1000)
    incomeDataUpdate=dataInput[1:incomeDataLength,]# take the first 22 obs as incomeDataUpdate??????????????????????????????????????????????????????????????
    #incomeDataUpdateTrain=incomeDataUpdate[1:110,]
#plot(cutoff,type="b")
    

    #### 逐点检验循环
    noChange=TRUE
maxDkn=0
    while(noChange==TRUE & incomeDataLength<=dataLength){
       
#noChange==TRUE & incomeDataLength<=dataLength
print(paste("incomeStepLength=",incomeDataLength))
plot(1:incomeDataLength, main=paste("RUN=",runIndex," incomeStepLength=",incomeDataLength,"  finished=",100*round(incomeDataLength/dataLength, digit=2),"%"))
        realSectionData=list()

        ## 分割数据集

    library(VineCopula)

        
        
         # for(i in 1:(dim(incomeDataUpdate)[1]-20-1)){
		
			 # realSectionData[[i]]=as.matrix(incomeDataUpdate[1:(20+i+1),])
			
         # }
		 
		 for(i in 1:(dim(incomeDataUpdate)[1]-wNr+1)){
		
			 realSectionData[[i]]=as.matrix(cbind(incomeDataUpdate,i))
			
         }
        ##        

        
        # 求检验统计量Dkn
        vDn=function(Mat){#20170115找不到对象'incomeDataUpdate' 解决
         tryCatch({                     #要解决此问题关键是输入变量section中必须包含incomeDataUpdate的数据
#Mat=realSectionData[[1]]
            #library(VineCopula) 
			#if(dim(Mat)[1]>170){Mat=Mat[(dim(Mat)[1]-170+1):dim(Mat)[1],]}# quarantine 只保留最后50个点 20170117

            #Mat=as.data.frame(section)
			sectionPt=Mat[1,dim(Mat)[2]]
			Mat=Mat[,-dim(Mat)[2]]
            #trainingSet=Mat[1:(sectionPt+49),]#read.csv("C:/Rcode-P2/sp500AfterPQRres200-tau03-20161222-1.csv")[,-1];trainingSet=pobs(trainingSet)################################此处需要修改？？？已解决20170115 配合618行
            
            #if(dim(Mat)[1]>70){Mat=Mat[(dim(Mat)[1]-70+1):dim(Mat)[1],]}# quarantine 只保留最后50个点 20170117
             # lengthTrainingSet=dim(trainingSet)[1]
             # lengthSimSet=dim(Mat)[1]-lengthTrainingSet
            # tempDn=numeric()
            #vineStructureTemp=RVineStructureSelect(trainingSet, c(3:5))# 所有test的模型都是基于这个vine structure建模的数据来自于training data 
            #for(i in 1:(lengthSimSet-1)){
				# library(VineCopula)
				# library(copula)
                #vineStructureTemp=RVineStructureSelect(trainingSet, c(1:5))# 所有test的模型都是基于这个vine structure建模的数据来自于training data 
                #vineStructureTemp=RVineStructureSelect(trainingSet, c(3:5))
                # listPara=list()
                # listPara[[1]]=1#
                # listPara[[2]]=vineStructureTemp
                # listPara[[3]]="ECP2"
                # listPara[[4]]="CvM"
                # listPara[[5]]=2
                # library(VineCopula)#in referred function, all variables must be cleared
                #tempDn[i]=testStat=RVineGofTest(Mat[-(1:(20+i-1)),],listPara[[2]],listPara[[3]],listPara[[4]],listPara[[5]])$CvM
                #tryCatch({
                #tempDn=testStat=RVineGofTest(Mat[-(1:(sectionPt+19)),],listPara[[2]],listPara[[3]],listPara[[4]],listPara[[5]])$CvM
                
                library(energy)
                tempDn=eqdist.e(Mat, c((sectionPt+wNr/2-1),(dim(Mat)[1]-(sectionPt+wNr/2-1))))#qrt
                #}, error=function(e){return(222)})#111
            #}
            #max(tempDn)
			#
            return(tempDn)
            }, error=function(e){return(222)})#222
        }
        #
        

        
        ### 做并行计算求出Dkn
        library(parallel)#对统计量矩阵计算做并行计算
        cl <- makeCluster(detectCores())  
        oneVectorSimDnHERE=parSapply(cl, X=realSectionData, FUN=vDn)# 并行计算函数
        stopCluster(cl)
        oneVectorSimDnHERE1=oneVectorSimDnHERE
        if(length(which(oneVectorSimDnHERE1==222))>0){oneVectorSimDnHERE1[which(oneVectorSimDnHERE1==222)]=mean(oneVectorSimDnHERE1[-(which(oneVectorSimDnHERE1==222))])}#---------------------------------------- handle error value 20170123
		maxDkn=max(oneVectorSimDnHERE1); epochMaxDkn=which(oneVectorSimDnHERE1==max(oneVectorSimDnHERE1))+wNr-1
        ###
		#maxDkn=epochMaxDkn=numeric()
		# for(i in 1:(N-1)){
		#i=3
			# maxDkn[i]=oneVectorSimDn[[i]][1]
			# epochMaxDkn[i]=which(oneVectorSimDn[[i]][-1]==max(oneVectorSimDn[[i]][-1]))+20-1
		# }
        #maxDkn=oneVectorSimDn; epochMaxDkn=which(oneVectorSimDn==max(oneVectorSimDn))+20-1
        ###

        
        #### 条件判断是否maxDkn是否超过control limit--------------------------------------------------------------------------------- permutation
        cutoff=maxEnergyVector(incomeDataUpdate,permNr)
        maxDknPValue= 1-sum(cutoff<maxDkn)/(length(cutoff)+1)
        # cutoff=maxEnergyVector(incomeDataUpdate,500)
        # percentile <- ecdf(cutoff)
        # maxDknPValue=percentile(maxDkn)
        #if(maxDkn<=cutoff[incomeDataLength-110] ){#???????????????????????????????????存在BUG
        if(maxDknPValue>(alpha)){
            noChange=TRUE; incomeDataLength=incomeDataLength+1; if(incomeDataLength>dim(dataInput)){print("endReached")}else{incomeDataUpdate=dataInput[1:incomeDataLength,]}#???????????????????????????????????????????????????????
        }else{
            noChange=FALSE; incomeDataLength=incomeDataLength
        }

    }# 两种情况下noChange循环跳出1，出现间断点，2，到达最后一个记录
    ####
    if(noChange==FALSE & incomeDataLength<=dataLength){
    #if(noChange==FALSE){
        dIndex=dIndex+1
        detectionTime[dIndex]=incomeDataLength# estimate out Detection time point
        estimatedCP[dIndex]=epochMaxDkn# estimate out change time point
        returnList=list()
        returnList[[1]]=detectionTime
        returnList[[2]]= estimatedCP
        returnList[[3]]=oneVectorSimDnHERE
        returnList[[4]]=oneVectorSimDnHERE1
        returnList[[5]]=maxDkn
        returnList[[6]]=cutoff
        returnList
        #incomeDataUpdate=dataInput[-(1:epochMaxDkn),]
        return(returnList)
        #incomeDataLength=incomeDataLength+1
    }

     if(noChange==TRUE &  incomeDataLength>dataLength){
     
        incomeDataLength=incomeDataLength-1
        returnList=list()
        returnList[[1]]=2600000
        returnList[[2]]=2600000
        returnList[[3]]=oneVectorSimDnHERE
        returnList[[4]]=oneVectorSimDnHERE1
        returnList[[5]]=maxDkn
        returnList[[6]]=cutoff
        returnList
        return(returnList)
     }    
     if(noChange==FALSE & incomeDataLength>dataLength){
        dIndex=dIndex+1
        detectionTime[dIndex]=incomeDataLength# estimate out Detection time point
        estimatedCP[dIndex]=epochMaxDkn# estimate out change time point
        incomeDataLength=incomeDataLength-1
        returnList=list()
        returnList[[1]]=detectionTime
        returnList[[2]]= estimatedCP
        returnList[[3]]=oneVectorSimDnHERE
        returnList[[4]]=oneVectorSimDnHERE1
        returnList[[5]]=maxDkn
        returnList[[6]]=cutoff        
        returnList
        return(returnList)
     }
}else{return(2020)}
          #  }, error=function(e){return("errorInUniCPM")})#222

}
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
############################################################################################################################################################################################
# norm.cop01 <- normalCopula(rep(0.1,6), dim = 4, dispstr = "un")
# norm.cop09 <- normalCopula(rep(0.9,6), dim = 4, dispstr = "un")
# norm.cop02 <- normalCopula(rep(0.1,6), dim = 4, dispstr = "un")

# norm.cop01 <- claytonCopula(iTau(claytonCopula(), 0.9), dim = 4)
# norm.cop02 <- joeCopula(iTau(joeCopula(), 0.9), dim = 4)
# norm.cop09 <- normalCopula(rep(0.1,6), dim = 4, dispstr = "un")

# rnorm.cop01 =rCopula(25, norm.cop01)
# rnorm.cop09 =rCopula(25, norm.cop09)
# rnorm.cop02 =rCopula(25, norm.cop02)


#data(Alsmelterdata)
#setComb90=dd=Alsmelterdata
#setComb90=dd3norm1t250
#setComb90=read.csv("C:/Rcode-P2/ETFresidual.csv")[1:1000,2:6]#--------------------------------------------------------# try ETFresidual data set
# rnorm.cop01 =matrix(rnorm(360, mean = 100, sd = 10),ncol=3)
# rnorm.cop02 =matrix(runif(360, min = 0, max = 1),ncol=3)
# rnorm.cop09 =matrix(rnorm(360, mean = 1, sd = 1),ncol=3)


setComb90=data1#rbind(rnorm.cop01,rnorm.cop09,rnorm.cop02)
# setComb90=rbind(rnorm.cop01,rnorm.cop09,rnorm.cop02)
plot(setComb90[,1],type="b")

#source("C:/Rcode-P2/hn20170201-hnqEnergy-test.r") # 导入目标函数 
#source("C:/Rcode-P2/uniqueVineCPMenergy-20170204-test10-ok3.r") # 导入目标函数 
dataStart=setComb90
#trainingSet; simNr=100; stepLength=10;confLev=.99; warmNr=20
#trainingSet=read.csv("C:/Rcode-P2/sp500AfterPQRres200-tau03-Mood-20170105-2.csv")[,-1];trainingSet=pobs(trainingSet)
simNr=1000
stepLength=110
confLev=.99
warmNr=wNr#20
runIndex=1


##################################################################################################
#dataStart=read.csv("C:/Rcode-P2/sp500AfterPQRres200-tau03-Mood-20170105-2.csv")[,-1]  #----------------------------------------------------输入区间数据
#source("C:/Rcode-P2/uniqueVineCPMenergy-20170204-test10-ok3.r") # 导入目标函数
#dataStart=setComb120;simNr=100;stepLength=50;confLev=.99;warmNr=20;confLev=.99;runIndex=1
library(copula)
dataInput=dataStart
cp=dtt=numeric();mVCPMindex=1#;confLev=.99;simNr=30;warmNr=20;stepLength=41
#uVCPM=uniqueVineCPM(dataInput,confLev=.99,simNr,warmNr,stepLength,runIndex)
uVCPM=uniqueVineCPMenergy1(dataInput,confLev,simNr,warmNr,stepLength,runIndex)
oneVectorSimDnHERE1=list()
oneVectorSimDnHERE=list()
oneVectorSimDnHERE[[1]]=uVCPM[[3]]
oneVectorSimDnHERE1[[1]]=uVCPM[[4]]
cp[1]=uVCPM[[2]];dtt[1]=uVCPM[[1]]
cp;dtt
#runIndex=0
warn2020multi=0

while((cp[length(cp)])!=2600000 & dim(dataStart[-(1:cp[length(cp)]),])[1]>=(warmNr+1) & warn2020multi!=2020){
#is.na(cp[length(cp)])!=TRUE & dim(dataStart[-(1:cp[length(cp)]),])[1]>=(warmNr+2)
runIndex=runIndex+1
    #dataInput=pobs(dataStart[-(1:cp[length(cp)]),])
if((cp[length(cp)]+warmNr)>length(setComb90[,1])){
warn2020multi=2020
}else{
    dataInput=dataStart[-(1:cp[length(cp)]),]
             #k=tryCatch({    
    uResult=uniqueVineCPMenergy1(dataInput,confLev=.99,simNr,warmNr,stepLength=(dim(dataInput)[1]-wNr/2+5),runIndex)
                #}, error=function(e){return("errorInTheUniCPM")})#222
                #if(k!="errorInTheUniCPM"){
    #if(uResult!=2020 & uResult!="errorInUniCPM"){
                        if((uResult[[1]])!=2600000){
                            mVCPMindex=1+mVCPMindex
                            cp[mVCPMindex]=uResult[[2]]+cp[mVCPMindex-1]
                            dtt[mVCPMindex]=uResult[[1]]+cp[mVCPMindex-1]
                            oneVectorSimDnHERE[[mVCPMindex]]=uResult[[3]]
                            oneVectorSimDnHERE1[[mVCPMindex]]=uResult[[4]]
                            # plot(setComb90[,1],type="b")
                            # for(i in 1:length(cp)){
                            # if(length(cp)==1){lines(x=rep(cp,3),y=seq(min(setComb90[,1]),max(setComb90[,1]),length.out=3),col="red")}else{lines(x=rep(cp[i],3),y=seq(min(setComb90[,1]),max(setComb90[,1]),length.out=3),col="red")}                            
                            # }     
                        }else{warn2020multi=2020;print(paste("restSetLessThan22"))}   
                        # list4print=list()
                        # list4print[[1]]=cp
                        # list4print[[2]]=dtt
                        
                        #print(list4print)
                        #runIndex=runIndex+1
                        print(paste("changePt",cp,"and detectionPt",dtt))
                    #}else{warn2020multi=2020;print(paste("restSetLessThan22"))}
    # }else{if(uResult==2020){
            # warn2020multi=2020;print(paste("restSetLessThan22"));return(warn2020multi)#qrt
            # }else{print(paste("errorInUniCPM"))}
            # }

 }                          
}

                            plot(setComb90[,1],type="b")
                            for(i in 1:length(cp)){
                            if(length(cp)==1){lines(x=rep(cp,3),y=seq(min(setComb90[,1]),max(setComb90[,1]),length.out=3),col="red")}else{lines(x=rep(cp[i],3),y=seq(min(setComb90[,1]),max(setComb90[,1]),length.out=3),col="red")}                            
                            } 
                            
retList=list()
retList[[1]]=dtt
retList[[2]]=cp

return(c(dtt,cp))
#((proc.time() - ptm)/60)[3]

}, error=function(e){return(1999)})#111

}                            
                        
# BUG LOG
# 20170116-1 到今天的问题在于取20个点做training然后第22点必然是间断点，然后基本上每次相隔20个点就是间断点，这是个问题
# 即每次第22个点都是间断点，control limit的算法可能有问题
# 20170116-2
# Error in if (maxDkn <= cutoff[incomeDataLength - 20 - 1]) { : 这种问题是cutoff数量不足，解决办法是加大stepLength
#  需要TRUE/FALSE值的地方不可以用缺少值

# Error in dataInput[1:incomeDataLength, ] : subscript out of bounds
# Error in 1:cp[length(cp)] : argument of length 0

# >  uniqueVineCPM(dataInput,confLev,simNr,warmNr=20,stepLength)
# Error in dataInput[1:incomeDataLength, ] : subscript out of bounds

# 程序调试结果1
# >  dataInput=setComb;simNr=1000;stepLength=100;confLev=.99;warmNr=20;confLev=.99
# >  uniqueVineCPM(dataInput,confLev,simNr,warmNr=20,stepLength)
# [[1]]
# [1] 114

# [[2]]
# [1] 111

# There were 32 warnings (use warnings() to see them)
# > 
# > ((proc.time() - ptm)/60)[3]
 # elapsed 
# 2558.399 

# 程序调试结果2
# fLev=.99;N=2
# > 
# > ### test a set with 3 diverse distribution
# >  ptm <- proc.time()
# >  library(copula)
# > library(VineCopula)
# > 
# > clayton.cop <- claytonCopula(2, dim = 4)
# > rClayton=rCopula(50, clayton.cop)
# > frank.cop <- frankCopula(15, dim = 4)
# > rFrank=rCopula(50, frank.cop)
# > joe.cop <- joeCopula(3, dim = 4)
# > rJoe=rCopula(50, frank.cop)
# > setComb=rbind(rClayton,rFrank,rJoe)
# > 
# >  dataInput=setComb[1:70,];simNr=1000;stepLength=100;confLev=.99;warmNr=20;confLev=.99
# >  uniqueVineCPM(dataInput,confLev,simNr,warmNr=20,stepLength)
# Error in dataInput[1:incomeDataLength, ] : subscript out of bounds
# > 
# > ((proc.time() - ptm)/60)[3]
 # elapsed 
# 2513.611 
# > 
# > 

### 20170123
# last point should not be the cp. need to be corrected
# Error in if (maxDkn <= cutoff[incomeDataLength - 20 - 1]) { : 
  # 需要TRUE/FALSE值的地方不可以用缺少值

##############################################################################################################################################################################################################################################################################################################################
##############################################################################################################################################################################################################################################################################################################################
############### <<<---- multiple vine cpm detection 20170116 测试###############################################################################################################################################################################################################################################################################################################
##############################################################################################################################################################################################################################################################################################################################
##############################################################################################################################################################################################################################################################################################################################
