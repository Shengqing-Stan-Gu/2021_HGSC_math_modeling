setwd("~/Desktop/2014-5 Modeling Ovarian Progression/2015-5-28 Training parameters for the model/")

observePDS = c(0,5,8,8,0,4,15,60)
observeNACT = c(5,27,17,24,5,5,0,17)
Observe = c(observePDS,observeNACT)
u = 10^(c(-11,-10,-9,-8,-7,-6,-5))
epsilon = 10^(c(-10,-9,-8,-7,-6,-5,-4))

Fit = matrix(nrow=length(epsilon),ncol=length(u))

for (k in length(epsilon):1) {
  for (j in 1:length(u)) {
    predictPDS = vector("numeric")
    predictNACT = vector("numeric")
    predict = read.csv(file=paste("./Parameter training rough/u",j," epsilon",k,".csv",sep=""))[,c(6,8)]
    predictPDS[1] = sum(predict[,1]<12)
    predictPDS[2] = sum(predict[,1]<24&predict[,1]>=12)
    predictPDS[3] = sum(predict[,1]<36&predict[,1]>=24)
    predictPDS[4] = sum(predict[,1]<48&predict[,1]>=36)
    predictPDS[5] = sum(predict[,1]<60&predict[,1]>=48)
    predictPDS[6] = sum(predict[,1]<72&predict[,1]>=60)
    predictPDS[7] = sum(predict[,1]<84&predict[,1]>=72)
    predictPDS[8] = sum(predict[,1]>=84)
    predictNACT[1] = sum(predict[,2]<12)
    predictNACT[2] = sum(predict[,2]<24&predict[,2]>=12)
    predictNACT[3] = sum(predict[,2]<36&predict[,2]>=24)
    predictNACT[4] = sum(predict[,2]<48&predict[,2]>=36)
    predictNACT[5] = sum(predict[,2]<60&predict[,2]>=48)
    predictNACT[6] = sum(predict[,2]<72&predict[,2]>=60)
    predictNACT[7] = sum(predict[,2]<84&predict[,2]>=72)
    predictNACT[8] = sum(predict[,2]>=84)
    Predict = c(predictPDS,predictNACT)
    Fit[length(epsilon)-k+1,j] = sum((Predict-Observe)^2)^0.5
  }
}

write.csv(cbind(length(epsilon):1,Fit),file="./Parameter training rough/Deviation matrix.csv")
