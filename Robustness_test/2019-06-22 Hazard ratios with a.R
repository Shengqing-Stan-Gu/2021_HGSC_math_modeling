library(survival)

set.seed(100)
N = 50
r = abs(rnorm(N,mean=2,sd=0.8))
d = r/10
rc = abs(rnorm(N,mean=0.2,sd=0.08))
dc = abs(rnorm(N,mean=4.9,sd=1))
aseq = c(0.4,0.6,0.84,1.1,1.3,1.5)
a = matrix(nrow=N,ncol=length(aseq))
for (i in 1:length(aseq)) {
  a[,i] = abs(rnorm(N,mean=aseq[i],sd=0.5*aseq[i]))
}
b = a/10
ac = a
bc = ac/5

epsilon = 10^(c(-8.4,-7.4,-6.4,-5.4,-4.4))
PDSM1 = 10^rnorm(N,mean=11.5,sd=0.4)
NACTM1 = PDSM1
M2 = 10^rnorm(N,mean=13,sd=0.4)
PDSM = PDSM1/(10^5.5)
NACTM = PDSM
u = 10^(c(-9.6,-8.6,-7.6,-6.6,-5.6))
uc = 10*u
deltat = 0.005
t = 4.2
gap = 1
gap1 = 0.3
randomM1 = 1e8
deltatao = log(1.1)/(a-b)
Cure = 120


setwd("~/Desktop/2014-5 Modeling Ovarian Progression/2018-08_Vary_Parameters_a/")
source("Transcendental_solution.R")
source("EstimateY.R")
source("PDS_treatment.R")
source("NACT_treatment.R")
source("Survivaltime.R")


SurvivalPDS = vector("numeric")
SurvivalNACT = vector("numeric")

for (l in 1:length(aseq)) {
  HRmatrix = matrix(nrow=length(epsilon),ncol=length(u))
  SurvDiff = matrix(nrow=length(epsilon),ncol=length(u))
  for (k in length(epsilon):1) {
    for (j in 1:length(u)) {
      Survivaltimes = read.csv(file=paste("./Survival curve/a",l," u",j," epsilon",k,".csv",sep=""))[,c(6,8)]
      SurvivalPDS = Survivaltimes[,1]
      SurvivalNACT = Survivaltimes[,2]
      Survival = c(SurvivalPDS,SurvivalNACT)
      Status = ifelse(Survival==125.2,0,1)
      Group = c(rep(0,N),rep(1,N))
      Combined = data.frame(Survival,Status,as.factor(Group))
      HRatio = coxph(Surv(Survival,Status) ~ Group, data=Combined)
      HRmatrix[length(epsilon)+1-k,j] = HRatio$coefficients
      SurvDiff[length(epsilon)+1-k,j] = median(SurvivalPDS)-median(SurvivalNACT)
    }
  }
  write.csv(HRmatrix,file=paste("./Hazard ratio/HR with a",l,".csv",sep=""))
  write.csv(SurvDiff,file=paste("./Survival difference/SurvDiff with a",l,".csv",sep=""))
}


CurvePDS = vector("numeric")
CurveNACT = vector("numeric")
Timefollow = seq(from=0,to=84,length.out=200)

for (l in 1:length(aseq)) {
  pdf(paste("./Survival curve/Survival with a",l,".pdf",sep=""),width=12,height=10)
  par(mfrow=c(length(epsilon),length(u)),mar=c(2,2,0.5,0.5))
  for (k in length(epsilon):1) {
    for (j in 1:length(u)) {
      Survivaltimes = read.csv(file=paste("./Survival curve/a",l," u",j," epsilon",k,".csv",sep=""))[,c(6,8)]
      SurvivalPDS = Survivaltimes[,1]
      SurvivalNACT = Survivaltimes[,2]
      for (i in 1:length(Timefollow)) {
        CurvePDS[i] = (length(SurvivalPDS)+1-rank(c(Timefollow[i],SurvivalPDS))[1])/length(SurvivalPDS)
        CurveNACT[i] = (length(SurvivalNACT)+1-rank(c(Timefollow[i],SurvivalNACT))[1])/length(SurvivalNACT)
      }
      plot(Timefollow,CurvePDS,type="l",xlim=c(0,85),ylim=c(0,1),xaxt='n',lwd=2)
      points(Timefollow,CurveNACT,type="l",lwd=2,col=2)
      axis(side=1,at=c(0,12,24,36,48,60,72,84))
    }
  }
  dev.off()
}



