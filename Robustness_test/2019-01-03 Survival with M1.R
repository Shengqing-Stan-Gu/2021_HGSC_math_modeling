set.seed(100)
N = 50
r = abs(rnorm(N,mean=2,sd=0.8))
d = r/10
rc = abs(rnorm(N,mean=0.2,sd=0.08))
dc = abs(rnorm(N,mean=4.9,sd=1))
a = abs(rnorm(N,mean=0.84,sd=0.42))
b = a/10
ac = a
bc = ac/5

epsilon = 10^(c(-8.4,-7.4,-6.4,-5.4,-4.4))
M1seq = c(10.5,11,11.5,12,12.4,12.8)
PDSM1 = matrix(nrow=N,ncol=length(M1seq))
for (i in 1:length(M1seq)) {
  PDSM1[,i] = 10^(rnorm(N,mean=M1seq[i],sd=0.4))
}
NACTM1 = PDSM1
M2 = 10^rnorm(N,mean=13,sd=0.4)
PDSM = PDSM1[,3]/(10^5.5)
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

setwd("~/Desktop/2014-5 Modeling Ovarian Progression/2019-01_Vary_Parameters_M1/")
source("Transcendental_solution.R")
source("EstimateY.R")
source("PDS_treatment.R")
source("NACT_treatment.R")
source("Survivaltime.R")

postPDSx = vector("numeric")
postPDSy = vector("numeric")
postNACTx = vector("numeric")
postNACTy = vector("numeric")
SurvivalPDS = vector("numeric")
ChemocountPDS = vector("numeric")
SurvivalNACT = vector("numeric")
ChemocountNACT = vector("numeric")

for (l in 1:length(M1seq)) {
  y0PDS = read.csv(file=paste("./Expected y at diagnosis/yexpected with M1_",l,".csv",sep=""))[,-1]
  y0NACT = y0PDS
  for (k in 1:length(epsilon)) {
    for (j in 1:length(u)) {
      for (i in 1:N) {
        postPDS = PDStreat(r=r[i],d=d[i],a=a[i],b=b[i],rc=rc[i],dc=dc[i],ac=ac[i],bc=bc[i],epsilon=epsilon[k],M=PDSM[i],M1=PDSM1[i,l],u=u[j],uc=uc[j],y0=y0PDS[i,j])
        postPDSx[i] = postPDS[1]
        postPDSy[i] = postPDS[2]
        if (y0NACT[i,j]*exp((ac[i]-bc[i])*t/2+(a[i]-b[i])*gap1)<M2[i]) {
          postNACT = NACTtreat(r=r[i],d=d[i],a=a[i],b=b[i],rc=rc[i],dc=dc[i],ac=ac[i],bc=bc[i],epsilon=epsilon[k],M=NACTM[i],M1=NACTM1[i,l],u=u[j],uc=uc[j],y0=y0NACT[i,j])
        } else {
          postNACT = c(0,y0NACT[i,j]*exp((ac[i]-bc[i])*t/2+(a[i]-b[i])*gap1))
        }
        postNACTx[i] = postNACT[1]
        postNACTy[i] = postNACT[2]
        if (postPDSx[i]+postPDSy[i]<M2[i]) {
          PDSresult = Survive(postPDSx[i],postPDSy[i],r=r[i],d=d[i],rc[i],dc[i],a[i],b[i],ac[i],bc[i],PDSM[i],PDSM1[i,l],M2[i],u[j],uc[j],randomM1,cure=Cure)
          SurvivalPDS[i] = PDSresult[1]+t+gap
          ChemocountPDS[i] = PDSresult[2]
        } else {
          SurvivalPDS[i] = t+gap
          ChemocountPDS[i] = 1
        }
        if (postNACTx[i]+postNACTy[i]<M2[i]) {
          NACTresult = Survive(postNACTx[i],postNACTy[i],r=r[i],d=d[i],rc[i],dc[i],a[i],b[i],ac[i],bc[i],NACTM[i],NACTM1[i,l],M2[i],u[j],uc[j],randomM1,cure=Cure)
          SurvivalNACT[i] = NACTresult[1]+t+gap+gap1
          ChemocountNACT[i] = NACTresult[2]
        } else {
          SurvivalNACT[i] = t+gap+gap1
          ChemocountNACT[i] = 1
        }
      }
      write.csv(cbind(postPDSx,postPDSy,postNACTx,postNACTy,SurvivalPDS,ChemocountPDS,SurvivalNACT,ChemocountNACT),file=paste("./Survival curve/M1",l," u",j," epsilon",k,".csv",sep=""))
    }
  }
}
