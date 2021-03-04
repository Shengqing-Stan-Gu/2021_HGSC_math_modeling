set.seed(100)
N = 100
r = abs(rnorm(N,mean=2,sd=0.8))
d = r/10
rc = abs(rnorm(N,mean=0.2,sd=0.08))
dc = abs(rnorm(N,mean=4.9,sd=1))
a = abs(rnorm(N,mean=0.84,sd=0.42))
b = a/10
ac = a
bc = ac/5

epsilon = 10^(-7.4)
PDSM1 = 10^rnorm(N,mean=11.5,sd=0.4)
NACTM1 = 10^rnorm(N,mean=12,sd=0.4)
M2 = 10^rnorm(N,mean=13,sd=0.4)
PDSM = PDSM1/(10^5.5)
NACTM = NACTM1/(10^6)
u = 10^(-7.6)
uc = 10*u
deltat = 0.005
t = 4.2
gap = 1
gap1 = 0.3
gap2 = 1.5
gap3 = 2
gap4 = 3
randomM1 = 1e8
deltatao = log(1.1)/(a-b)
Cure = 120

cycle = 2.1

setwd("~/Desktop/2014-5 Modeling Ovarian Progression/2015-6-2 Model alternative strategies/")
source("Transcendental_solution.R")
source("EstimateY.R")
source("PDS_treatment.R")
source("NACT_treatment.R")
source("Survivaltime.R")

y0PDS = read.csv("./Expected y at diagnosis/yexpectedPDS.csv")[,-1]
y0NACT = read.csv("./Expected y at diagnosis/yexpectedNACT.csv")[,-1]
postPDSx = vector("numeric")
postPDSy = vector("numeric")
postNACTx = vector("numeric")
postNACTy = vector("numeric")
SurvivalPDS = vector("numeric")
ChemocountPDS = vector("numeric")
SurvivalNACT = vector("numeric")
ChemocountNACT = vector("numeric")

for (i in 1:N) {
  if (y0NACT[i]*exp((ac[i]-bc[i])*t/2+(a[i]-b[i])*gap1)<M2[i]) {
    postNACT = NACTtreat(r=r[i],d=d[i],a=a[i],b=b[i],rc=rc[i],dc=dc[i],ac=ac[i],bc=bc[i],epsilon=epsilon,M=NACTM[i],M1=NACTM1[i],u=u,uc=uc,y0=y0NACT[i],cycle=cycle,gap=gap)
  } else {
    postNACT = c(0,y0NACT[i]*exp((ac[i]-bc[i])*t/2+(a[i]-b[i])*gap1))
  }
  postNACTx[i] = postNACT[1]
  postNACTy[i] = postNACT[2]
  
  if (postNACTx[i]+postNACTy[i]<M2[i]) {
    NACTresult = Survive(postNACTx[i],postNACTy[i],r=r[i],d=d[i],rc[i],dc[i],a[i],b[i],ac[i],bc[i],PDSM[i],PDSM1[i],M2[i],u,uc,randomM1,cure=Cure)
    SurvivalNACT[i] = NACTresult[1]+t+gap+gap1
    ChemocountNACT[i] = NACTresult[2]
  } else {
    SurvivalNACT[i] = t+gap+gap1
    ChemocountNACT[i] = 1
  }
}

write.csv(cbind(postNACTx,postNACTy,SurvivalNACT,ChemocountNACT),file="./Survival NACT chemo delay/Survival 1 month optimal.csv")


for (i in 1:N) {
  if (y0NACT[i]*exp((ac[i]-bc[i])*t/2+(a[i]-b[i])*gap1)<M2[i]) {
    postNACT = NACTtreat(r=r[i],d=d[i],a=a[i],b=b[i],rc=rc[i],dc=dc[i],ac=ac[i],bc=bc[i],epsilon=epsilon,M=NACTM[i],M1=NACTM1[i],u=u,uc=uc,y0=y0NACT[i],cycle=cycle,gap=gap2)
  } else {
    postNACT = c(0,y0NACT[i]*exp((ac[i]-bc[i])*t/2+(a[i]-b[i])*gap1))
  }
  postNACTx[i] = postNACT[1]
  postNACTy[i] = postNACT[2]
  
  if (postNACTx[i]+postNACTy[i]<M2[i]) {
    NACTresult = Survive(postNACTx[i],postNACTy[i],r=r[i],d=d[i],rc[i],dc[i],a[i],b[i],ac[i],bc[i],PDSM[i],PDSM1[i],M2[i],u,uc,randomM1,cure=Cure)
    SurvivalNACT[i] = NACTresult[1]+t+gap+gap1
    ChemocountNACT[i] = NACTresult[2]
  } else {
    SurvivalNACT[i] = t+gap2+gap1
    ChemocountNACT[i] = 1
  }
}

write.csv(cbind(postNACTx,postNACTy,SurvivalNACT,ChemocountNACT),file="./Survival NACT chemo delay/Survival 1.5 month optimal.csv")


for (i in 1:N) {
  if (y0NACT[i]*exp((ac[i]-bc[i])*t/2+(a[i]-b[i])*gap1)<M2[i]) {
    postNACT = NACTtreat(r=r[i],d=d[i],a=a[i],b=b[i],rc=rc[i],dc=dc[i],ac=ac[i],bc=bc[i],epsilon=epsilon,M=NACTM[i],M1=NACTM1[i],u=u,uc=uc,y0=y0NACT[i],cycle=cycle,gap=gap3)
  } else {
    postNACT = c(0,y0NACT[i]*exp((ac[i]-bc[i])*t/2+(a[i]-b[i])*gap1))
  }
  postNACTx[i] = postNACT[1]
  postNACTy[i] = postNACT[2]
  
  if (postNACTx[i]+postNACTy[i]<M2[i]) {
    NACTresult = Survive(postNACTx[i],postNACTy[i],r=r[i],d=d[i],rc[i],dc[i],a[i],b[i],ac[i],bc[i],PDSM[i],PDSM1[i],M2[i],u,uc,randomM1,cure=Cure)
    SurvivalNACT[i] = NACTresult[1]+t+gap+gap1
    ChemocountNACT[i] = NACTresult[2]
  } else {
    SurvivalNACT[i] = t+gap3+gap1
    ChemocountNACT[i] = 1
  }
}

write.csv(cbind(postNACTx,postNACTy,SurvivalNACT,ChemocountNACT),file="./Survival NACT chemo delay/Survival 2 month optimal.csv")


for (i in 1:N) {
  if (y0NACT[i]*exp((ac[i]-bc[i])*t/2+(a[i]-b[i])*gap1)<M2[i]) {
    postNACT = NACTtreat(r=r[i],d=d[i],a=a[i],b=b[i],rc=rc[i],dc=dc[i],ac=ac[i],bc=bc[i],epsilon=epsilon,M=NACTM[i],M1=NACTM1[i],u=u,uc=uc,y0=y0NACT[i],cycle=cycle,gap=gap4)
  } else {
    postNACT = c(0,y0NACT[i]*exp((ac[i]-bc[i])*t/2+(a[i]-b[i])*gap1))
  }
  postNACTx[i] = postNACT[1]
  postNACTy[i] = postNACT[2]
  
  if (postNACTx[i]+postNACTy[i]<M2[i]) {
    NACTresult = Survive(postNACTx[i],postNACTy[i],r=r[i],d=d[i],rc[i],dc[i],a[i],b[i],ac[i],bc[i],PDSM[i],PDSM1[i],M2[i],u,uc,randomM1,cure=Cure)
    SurvivalNACT[i] = NACTresult[1]+t+gap+gap1
    ChemocountNACT[i] = NACTresult[2]
  } else {
    SurvivalNACT[i] = t+gap4+gap1
    ChemocountNACT[i] = 1
  }
}

write.csv(cbind(postNACTx,postNACTy,SurvivalNACT,ChemocountNACT),file="./Survival NACT chemo delay/Survival 3 month optimal.csv")


