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
randomM1 = 1e8
deltatao = log(1.1)/(a-b)
Cure = 120

setwd("~/Desktop/2014-5 Modeling Ovarian Progression/2015-6-2 Model alternative strategies/")
source("Transcendental_solution.R")
source("EstimateY.R")

yexpectedPDS = vector("numeric")
yexpectedNACT = vector("numeric")

## Make sure to create folder "Expected y at diagnosis" to store results

for (i in 1:N) {
  yexpectedPDS[i] = EstimateY(r=r[i],d=d[i],a=a[i],b=b[i],u=u,M1=PDSM1[i],deltatao=deltatao[i],threshold=randomM1)
  yexpectedNACT[i] = EstimateY(r=r[i],d=d[i],a=a[i],b=b[i],u=u,M1=NACTM1[i],deltatao=deltatao[i],threshold=randomM1)
}


write.csv(yexpectedPDS,file="./Expected y at diagnosis/yexpectedPDS.csv")
write.csv(yexpectedNACT,file="./Expected y at diagnosis/yexpectedNACT.csv")
write.csv(cbind(r,a,rc,dc,PDSM1,NACTM1,M2),file="./Expected y at diagnosis/parameters used.csv")
