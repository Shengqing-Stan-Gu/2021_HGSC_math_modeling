r = 2
d = 0.2
a = 0.84
b = 0.084
u = 10^(-7.6)
M1 = 10^11.5
randomM1 = 1e8
deltatao = log(1.1)/(a-b)

setwd("~/Desktop/2014-5 Modeling Ovarian Progression/2015-5-28 Training parameters for the model/")
source("Transcendental_solution.R")
source("EstimateY.R")

yexpected = vector("numeric")

## Make sure to create folder "Expected y at diagnosis" to store results

rseq = seq(from=1,to=3,length.out=41)
for (i in 1:length(rseq)) {
  yexpected[i] = EstimateY(r=rseq[i],d=d,a=a,b=b,u=u,M1=M1,deltatao=deltatao,threshold=randomM1)
}
pdf("./Y dependence/Expected y with r.pdf",width=4.5,height=4)
plot(rseq,yexpected,type="l",lwd=2,ylim=c(0,100000),xlab="",ylab="")
dev.off()
write.csv(cbind(rseq,yexpected),file="./Y dependence/Expected y with r.csv")
yexpected = read.csv("./Y dependence/Expected y with r.csv")[,3]

dseq = seq(from=0.1,to=1,length.out=41)
for (i in 1:length(dseq)) {
  yexpected[i] = EstimateY(r=r,d=dseq[i],a=a,b=b,u=u,M1=M1,deltatao=deltatao,threshold=randomM1)
}
pdf("./Y dependence/Expected y with d.pdf",width=4.5,height=4)
plot(dseq,yexpected,type="l",lwd=2,ylim=c(0,100000),xlab="",ylab="")
dev.off()
write.csv(cbind(dseq,yexpected),file="./Y dependence/Expected y with d.csv")
yexpected = read.csv("./Y dependence/Expected y with d.csv")[,3]

aseq = seq(from=0.4,to=1.6,length.out=41)
for (i in 1:length(aseq)) {
  yexpected[i] = EstimateY(r=r,d=d,a=aseq[i],b=b,u=u,M1=M1,deltatao=log(1.1)/(aseq[i]-b),threshold=randomM1)
}
pdf("./Y dependence/Expected y with a.pdf",width=4.5,height=4)
plot(aseq,yexpected,type="l",lwd=2,ylim=c(0,100000),xlab="",ylab="")
dev.off()
write.csv(cbind(aseq,yexpected),file="./Y dependence/Expected y with a.csv")
yexpected = read.csv("./Y dependence/Expected y with a.csv")[,3]

bseq = seq(from=0.042,to=0.42,length.out=41)
for (i in 1:length(bseq)) {
  yexpected[i] = EstimateY(r=r,d=d,a=a,b=bseq[i],u=u,M1=M1,deltatao=log(1.1)/(a-bseq[i]),threshold=randomM1)
}
pdf("./Y dependence/Expected y with b.pdf",width=4.5,height=4)
plot(bseq,yexpected,type="l",lwd=2,ylim=c(0,100000),xlab="",ylab="")
dev.off()
write.csv(cbind(bseq,yexpected),file="./Y dependence/Expected y with b.csv")
yexpected = read.csv("./Y dependence/Expected y with b.csv")[,3]

useq = 10^(seq(from=-9,to=-7,length.out=41))
for (i in 1:length(useq)) {
  yexpected[i] = EstimateY(r=r,d=d,a=a,b=b,u=useq[i],M1=M1,deltatao=deltatao,threshold=randomM1)
}
pdf("./Y dependence/Expected y with u.pdf",width=4.5,height=4)
plot(log10(useq),yexpected,type="l",lwd=2,ylim=c(0,100000),xlab="",ylab="")
dev.off()
write.csv(cbind(useq,yexpected),file="./Y dependence/Expected y with u.csv")
yexpected = read.csv("./Y dependence/Expected y with u.csv")[,3]

M1seq = 10^(seq(from=11,to=12.5,length.out=41))
for (i in 1:length(M1seq)) {
  yexpected[i] = EstimateY(r=r,d=d,a=a,b=b,u=u,M1=M1seq[i],deltatao=deltatao,threshold=randomM1)
}
pdf("./Y dependence/Expected y with M1.pdf",width=4.5,height=4)
plot(log10(M1seq),yexpected,type="l",lwd=2,ylim=c(0,100000),xlab="",ylab="")
dev.off()
write.csv(cbind(M1seq,yexpected),file="./Y dependence/Expected y with M1.csv")
yexpected = read.csv("./Y dependence/Expected y with M1.csv")[,3]
