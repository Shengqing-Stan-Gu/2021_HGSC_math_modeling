setwd("~/Desktop/2014-5 Modeling Ovarian Progression/2015-4 Training parameters for the model/")

r = 2
d = 0.2
a = 0.84
b = 0.084
u = 10^(-7.6)
Diameter = seq(from=0.1,to=10,length.out=199)
M1 = Diameter^3/6*3.1416*10^9

py = vector("numeric")

for (i in 1:20) {
  x = 1:(M1[i]-1)
  pnoy = exp(-(1-b/a)*u*(x-1)/(1-d/r))
  pyes = pnoy*(1-exp(-(1-b/a)*u/(1-d/r)))
  py[i] = sum(pyes)
}

py[21:199] = 1

pdf("./Probability of resistance/Probability.pdf",width=5,height=4)
plot(Diameter,py,type="l",xlim=c(0,10),ylim=c(0,1),xlab="",ylab="")
dev.off()

write.csv(py,"./Probability of resistance/Probability.csv")
