## Function to return survival time and number of chemo iterations

source("Transcendental_solution.R")
source("EstimateY.R")

Survive = function (xpost,ypost,r,d,rc,dc,a,b,ac,bc,M,M1,M2,u,uc,randomeM1,cure,deltat=0.005,t=4.2) {
  Survivaltime = 0
  Chemoround = 1
  
  if ((xpost==0)&&(ypost==0)) {
    return(c(cure,1))
  }
  
  else if (xpost==0) {
    Chemocount=1
    if (ypost>M1) {
      Survivaltime = log(M2/ypost)/(a-b)
      Chemofail = 1
    }
    else {
      Chemocount = 2
      Survivaltime = log(M1/ypost)/(a-b)
      if (M1*exp((ac-bc)*t)>M2) {
        Survivaltime = Survivaltime+log(M2/M1)/(ac-bc)
      }
      else {
        ypdsnum = M1*exp((ac-bc)*t)
        Survivaltime = Survivaltime+t+log(M2/ypdsnum)/(a-b)
      }
    }
    return(c(Survivaltime,Chemocount))
  }
  
  else if (ypost==0) {
    Chemofail = 0
    Chemocount = 1
    xpdsnum = vector("numeric")
    ypdsnum = vector("numeric")
    ypdsnum[1] = EstimateY(r,d,a,b,u,M1,deltatao=log(1.1)/(a-b),threshold=randomM1)
    xpdsnum[1] = M1-ypdsnum[1]
    Survivaltime = log(xpdsnum[1]/xpost)/(r-d)
    
    while(Chemofail==0) {
      total = xpdsnum[1]+ypdsnum[1]
      for (i in 2:(t/deltat)) {
        if (total>M2) {
          Chemofail=1
          break
        }
        else {
          xpdsnum[i] = xpdsnum[i-1]*(1+(rc*(1-uc)-dc)*deltat)
          ypdsnum[i] = ypdsnum[i-1]*(1+(ac-bc)*deltat)+xpdsnum[i-1]*rc*uc*deltat
          total = xpdsnum[i]+ypdsnum[i]
        }
      }
      Chemocount = Chemocount+1
      Survivaltime = Survivaltime+i*deltat
      if (Chemofail == 0) {
        if (total>(xpdsnum[1]+ypdsnum[1])) {
          Chemofail=1
          while(total<M2) {
            ypdsnum[i] = ypdsnum[i]*(1+(a-b)*deltat)+xpdsnum[i]*r*u*deltat
            xpdsnum[i] = xpdsnum[i]*(1+(r*(1-u)-d)*deltat)
            Survivaltime = Survivaltime+deltat
            total = xpdsnum[i]+ypdsnum[i]
          }
        }
        else {
          while (total<M1) {
            ypdsnum[i] = ypdsnum[i]*(1+(a-b)*deltat)+xpdsnum[i]*r*u*deltat
            xpdsnum[i] = xpdsnum[i]*(1+(r*(1-u)-d)*deltat)
            Survivaltime = Survivaltime+deltat
            total = xpdsnum[i]+ypdsnum[i]
          }
        }
      }
      if (Survivaltime>cure) {
        Chemofail=1
      }
      xpdsnum[1] = xpdsnum[i]
      ypdsnum[1] = ypdsnum[i]
    }
    return(c(Survivaltime,Chemocount))
  }
  
  else {
    Chemofail = 0
    Chemocount = 1
    xpdsnum = vector("numeric")
    ypdsnum = vector("numeric")
    xpdsnum[1] = xpost
    ypdsnum[1] = ypost
    total = xpdsnum[1] + ypdsnum[1]
    Survivaltime = 0
    if (total>M1) {
      Chemofail=1
      while(total<M2) {
        ypdsnum[1] = ypdsnum[1]*(1+(a-b)*deltat)+xpdsnum[1]*r*u*deltat
        xpdsnum[1] = xpdsnum[1]*(1+(r*(1-u)-d)*deltat)
        Survivaltime = Survivaltime+deltat
        total = xpdsnum[1]+ypdsnum[1]
      }
    }
    else {
      while(total<M1) {
        ypdsnum[1] = ypdsnum[1]*(1+(a-b)*deltat)+xpdsnum[1]*r*u*deltat
        xpdsnum[1] = xpdsnum[1]*(1+(r*(1-u)-d)*deltat)
        Survivaltime = Survivaltime+deltat
        total = xpdsnum[1]+ypdsnum[1]
      }
    }
    while(Chemofail==0) {
      total = xpdsnum[1]+ypdsnum[1]
      for (i in 2:(t/deltat)) {
        if (total>M2) {
          Chemofail=1
          break
        }
        else {
          xpdsnum[i] = xpdsnum[i-1]*(1+(rc*(1-uc)-dc)*deltat)
          ypdsnum[i] = ypdsnum[i-1]*(1+(ac-bc)*deltat)+xpdsnum[i-1]*rc*uc*deltat
          total = xpdsnum[i]+ypdsnum[i]
        }
      }
      Chemocount = Chemocount+1
      Survivaltime = Survivaltime+i*deltat
      if (Chemofail == 0) {
        if (total>(xpdsnum[1]+ypdsnum[1])) {
          Chemofail=1
          while(total<M2) {
            ypdsnum[i] = ypdsnum[i]*(1+(a-b)*deltat)+xpdsnum[i]*r*u*deltat
            xpdsnum[i] = xpdsnum[i]*(1+(r*(1-u)-d)*deltat)
            Survivaltime = Survivaltime+deltat
            total = xpdsnum[i]+ypdsnum[i]
          }
        }
        else {
          while (total<M1) {
            ypdsnum[i] = ypdsnum[i]*(1+(a-b)*deltat)+xpdsnum[i]*r*u*deltat
            xpdsnum[i] = xpdsnum[i]*(1+(r*(1-u)-d)*deltat)
            Survivaltime = Survivaltime+deltat
            total = xpdsnum[i]+ypdsnum[i]
          }
        }
      }
      if (Survivaltime>cure) {
        Chemofail=1
      }
      xpdsnum[1] = xpdsnum[i]
      ypdsnum[1] = ypdsnum[i]
    }
    return(c(Survivaltime,Chemocount))
  }
}