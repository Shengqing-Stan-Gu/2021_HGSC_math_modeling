## Function that returns the number of sensitive (x) and resistant (y) cells after PDS treatment.

source("Transcendental_solution.R")
source("EstimateY.R")

PDStreat = function (r,d,a,b,rc,dc,ac,bc,epsilon,M,M1,u,uc,y0,deltat=0.005,t=4.2,gap=1) {
	  x0 = M1-y0
	  xpds = vector("numeric")
  	ypds = vector("numeric")
    
	  ## Model the effect of surgery
	  ## ypds0 here will be used to calculate the probability that chemo-resistant are depleted by surgery
    xpds0 = x0/M1*(M+epsilon*M1)
    ypds0 = y0/M1*(M+epsilon*M1)
	  ## Model the period of waiting time between surgery and adjuvant chemo
	  ## From here till the end of adjuvant chemo, ypds will mean the resistant cells that are converted from sensitive cells
    xpds[1] = xpds0
    ypds[1] = 0
    for (j in 2:(gap/deltat)) {
      xpds[j] = xpds[j-1]*(1+(r*(1-u)-d)*deltat)
      ypds[j] = ypds[j-1]+deltat*r*u*xpds[j-1]
    }
    ## Model the period of adjuvant chemo
    for (j in (gap/deltat+1):((t+gap)/deltat)) {
      xpds[j] = xpds[j-1]*(1+(rc*(1-uc)-dc)*deltat)
      ypds[j] = ypds[j-1]+deltat*rc*uc*xpds[j-1]
    }
  	ypds[j] = ypds[j]*(1-bc/ac)
  	
  	## Probability of no y right after debulking
  	pnoydeb = exp(-ypds0)
  	## Probability of y persist right after debulking
  	pwydeb = 1-pnoydeb
  	## Probability of eradication of x by chemo
  	pnoxchemo = exp(-xpds[j])
  	## Probability of x persist through chemo
  	pwxchemo = 1-pnoxchemo
  	## Probability of no mutation during chemo
  	pnoychemo = exp(-ypds[j])
  	## Probability of mutation during chemo
  	pwychemo = 1-pnoychemo
  	
  	p = runif(n=1,min=0,max=1)
  	
    ## No tumor left
  	if (p<pnoydeb*pnoxchemo*pnoychemo) {
  		return(c(0,0))
  	}
    ## No x left, no y from surgery, but y generated during adjuvant chemo
  	else if (p<pnoydeb*pnoxchemo) {
  		return(c(0,ypds[j]*exp((ac-bc)*(t+gap)/2)/pwychemo))
  	}
	  ## No y, but x is left
  	else if (p<pnoydeb*pnoxchemo+pnoydeb*pwxchemo*pnoychemo) {
  		return(c(xpds[j]/pwxchemo,0))
  	}
	  ## x is left, no y from surgery, but y generated during adjuvant chemo
  	else if (p<pnoydeb) {
  		return(c(xpds[j]/pwxchemo,ypds[j]*exp((ac-bc)*(t+gap)/2)/pwychemo))
  	}
    ## y is left from debulking, no x
  	else if (p<pnoydeb+pwydeb*pnoxchemo) {
  		return(c(0,ypds0/pwydeb*exp((a-b)*gap)*exp((ac-bc)*t)+ypds[j]*exp((ac-bc)*(t+gap)/2)))
  	}
    ## y is left from debulking, and x is left after chemo
  	else {
  		return(c(xpds[j]/pwxchemo,ypds0/pwydeb*exp((a-b)*gap)*exp((ac-bc)*t)+ypds[j]*exp((ac-bc)*(t+gap)/2)))
  	}
}

