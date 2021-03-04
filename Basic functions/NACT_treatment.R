## Function that returns the number of sensitive (x) and resistant (y) cells after PDS treatment.

source("Transcendental_solution.R")
source("EstimateY.R")

NACTtreat = function (r,d,a,b,rc,dc,ac,bc,epsilon,M,M1,u,uc,y0,deltat=0.005,t=4.2,cycle=2.1,gap=1,gap1=0.3) {
    xnact = vector("numeric")
  	ynact = vector("numeric")
  	xnact[1] = M1-y0
  	ynact[1] = y0
    cycle = ifelse(dc<5.2,2.805,2.1)
  	## Model the period of neo-adjuvant chemo
    for (j in 2:(cycle/deltat)) {
  		xnact[j] = xnact[j-1]*(1+(rc*(1-uc)-dc)*deltat)
  		ynact[j] = (1+(ac-bc)*deltat)*ynact[j-1]+deltat*rc*uc*xnact[j-1]
  	}
    ## Model the period of waiting time between chemo and surgery
    for (j in (cycle/deltat+1):((cycle+gap1)/deltat)) {
      xnact[j] = xnact[j-1]*(1+(r*(1-u)-d)*deltat)
      ynact[j] = (1+(a-b)*deltat)*ynact[j-1]+deltat*r*u*xnact[j-1]
    }
    ## Model the effect of surgery
    ## ynact0 here will be used to calculate the probability that chemo-resistant are depleted
  	xnact0 = xnact[j]/(xnact[j]+ynact[j])*M + epsilon*xnact[j]
  	ynact0 = ynact[j]/(xnact[j]+ynact[j])*M + epsilon*ynact[j]
    ## Model the period of waiting time between surgery and adjuvant chemo
    ## From here till the end of adjuvant chemo, ynact will mean the resistant cells that are converted from sensitive cells
  	xnact[j+1] = xnact0+deltat*(r*(1-u)-d)*xnact0
  	ynact[j+1] = deltat*rc*uc*xnact0
    for (j in ((cycle+gap1)/deltat+2):((cycle+gap1+gap)/deltat)) {
      xnact[j] = xnact[j-1]*(1+(r*(1-u)-d)*deltat)
      ynact[j] = ynact[j-1]+deltat*r*u*xnact[j-1]
    }
    ## Model the period of adjuvant chemo
  	for (j in ((cycle+gap1+gap)/deltat+1):((cycle+gap1+gap+t/2)/deltat)) {
    	xnact[j] = xnact[j-1]*(1+(rc*(1-uc)-dc)*deltat)
  		ynact[j] = ynact[j-1]+deltat*rc*uc*xnact[j-1]
  	}
  	ynact[j] = ynact[j]*(1-bc/ac)
  	
  	## Probability of no y right after debulking
  	pnoydeb = exp(-ynact0)
  	## Probability of y persist right after debulking
  	pwydeb = 1-pnoydeb
  	## Probability of eradication of x by chemo
  	pnoxchemo = exp(-xnact[j])
  	## Probability of x persist through chemo
  	pwxchemo = 1-pnoxchemo
  	## Probability of no mutation during chemo
  	pnoychemo = exp(-ynact[j])
  	## Probability of mutation during chemo
  	pwychemo = 1-pnoychemo
  	
  	p = runif(n=1,min=0,max=1)
  	
    ## No tumor left
  	if (p<pnoydeb*pnoxchemo*pnoychemo) {
  		return(c(0,0))
  	}
    ## No x left, no y from surgery, but y generated during adjuvant chemo
  	else if (p<pnoydeb*pnoxchemo) {
  		return(c(0,ynact[j]*exp((ac-bc)*(t-cycle+gap)/2)/pwychemo))
  	}
    ## No y, but x is left
  	else if (p<pnoydeb*pnoxchemo+pnoydeb*pwxchemo*pnoychemo) {
  		return(c(xnact[j]/pwxchemo,0))
  	}
    ## x is left, no y from surgery, but y generated during adjuvant chemo
  	else if (p<pnoydeb) {
  		return(c(xnact[j]/pwxchemo,ynact[j]*exp((ac-bc)*(t-cycle+gap)/2)/pwychemo))
  	}
    ## y is left from debulking, no x
  	else if (p<pnoydeb+pwydeb*pnoxchemo) {
  		return(c(0,ynact0/pwydeb*exp((a-b)*gap)*exp((ac-bc)*(t-cycle))+ynact[j]*exp((ac-bc)*(t-cycle+gap)/2)))
  	}
    ## y is left from debulking, and x is left after chemo
  	else {
  		return(c(xnact[j]/pwxchemo,ynact0/pwydeb*exp((a-b)*gap)*exp((ac-bc)*(t-cycle))+ynact[j]*exp((ac-bc)*(t-cycle+gap)/2)))
  	}
}

