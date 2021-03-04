## Function to calculate expected number of resistant cells at diagnosis, "y"
## threshold means the time with number of x cells, by which time we still consider random emergence of y

source("Transcendental_solution.R")

EstimateY = function(r,d,a,b,u,M1,deltatao,threshold) {
	if (M1<threshold) {
		threshold = M1
	}
	x = 1:(threshold-1)
	pnoy = exp(-(1-b/a)*u*(x-1)/(1-d/r))
	pyes = pnoy*(1-exp(-(1-b/a)*u/(1-d/r)))
	tao = vector("numeric")
	i = 1
	while (i < threshold) {
  		if (i*exp((r-d)*deltatao)<i+1) {
    		tao[i] = solvetd(x=i,r,d,a,b,M1,threshold)
    		i = i+1
  		}
  		else if (i*exp((r-d)*deltatao)>threshold-1) {
    		t = solvetd(x=(i+threshold-1)/2,r,d,a,b,M1,threshold)
    		tao[i:(threshold-1)] = t
    		i = threshold+1
  		}
  		else if (i*exp((r-d)*deltatao)>i+1) {
    		j = floor(i*exp((r-d)*deltatao))
    		t = solvetd(x=(i+j)/2,r,d,a,b,M1,threshold)
    		tao[i:j] = t
    		i = j+1
  		}
	}
	yexpect = sum(pyes*exp((a-b)*tao))
	total = threshold
	xexpect = total-yexpect
	dt = 0.005
	while (total<M1) {
		yexpect = yexpect+yexpect*(a-b)*dt+r*u*xexpect*(1-b/a)*dt
		xexpect = xexpect*(1+(r*(1-u)-d)*dt)
		total = yexpect+xexpect
	}
	return(yexpect)
}