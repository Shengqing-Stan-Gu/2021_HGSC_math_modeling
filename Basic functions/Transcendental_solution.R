## Function to solve transcendental equation numerically, with error rate <5%
## Numerically solve t for x*exp((r-d)t) + exp((a-b)t) = M1

solvetd <- function(x,r,d,a,b,M1,randomM1) {
  lower = log(randomM1/(x+1))/max((r-d),(a-b))
  upper = log(randomM1/(x+1))/min((r-d),(a-b))
  middle = (lower+upper)/2
  M = x*exp((r-d)*middle)+exp((a-b)*middle)
  while (abs(M-randomM1)/randomM1 > 0.05) {
    if (M>randomM1) {
      upper = middle
      middle = (lower+upper)/2
      M = x*exp((r-d)*middle)+exp((a-b)*middle)
    }
    else if (M<randomM1) {
      lower = middle
      middle = (lower+upper)/2
      M = x*exp((r-d)*middle)+exp((a-b)*middle)
    }
  }
  return(middle)
}