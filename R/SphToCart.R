SphToCart <- function(ang){
# Purpose: Convert spherical coords to cartesian system

P = length(ang)+1
x = rep(0,P)
x[1] = cos(ang[1])

for (p in 2:(P-1)){
  x[p] = tan(ang[p-1])*cos(ang[p])*x[p-1]
}
x[P] = x[P-1]*tan(ang[P-1])
return(x)
}