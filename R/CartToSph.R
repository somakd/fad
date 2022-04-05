CartToSph <- function(x){
# Purpose: Convert cartesian with norm 1 to spherical coord system

x = x/sqrt(sum(x^2)) # adjust for small error;
P = length(x)
ang = rep(0,P-1)
ang[1] = acos(pmin(pmax(x[1],-1.0), 1.0))

for (p in 2:(P-1)){
  ang[p] = acos(pmin(pmax(x[p]/x[p-1]/tan(ang[p-1]),-1.0), 1.0))
}
if (x[P] < 0){
  ang[P-1] = 2*pi - ang[P-1]
}
return(ang)
}
