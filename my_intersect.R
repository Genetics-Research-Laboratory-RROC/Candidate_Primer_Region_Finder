my_intersect <- function(x,y){
  z = intersect(x,y)
  kx = match(z,x)
  ky = match(z,y)
  return (list(z,kx,ky))
}
  