procesogamma=function(c,beta,tiempo,N=20,b=1) 
{
  if (length(tiempo)==1)  tiempo=seq(0,tiempo,length=N)
  dt=c(0,diff(tiempo^b))
  x=0
  for (i in 2:length(tiempo))
  {
    x[i]=x[i-1]+rgamma(1,c*dt[i],scale=beta)
  }
  matrix(c(tiempo,x),ncol=2)
}