# Análisis con b conocido
# NOTA: Las funciones necesarias se encuentran al final del 
# script

library(stats4) # mle - Estimación método de máxima
#                verosimilitud
library(lme4) # Nelder_Mead
library(GenSA) # GenSA - Simulated annealing con  
#               restricciones para los parámetros
library(GA) # ga - Algoritmo genético con restricciones para  
#            los parámetros

data=as.matrix(read.table("fatiguecrack_data.txt",sep=','))

colnames(data)=NULL
t=data[1,]
data=data[2:nrow(data),]


delta=t(apply(data,1,diff))
b=1

w=diff(t^b)

n_params=2

# Representación de las trayectorias de deterioro
{
  par(ask=FALSE)
  cols=rainbow(nrow(data)-1)
  plot(t,data[1,],type='b',ylim=c(0,max(data)*1.25),lwd=2.5,
       xlab="Cientos de miles de ciclos",
       ylab="Deterioro acumulado")
  for (i in 2:nrow(data))
  {
    lines(t,data[i,],type='b',lwd=2.5,col=cols[i-1])
  }
  abline(h=0.4,lwd=3)
  text(0.15, 0.43, "Nivel crítico de deterioro",cex=0.85)
  legend("topleft",legend=paste("Dispositivo",1:nrow(data)),
         ncol=2,y.intersp=0.5,cex=0.8,col=c("black",cols),
         lty=1,bty ="n")
}


# Método de los momentos - válido cuando b=1
if (b==1)
{
  beta_MM=numeric()
  c_MM=numeric()
  for (i in 1:nrow(data))
  {
    beta_MM[i]=sum( (delta[i,]-data[i,ncol(data)]*w/
                      t[length(t)]^b )^2) /
               (data[i,ncol(data)]*(1-sum(w^2)/(sum(w)^2)))
    c_MM[i]=data[i,ncol(data)]/(t[length(t)]^b *beta_MM[i])
  }
  
  # Comparamos los resultados gráficamente simulando un  
  # proceso gamma con los parámetros estimados
  par(ask=TRUE)
  p.MM=matrix(nrow=nrow(data),ncol=ncol(data))
  for (i in 1:nrow(data))
  {
    p.MM[i,]=procesogamma(c_MM[i],beta_MM[i],t,b=b)[,2]
    plot(t,data[i,],type='s', ylim=c(0,0.6), lwd=3,
         main='Estimación por el Método de los Momentos',
         xlab="Cientos de miles de ciclos",
         ylab="Deterioro acumulado (mm)")
    lines(t,p.MM[i,],type='s',lwd=2.5,col='red')
    legend('topleft',
           legend=c(paste('Datos dispositivo', i),
                    'Simulación con los parámetros estimados'),
           col=c('black','red'),lty=1)
  }
}



# Método de máxima verosimilitud
c_MV=numeric()
beta_MV=numeric()
for (i in 1:nrow(data))
{
  EMV=mle(NegLogLikb,start=list(c=1,beta=1),
          fixed=list(b=b,fila=i))
  c_MV[i]=as.numeric(coef(EMV)[1])
  beta_MV[i]=as.numeric(coef(EMV)[2])
}


# Comparamos los resultados gráficamente simulando un proceso
# gamma con los parámetros estimados
par(ask=TRUE)
p.MV=matrix(nrow=nrow(data),ncol=ncol(data))
for (i in 1:nrow(data))
{
  p.MV[i,]=procesogamma(c_MV[i],beta_MV[i],t,b=b)[,2]
  
  plot(t,data[i,],type='s',lwd=3, ylim=c(0,0.6),
       main='Estimación por el Método de Máxima Verosimilitud',
       xlab="Cientos de miles de ciclos",
       ylab="Deterioro acumulado (mm)")
  lines(t,p.MV[i,],type='s',lwd=2.5,col='red')
  legend('topleft',
         legend=c(paste('Datos dispositivo', i),
                  'Simulación con los parámetros estimados'),
         col=c('black','red'),lty=1)
}


# Métodos de optimización
# ag - Algoritmo genético
# nm - Método de Nelder-Mead
# sa - Simulated Annealing
# Optimizamos cada fila
opt_ag=matrix(,nrow=nrow(data),ncol=n_params)
opt_nm=matrix(,nrow=nrow(data),ncol=n_params)
opt_sa=matrix(,nrow=nrow(data),ncol=n_params)
dispositivo=matrix(,nrow=nrow(data),ncol=3*n_params)
for (fila in 1:nrow(data))
{
  # Algoritmo genético
  aux=ga(type='real-valued',fitness=LogLik,
         lower=rep(0.001,n_params),upper=rep(1000,n_params),
         optim=TRUE)
  opt_ag[fila,]=aux@solution
  
  # Nelder-Mead
  opt_nm[fila,]=Nelder_Mead(function(x) -LogLik(x),
                            rep(5,n_params),
                            lower=rep(0.001,n_params))$par
  
  # Simulated annealing
  opt_sa[fila,]=GenSA(fn=function(x) -Lik(x),
                      lower=rep(0.001,n_params),
                       upper=rep(100,n_params),
                      control=list(max.time=30))$par
  
  dispositivo[fila,]=c(opt_ag[fila,],opt_nm[fila,],opt_sa[fila,])
}

# Una vez estimados los parámetros simulamos un proceso gamma 
# por cada uno de los métodos de optimización y comparamos 
# gráficamente los resultados
cols=rainbow(5)
p.opt=numeric() 
      # columna 1 - dispositivo / columna 2 - método
                                # (1-ag / 2-nm / 3-sa)
for (i in 1:nrow(data))
{
  plot(t,data[i,],type='s',ylim=c(0,max(data)*1.25),lwd=3,
       main=paste("Dispositivo",i),
       xlab="Cientos de miles de ciclos",
       ylab="Deterioro acumulado (mm)")
  for (j in seq(1,3*n_params,n_params))
  {
    p=procesogamma(dispositivo[i,j],dispositivo[i,j+1],t,
                   length(t))[,2]
    lines(t,p,type='s',lwd=2.5,col=cols[ceiling(j/n_params)])
    p.opt=rbind(p.opt,c(i,ceiling(j/n_params),p))
  }
  
  legend('topleft',legend=c(paste('Datos dispositivo', i),
                            'Algoritmo genético',
                            'Método Nelder-Mead',
                            'Simulated annealing'),
         col=c('black',cols),lty=1,cex=0.75)
}

# Para comparar los métodos de los momentos y máxima 
# verosimilitud gráficamente con los métodos de optimización
cols=rainbow(5)
for (i in 1:nrow(data))
{
  plot(t,data[i,],type='s',ylim=c(0,max(data)*1.25),lwd=3,
       main=paste("Dispositivo",i),
       xlab="Cientos de miles de ciclos",
       ylab="Deterioro acumulado (mm)")
  for (j in seq(1,3*n_params,n_params))
  {
    aux=p.opt[(p.opt[,1]==i),][,3:ncol(p.opt)] #pgamma 
                                              # dispositivo i
    lines(t,aux[ceiling(j/n_params),],type='s',lwd=2.5,
          col=cols[ceiling(j/n_params)]) 
  }
  
  lines(t,p.MV[i,],type='s',lwd=2.5,col=cols[4])
  if (b==1) 
  {
    lines(t,p.MM[i,],type='s',lwd=2.5,col=cols[5])
    legend('topleft',
           legend=c(paste('Datos dispositivo', i),
                    'Algoritmo genético',
                    'Método Nelder-Mead',
                    'Simulated annealing',
                    'Método de máxima verosimilitud',
                    'Método de los momentos'),
           col=c('black',cols),lty=1,cex=0.75)
  } 
  else
  {
    legend('topleft',
           legend=c(paste('Datos dispositivo', i),
                    'Algoritmo genético',
                    'Método Nelder-Mead',
                    'Simulated annealing',
                    'Método de máxima verosimilitud'),
           col=c('black',cols[1:4]),lty=1,cex=0.75)
  }
}


# Comparamos los métodos mediante el criterio de información
# de Akaike con la corrección para muestras finitas
AIC=matrix(,nrow=nrow(data),ncol=5)
mejor=matrix(,nrow=nrow(data),ncol=n_params+1)
for (fila in 1:nrow(data))
{
  k=1
  for (j in seq(1,3*n_params,n_params))
  {
    AIC[fila,k]=2*n_params-2*log(Lik(c(dispositivo[fila,j],
                                       dispositivo[fila,j+1])))
    AIC[fila,k]=AIC[fila,k]+(2*n_params^2+2*n_params) /
                      (ncol(data)-n_params-1)
    k=k+1
  }
  
  AIC[fila,4]=2*n_params-2*log(Lik(c(c_MV[fila],beta_MV[fila])))
  AIC[fila,4]=AIC[fila,4]+(2*n_params^2+2*n_params) /
                    (ncol(data)-n_params-1)
  if (b==1)
  {
    AIC[fila,5]=2*n_params-2*log(Lik(c(c_MM[fila],beta_MM[fila])))
    AIC[fila,5]=AIC[fila,5]+(2*n_params^2+2*n_params) /
                      (ncol(data)-n_params-1)
  }
  
  
  # El mejor modelo para cada dispositivo
  aux=which.min(AIC[fila,])
  if (aux<4)
  {
    mejor[fila,]=c(dispositivo[fila,aux*n_params-1],
                dispositivo[fila,aux*n_params],aux)
  }
  else if (aux==4) 
  {
    mejor[fila,]=c(c_MV[fila],beta_MV[fila],aux)
  }
  else if (aux==5)
  {
    mejor[fila,]=c(c_MM[fila],beta_MM[fila],aux)
  }
}
aux=factor(x=mejor[,3],levels=c(1:5),
           labels=c('Algoritmo genético','Nelder-Mead',
                    'Simulated annealing',
                    'Método de máxima verosimilitud',
                    'Método de los momentos'))
mejor=data.frame(c=mejor[,1],beta=mejor[,2],Método=aux)

# Representamos para cada dispositivo los datos medidos junto
# con una simulación con los parámetros obtenidos con el 
# mejor método
cols=rainbow(5)
for (i in 1:nrow(data))
{
  plot(t,data[i,],type='s',ylim=c(0,max(data)*1.25),lwd=3,
       main=paste("Dispositivo",i),
       xlab="Cientos de miles de ciclos",
       ylab="Deterioro acumulado (mm)")
  switch (as.character(mejor[i,3]),
    'Algoritmo genético' = 
      {p=p.opt[(p.opt[,1]==i),][,3:ncol(p.opt)] # pgamma 
                                              # dispositivo i
       lines(t,p[1,],type='s',lwd=2.5,col=cols[1])
       aux=1},
    'Nelder-Mead' = 
      {p=p.opt[(p.opt[,1]==i),][,3:ncol(p.opt)] # pgamma 
                                              # dispositivo i
       lines(t,p[2,],type='s',lwd=2.5,col=cols[2])
       aux=2},
    'Simulated annealing' = 
      {p=p.opt[(p.opt[,1]==i),][,3:ncol(p.opt)] # pgamma 
                                              # dispositivo i
       lines(t,p[3,],type='s',lwd=2.5,col=cols[3])
       aux=3},
    'Método de máxima verosimilitud' = 
      {lines(t,p.MV[i,],type='s',lwd=2.5,col=cols[4])
        aux=4},
    'Método de los momentos' = 
      {lines(t,p.MM[i,],type='s',lwd=2.5,col=cols[5])
        aux=5}
  )
  
  legend('topleft',
         legend=c(paste('Datos dispositivo', i),
                  as.character(mejor[i,3])),
         col=c('black',cols[aux]),lty=1,cex=0.75)
}

# Escribimos los parámetros estimados de los modelos con menor
# AIC en un archivo de texto
write.table(mejor[,1:2],file='fc_mejor.txt',sep='\t',
            col.names=F,row.names=F)

#-----------------------------------------------------------

##------------------FUNCIONES-------------------------------
#-----------------------------------------------------------
# Simulación de un proceso gamma
# t puede ser tanto el tiempo máximo como el vector de 
# tiempos
procesogamma=function(c,beta,t,N=10,b=1) 
{
  if (length(t)==1)  t=seq(0,t,length=N)
  dt=c(0,diff(t^b))
  x=0
  for (i in 2:length(t))
  {
    x[i]=x[i-1]+rgamma(1,c*dt[i],scale=beta)
  }
  matrix(c(t,x),ncol=2)
}


# Funciones de verosimilitud y negativas de verosimilitud y 
# logverosimilitud que queremos minimizar o maximizar
Lik=function(x)
{
  c=x[1]
  beta=x[2]
  if (length(x)==3) b=x[3]
  else b=1
  ( prod( delta[fila,]^(c*(diff(t^b))-1)*
            exp(-delta[fila,]/beta) /
         (gamma(c*diff(t^b))*beta^(c*diff(t^b))) ) )
}

LogLik=function(x)
{
  c=x[1]
  beta=x[2]
  if (length(x)==3) b=x[3]
  else b=1
 ( sum(c*(diff(t^b))*log(delta[fila,]/(beta+0.0001))) 
    - sum(log(delta[fila,]*gamma(c*diff(t^b)))) 
    - data[fila,ncol(data)]/(beta+0.0001) )
}

# b no conocido
NegLogLikb=function(c,beta,b=1,fila=1)
{
  -(sum(c*(diff(t^b))*log(delta[fila,]/(beta+0.0001))) 
    - sum(log(delta[fila,]*gamma(c*(diff(t^b))))) 
    - data[fila,ncol(data)]/(beta+0.0001) )
}




