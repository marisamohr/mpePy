
simMFBM <- function(n=50, H=c(.2,.2),
                   sig=c(1,1),rho=matrix(c(1,.5,.5,1),nc=2),eta=matrix(c(0,.5,-.5,0),nc=2),
                   plot=FALSE,print=TRUE,choix=NULL,forceEta=FALSE,nbSamp=1){ 
  
  ## Simulation of a multivariate fractional Brownian motion.
  ## nbSamp : number of sample paths if =1 a (n,p) matrix is returned
  ## if >1 a list of (n,p) matrices is returned 
  ##  notation: causal, wb, eta stand for the causal MFBM, the well-balanced MFBM, 
  ## a general MFBM with prescribed eta matrix. 
  ## forceEta=TRUE means that eta_jk is changed by eta_jk/(1-H[j]-H[k]) only when H[j]+H[k] !=1
  ## so that this is easy to check the semidefinite positiveness using ellipsis in 
  ## See the reference for details.
  ##
  ## Reference. Amblard, Coeurjolly, Philippe and Lavancier (2012)
  ## Basic properties of the multivariate fractional Brownian motion.
  ## Bulletin de la SMF, SÃ©minaires et CongrÃ©s, 28, 65-87.
  
  p<-length(H)
  G<-matrix(0,nr=p,nc=p)
  if (is.null(choix)) choix<-"eta"
  switch(choix,
         "eta"={
           if (forceEta){
             for (j in (1:p)) { for (k in (1:p)){
               if((H[j]+H[k])!=1 ) eta[j,k]<-eta[j,k]/(1-H[j]-H[k])
               
             }}
           }
         }, 
         "causal"={
           eta<-matrix(0,nc=p,nr=p)
           for (j in (1:p)) 
             for (k in (1:p)){
               Hjk<-H[j]+H[k]
               if (Hjk!=1){ eta[j,k]<- rho[j,k]*(cos(pi*H[j])-cos(pi*H[k]))/(cos(pi*H[j])+cos(pi*H[k]))	}
               else{ eta[j,k]<-rho[j,k]*2/pi/tan(pi*H[j])}
             }
           diag(eta)<-0
         },
         "wb"={
           eta<-matrix(0,nc=p,nr=p)
         })	
  
  ## computation of the matrix G (existence proposition)
  for (j in (1:p)){for (k in (1:p)){
    Hjk<-H[j]+H[k]
    if (Hjk!=1){
      pre<-rho[j,k]*sin(pi/2*Hjk) 
      pim<- -eta[j,k]*cos(pi/2*Hjk)
    }
    else{
      pre<-rho[j,k];pim<- -pi/2*eta[j,k]
    }
    G[j,k]<-complex(real=pre,imaginary=pim)
    G[j,k]<-gamma(H[j]+H[k]+1)*G[j,k]
  }}
  ltmp<-sum(eigen(G)$values<0)
  #print(H)
  #ltmp<-existMFBM(H=H,rho=rho,eta=eta,forceEta=forceEta,choix=choix)  
  txt<-paste("cannot be simulated for these parameters :",ltmp," negative eigenvalues")
  if (print & (ltmp==0)) cat('No negative eigenvalues for G \n')
  if (ltmp>0) stop(txt)
  
  xlogx<-function(h){
    res<-NULL
    for (e in h) { if (e==0) res<-c(res,0) else res<-c(res,h*log(abs(h)))}
    res
  }
  
  gammauv<-function(h,Hj,Hk,rhojk,etajk){
    #pour h>=1
    Hjk<-Hj+Hk
    wjk<-function(h,rhojk,etajk){
      tmp2<-NULL
      for (e in h){
        if (Hjk!=1) {
          tmp<-(rhojk-etajk*sign(e))*abs(e)^Hjk 
        }	else { tmp<-rhojk*abs(e)-etajk*xlogx(e)}
        tmp2<-c(tmp2,tmp)
      }
      tmp2
    }
    .5*(wjk(h-1,rhojk,etajk)-2*wjk(h,rhojk,etajk)+wjk(h+1,rhojk,etajk))
  }
  cuvj <- function(u,v, m){
    ## provides the first line of the matrix C_m(H_u,H_v) de Wood et Chan version p>1 !
    Hu<-H[u]
    Hv<-H[v]
    Huv<-Hu+Hv
    z0<-rho[u,v]
    z1<-gammauv(1:(m/2-1),Hu,Hv,rho[u,v],eta[u,v])
    z2<-gammauv(m/2,Hu,Hv,rho[u,v],eta[u,v])+gammauv(m/2,Hv,Hu,rho[v,u],eta[v,u])
    z2<-z2/2
    z3<-gammauv(m-(m/2+1):(m-1),Hv,Hu,rho[v,u],eta[v,u])	
    c(z0,z1,z2,z3)
  }
  m<-2^(trunc(log(n)/log(2))+2)
  tab<-B<-array(0,dim=c(p,p,m))
  for (u in 1:p){ for (v in 1:p){
    vpCuv <- cuvj( u, v,m)
    vpCuv <- (fft(c(vpCuv), inverse = FALSE))
    vpCuv<-(vpCuv)
    tab[u,v,]<-vpCuv
  }}
  res<-list()
  for (iSamp in (1:nbSamp)){
    # cat('\r i=',iSamp)
    W<-matrix(0,nr=m,ncol=p)
    ## Simulations of U and V and storage in an array
    Z<-matrix(0,nr=m,nc=p)
    Z[1,]<-rnorm(p)/sqrt(m)
    Z[m/2+1,]<-rnorm(p)/sqrt(m)
    
    trace<-NULL
    for (j in (0:(m-1))){
      if ((j>0) & (j<=(m/2-1))){
        U<-rnorm(p); V<-rnorm(p)
        Z[j+1,]<-complex(real=U,imaginary=V) /sqrt(2*m)
        Z[m-j+1,]<-complex(real=U,imaginary=(-V)) /sqrt(2*m)
      } 
      A<-tab[,,j+1]
      tmp<-eigen(A)
      vpA<-Re(tmp$values)
      vecpA<-	tmp$vec
      vecpAinv<-solve((tmp$vec))
      vpAPlus<-apply(cbind(vpA,rep(0,length(vpA))),1,max)
      vpANeg<-apply(cbind(-vpA,rep(0,length(vpA))),1,max)
      trace<-c(trace,sum(vpANeg)/sum(vpA))
      B[,,j+1]<-vecpA %*% diag(sqrt((vpAPlus))) %*% vecpAinv
      W[j+1,]<-B[,,j+1] %*% Z[j+1,]
    }
    X<-Re(mvfft(W,inverse=F)[1:n,])
    vFBM<-apply(X,2,cumsum)
    if (print) cat('  Mean trace:',mean(trace),'\n')
    res[[iSamp]]<-vFBM
  }
  if (nbSamp==1) {  res2<- res[[1]]} else res2<-res	
  res2
  
}

print('...start executing R script')

args = commandArgs(trailingOnly=TRUE)
# print(args)
arg_n <- strtoi(args[1])
arg_h1 <- as.double(args[2]) 
arg_h2 <- as.double(args[3]) 
arg_h3 <- as.double(args[4]) 
arg_h4 <- as.double(args[5]) 
arg_h5 <- as.double(args[6]) 
arg_rho <- as.double(args[7])
arg_csv <- args[8]


H = c()
sig = c()
roh = c()
eta = c()
dim = 0
if (is.na(arg_h1) == FALSE) {
  # H[0] <- arg_h1
  H <- append(H, arg_h1)
  sig <- append(sig, 1)
  dim = dim + 1
}
if (is.na(arg_h2) == FALSE) {
  # H[1] <- arg_h2
  H <- append(H, arg_h2)
  sig <- append(sig, 1)
  dim = dim + 1
}
if (is.na(arg_h3) == FALSE) {
  # H[2] <- arg_h3
  H <- append(H, arg_h3)
  sig <- append(sig, 1)
  dim = dim + 1
}
if (is.na(arg_h4) == FALSE) {
  # H[3] <- arg_h4
  H <- append(H, arg_h4)
  sig <- append(sig, 1)
  dim = dim + 1
}
if (is.na(arg_h5) == FALSE) {
  # H[4] <- arg_h5
  H <- append(H, arg_h5)
  sig <- append(sig, 1)
  dim = dim + 1
}
# print(dim)
for (x in 0:(dim-1)) {
  for (y in 0:(dim-1)) {
    if (x == y) {
      roh <- append(roh, 1)
    } else {
      roh <- append(roh, arg_rho)
    }
  }
}

for (x in 0:(dim-1)) {
  for (y in 0:(dim-1)) {
    if (x == y) {
      eta <- append(eta, 0)
    } else {
      h1 = as.double(args[2 + x])
      h2 = as.double(args[2 + y])
      var = (0.1 / ( 1 - h1 - h2))
      #var = as.double(0)
      if (y > x) {
        var = var * -1
      }
      eta <- append(eta, var)
    }
  }
}

paste(H,collapse=" ")
# paste(sig,collapse=" ")
# paste(roh,collapse=" ")
# paste(eta,collapse=" ")


mfBm <- simMFBM(
  n=arg_n, 
  H=H,
  sig=sig,
  rho=matrix(roh,nc=dim),
  eta=matrix(eta,nc=dim),
  plot=FALSE,
  print=TRUE,
  choix=NULL,
  forceEta=FALSE,
  nbSamp=1
  )

  # mfBm <- simMFBM(
  # n=arg_n, 
  # H=H,
  # sig=sig,
  # rho=matrix(roh,nc=dim),
  # eta=matrix(c(0,0.1/(1-arg_h1-arg_h2),0.1/(1-arg_h1-arg_h3),0.1/(1-arg_h1-arg_h4),0.1/(1-arg_h1-arg_h5),-0.1/(1-arg_h2-arg_h1),0,0.1/(1-arg_h2-arg_h3),0.1/(1-arg_h2-arg_h4),0.1/(1-arg_h2-arg_h5),-0.1/(1-arg_h3-arg_h1),-0.1/(1-arg_h3-arg_h2),0,0.1/(1-arg_h3-arg_h4),0.1/(1-arg_h3-arg_h5),-0.1/(1-arg_h4-arg_h1),-0.1/(1-arg_h4-arg_h2),-0.1/(1-arg_h4-arg_h3),0,0.1/(1-arg_h4-arg_h5),-0.1/(1-arg_h5-arg_h1),-0.1/(1-arg_h5-arg_h2),-0.1/(1-arg_h5-arg_h3),-0.1/(1-arg_h5-arg_h4),0),nc=5),
  # plot=FALSE,
  # print=TRUE,
  # choix=NULL,
  # forceEta=FALSE,
  # nbSamp=1
  # )

# mfBm <- simMFBM(n=arg_n, H=c(.5,.5),
#                    sig=c(1,1),rho=matrix(c(1,.5,.5,1),nc=2),eta=matrix(c(0,.5,-.5,0),nc=2),
#                    plot=FALSE,print=TRUE,choix=NULL,forceEta=FALSE,nbSamp=1)

write.table(mfBm, file=arg_csv, row.names=FALSE, col.names=FALSE, sep=",")

print('...done executing R script')