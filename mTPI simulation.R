
## This is the R code for performing simulation studies based on the mTPI method (Ji et al, 2009).
## The users only need to specify an equivalence interval (pT-epsi1, pT+epsi2)
## where pT is the target tox probability for the MTD and epsi1 and epsi2 are two small values 
## so that any dose in the equivalence interval will be considered potential candidates for the true MTD. 
## The design is then upon the specification of D, the number of candidate doses, sample size, cohort size, starting dose, and simN, the number of simulations.

##library(SAGx)
#system("rm -f simresult.txt")

## The mPTI algorithm and program were developed by Dr. Yuan Ji from MD Anderson Cancer Center in 2009.
## I modified the code and conducted the simulation study for my research project.

dd.dir <- "Users/Megan Wang/Documents/Research/"  # Working directory

setwd(dd.dir) 			# sets working directory

set.seed(6)

############################################### Functions needed ####################################

## pava is the pool adjacent violator algorithm to perform isotonic transformation for the posterior means later
pava <- function (x, wt = rep(1, length(x))) 
{
  n <- length(x)
  if (n <= 1) 
    return(x)
  if (any(is.na(x)) || any(is.na(wt))) {
    stop("Missing values in 'x' or 'wt' not allowed")
  }
  lvlsets <- (1:n)
  repeat {
    viol <- (as.vector(diff(x)) < 0)
    if (!(any(viol))) 
      break
    i <- min((1:(n - 1))[viol])
    lvl1 <- lvlsets[i]
    lvl2 <- lvlsets[i + 1]
    ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
    x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
    lvlsets[ilvl] <- lvl1
  }
  x
}

## betavar computes variances of beta distributions 

betavar<-function(a,b){
  resp <- a*b/((a+b)^2*(a+b+1))
  return(resp)
}

###########################################################################################################################





################## simulation setup ############################################

simN<-100000  ## 100000 simulations

pT<-.3; D<-5; sampsize <- 30; csize <-3; startdose <- 1 ## target tox probility is pT, with D doses, cohorts size is csize, and max sample size is sampsize

eps1<-.05; eps2<-.05; ## (pT-epsi1, pT+epsi2) is the equivalence interval 

a<-1; b<-1; ## the prior is beta(1,1). Keep a, b fixed at 1 -- essential for this design to work.

xi  <- 0.8  ## the cutoff probability for excessive toxicity.

xi2 <- 1.0

write("******************* simulation with the following (a, b, eps1, eps2) and pT values *************", "simresult.txt", append=T)
write(c(a,b, eps1, eps2,pT, sampsize), ncolumns=6, "simresult.txt", append=T) 
write("**********************", "simresult.txt", append=T)




######################################## Set up scenarios ######################
for(sc in c(1:5)){

  if (sc==1){
    p <- c(0.3, 0.46, 0.5, 0.54,0.58)
  }
  if (sc==2){
    p <- c(0.16, 0.30, 0.47, 0.54, 0.60)
  }
  if (sc==3){
    p <- c(0.04, 0.15, 0.30, 0.48, 0.68)
  }
  if (sc==4){
    p <- c(0.02, 0.07, 0.12, 0.3, 0.45)
  }
  if (sc==5){
    p <- c(0.02, 0.06, 0.10, 0.13, 0.3)
  }

  datan<-matrix(rep(0,simN*D),ncol=D)
  datax<-matrix(rep(0,simN*D),ncol=D)
  
  rez<-rep(0,simN)

############# End of scenarios ###############


  ####################################### Start simulations ###############
  for(sim in 1:simN){
    
    x<-rep(0,D); n<-rep(0,D)
    pa<-rep(0, D); pb<-rep(0,D)
    q <- rep(0,3)
    d<-startdose; st<-0; nodose<-0; maxdose<-1; toxdose<-D+1; seldose<-0
    
    
    while(st==0){  ## st = 1 indicates the trial must be terminated
      maxdose<-max(maxdose, d)

      ### generate random toxicity response
      xx <- 0
      for(i in 1:csize){
        ttt <- runif(1)
        if(ttt < p[d]) xx <- xx + 1
      }
      
      x[d] <- x[d] + xx; n[d] <- n[d] + csize
           
      #### Update posterior beta distribution
      pa[d]<-x[d]+a; pb[d]<-n[d]-x[d]+b
      pa[d+1]<-x[d+1]+a; pb[d+1]<-n[d+1]-x[d+1]+b
      
      ###Compute the indicator T_{i+1} to see if the next dose is too toxic
      if(d<D){
        temp<-1-pbeta(pT, pa[d+1], pb[d+1])
        if(temp>xi) {tt<-1; toxdose<-d+1} else {tt<-0}
      }

      
      ##Compute the UPM for three intervals defined by the equivalence interval      
      q[1] <- (1-pbeta(eps2+pT, pa[d], pb[d]))/(1-eps2-pT)
      q[2] <- (pbeta(eps2+pT, pa[d], pb[d]) - pbeta(pT-eps1, pa[d], pb[d]))/(eps2+eps1)
      q[3] <- (pbeta(pT-eps1, pa[d], pb[d])/(pT-eps1))*(1-tt)
      
      	  
      ### implement the dose-assignment rules based on the UPM
      if(d==1){
        ## if the first dose is too toxic, the trial will be terminated
        if((1-pbeta(pT, pa[d], pb[d]))>xi){st=1; nodose<-1;} 
        else{
          if((q[2]>q[1])&&(q[2]>q[3])){d<-d}
          if((q[3]>q[1])&&(q[3]>q[2])){d<-d+1}
        }
      }
      else{
        if(d==D){
          if((q[1]>q[2])&&(q[1]>q[3])){d<-d-1}
          if((q[2]>q[1])&&(q[2]>q[3])){d<-d}
        }
        else{
          if((d>1)&&(d<D)){
            if((q[1]>q[2])&&(q[1]>q[3])){d<-d-1}
            if((q[2]>q[1])&&(q[2]>q[3])){d<-d}
            if((q[3]>q[1])&&(q[3]>q[2])){d<-d+1}
          }
        }
      }
      total<-sum(n)
      if(total >= sampsize){st<-1}
    }
    
    ### compute the posterior mean from the beta distribution
    if(nodose==0){
      tdose<-min(maxdose, toxdose-1)
      pp<-rep(-100,tdose)
      pp.var<-rep(0, tdose)
      
      for(i in 1:tdose){
        pp[i] <- (x[i]+.005)/(n[i]+.01); pp.var[i] <- betavar(x[i]+0.005, n[i]-x[i]+0.005) ### here adding 0.005 is to use beta(0.005, 0.005) for estimating the MTD, which is different from the dose-finding.
      }
     
      pp<-pava(pp, wt=1/pp.var)  ## perform the isotonic transformation using PAVA with weights being posterior variances
      
      for(i in 2:tdose){
        pp[i] <- pp[i] + i*1E-10 ## by adding an increasingly small number to tox prob at higher doses, it will break the ties and make the lower dose level the MTD if the ties were larger than pT or make the higher dose level the MTD if the ties are smaller than pT
      }
      seldose<-sort(abs(pp-pT), index.return=T)$ix[1]
      ##seldose is the final MTD that is selected based on the order-transofromed posterior means
    }
    rez[sim] <- seldose;
    for(i in 1:D){
      datan[sim,i] <- n[i]
      datax[sim,i] <- x[i]
    }
  }
  
  ##rez[is.na(rez)]<-0
  aaa<-rep(0,D)


  ################## output results ################################
  for(i in 1:D){
    aaa[i] <- sum(rez==i)/simN ### aaa is the propotion of be selected as the MTD
  }
  write(sc, "simresult.txt", append=T) ## Scenario number
  write(p, "simresult.txt", ncolumns = 8, append=T) ## the true tox probabilities
  write(aaa, ncolumns = 8, "simresult.txt", append=T) ## aaa is the vector containing the percentages of selection
  write(apply(datan,2,mean), ncolumns = 8, "simresult.txt", append=T) ## prints out the average number of patients at each dose
  write(sum(apply(datan,2,mean)),  "simresult.txt", append=T) ## prints out the average total sample size
  write(sum(datax)/sum(datan), "simresult.txt", append=T)  ## prints ot the overall toxicity
  write("**********************", "simresult.txt", append=T)
}


