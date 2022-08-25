################################################################################
### Estimation Code ############################################################
### 1) Under the known strong hierarchy structure condition ####################
### 2) Use Expert knowledge for prior distribution #############################
### * Original Code (Fully data based) was suggested by Chung ##################
################################################################################


# One-to-one correspondence between class number(0 ~ 2^K-1) <-> attribute pattern
as.binary = function(x){
  ans = NULL
  while(any(x!=0)){
    ans = cbind(x%%2,ans)
    x = floor(x/2)
  }
  ans
}


# Bayesian Estimation - Main function (Under strong hierarchy condition)
# [Add or Modified] 
# | QI      : Q-matrix information matrix suggest by test developers 
# |           Each QI entries mean the probability of 'corresponding q-entry is equal to 1'
# | A       : Adjacent matrix that indicate the hierarchy information 
# | Qnum    : Each q-vectors will be saved or extracted as equivalence class numbers
# | phi     : pmf for each equivalence classes of q-vectors
# |           prior for phi - Dirichlet(d)
# | d       : parameters for Dirichlet prior for phi, determined by QI
# | lambda  : effect size of QI for phi (lambda reflects belief in QI matrix)
# | alpha   : examinees' attribute patterns (Only possible attribute patterns under strong H)
# | pi      : pmf for each possible attribute patterns 
# | Sample.Q  : Function for choosing equivalence class numbers for each q-vectors
# [Others]  
# | Y : Response
# | N : number of examinees, J : number of items, K : number of attributes
# | g & g.start : (Initial) guess parameters
# | s & s.start : (Initial) slip parameters 
# | pi.start : Initial pi (length:2^K)
# | niter : number of iterations

EstQSHMCMC = function(Y, QI=NULL, A=NULL, g.start=NULL, s.start=NULL, pi.start=NULL, lambda=1, niter){
  
  N = nrow(Y)
  J = ncol(Y)
  K = ncol(QI)
  
  Y = as.matrix(Y)
  
  if(is.null(QI) || is.null(A))
    stop('User must supply Q-matrix Info matrix and Hirarchy Adjacent matrix!\n\n')
  
  
  # Possible Attribute Patterns under Strong Hierarchy Structure (Columns of Access Matrix & 0 vector)
  R=diag(K)
  for(k in 1:K){
    R=R+R%*%A
  }
  R=ifelse(R>0,1,0) # Access Matrix
  pos.a = subset(t(R),rowSums(t(R))!=0)
  pos.a = rbind(c(0,0,0),pos.a)
  nclass = nrow(pos.a) # Number of possible attribute patterns (= number of equivalence classess)
  all.a = as.binary(0:(2^K-1)) #2^K by K : all class <->  attribute pattern
  eq.a = ifelse(all.a%*%t(R)>0,1,0) #2^K by K : all class <-> attribute pattern including precedence attribute 
  class.all = ifelse(tcrossprod(eq.a,pos.a)+tcrossprod(1-eq.a,1-pos.a)==K,1,0) #2^K by nclass : Each attribute patterns in which equivalence class 
  
  
  #  Initial pi - If not given : from Dirichlet(1,1,...,1)
  if(is.null(pi.start)){
    pi = rgamma(nclass, rep(1,nclass)) # Use rgamma function for Dirichlet dis.  
    pi = pi/sum(pi)
  }
  else
    pi = pi.start%*%class.all
  
  
  # Initial g, s - If not given : from unif(0.1, 0.3)
  if(is.null(g.start))
    g = runif(J, 0.1, 0.3)   
  else
    g = g.start
  if(is.null(s.start))
    s = runif(J, 0.1, 0.3)   
  else
    s = s.start
  
  
  # Initial Q : determined by QI (cut-off = 0.5)
  Q = ifelse(QI>=0.5, 1, 0) 
  Qnum = as.vector(Q%*%(2^((K-1):0)))+1 #Convert each q-vec to corresponding number 1~2^K 
  Qnum = sapply(1:J,function(j){which(class.all[Qnum[j],]==1)})-1 #Vector of equivalence class numbers of q-vectors 
  
  
  # Initial phi : determined by QI 
  phi=matrix(NA,J,nclass-1)
  for(j in 1:J){
    qi_j=tcrossprod(rep(1,2^K),QI[j,])
    p_j=apply(all.a*qi_j+(1-all.a)*(1-qi_j),1,prod)
    P_j=p_j%*%class.all
    P_j=P_j/(1-P_j[1])
    phi[j,]=P_j[-1]
  }
  d=lambda*phi 
    
  Yt = t(Y)
  alpha.out = NULL; pi.out = pi; g.out = g; s.out = s; Qnum.out = Qnum; phi.out = phi; 
  
  
  for(ii in 1:niter){
    #Update alpha
    etajm = tcrossprod(Q,pos.a) 
    natt = apply(Q,1,sum) 
    etajm = (etajm == natt)
    pp = g*(1-etajm) + (1-s)*etajm 
    ll = Y %*% log(pp) + (1-Y)%*%log(1-pp) 
    ll = sweep(ll,2,log(pi),'+') 
    pp = exp(ll) 
    pp = apply(pp,1,cumsum) 
    pp = sweep(pp,2,pp[nclass,],'/')
    u = runif(N) 
    alphanum = apply(sweep(pp,2,u,'<'),2,sum)+1 
    alpha = pos.a[alphanum,]  
    alpha.out = cbind(alpha.out,alpha)  
    
    #pi Update
    cc = apply(outer(alphanum,1:nclass,'=='),2,sum) 
    pi = rgamma(nclass, 1+cc) 
    pi = pi/sum(pi)
    pi.out = rbind(pi.out,pi) 
    
    #g,s Update 
    etaim = tcrossprod(Q,alpha) 
    etaim = (etaim == natt) 
    ga = apply((1-etaim)*Yt,1,sum) 
    gb = apply((1-etaim)*(1-Yt),1,sum) 
    sa = apply(etaim*(1-Yt),1,sum) 
    sb = apply(etaim*Yt,1,sum)     
    g = qbeta(runif(J, 0,pbeta(1-s,1+ga,1+gb)),1+ga,1+gb) 
    s = qbeta(runif(J,0,pbeta(1-g,1+sa,1+sb)),1+sa,1+sb)  
    g.out = rbind(g.out,g) 
    s.out = rbind(s.out,s) 
    
    #Q Update > Use sample.Q function 
    Qnum = sapply(1:J,function(j){sample.Q(pos.a[-1,],g[j],s[j],Y[,j],alpha,phi[j,])}) 
    Q = t(sapply(1:J,function(j){pos.a[Qnum[j]+1,]}))
    
    #phi Update 
    d1 = d + diag(K)[Qnum,] #parameters for Dirichlet posterior 
    phi = t(sapply(1:J,function(j){rgamma(nclass-1, d1[j,])}))
    phi = phi/rowSums(phi)
    phi.out = cbind(phi.out, phi) 
  }
  alpha.out = array(alpha.out, c(N,K,niter)) 
  phi.out = array(phi.out, c(J,K,niter)) 
  out = list(alpha=alpha.out, pi=pi.out, s=s.out, g=g.out, Qnum=Qnum.out, phi=phi.out)
  class(out) = 'cdmcmc'
  out
}

# Function for updating j-th q-vector (output : equivalence class number)
sample.Q = function(pos.q, g_j, s_j, Y_j, alpha, phi_j){
  natt = apply(pos.q,1,sum) 
  cc = tcrossprod(pos.q,alpha)  
  etaim = (cc==natt) 
  phi_j[phi_j<1e-8] = 1e-8 
  phi_j[phi_j>1-1e-8] = 1-1e-8
  ga = (1-etaim)%*%Y_j 
  gb = (1-etaim)%*%(1-Y_j) 
  sa = etaim%*%(1-Y_j) 
  sb = etaim%*%Y_j     
  pm = ga*log(g_j) + gb*log(1-g_j) + sa*log(s_j) + sb*log(1-s_j)
  pm = pm + log(phi_j)
  pm = pm - max(pm)
  pm = as.vector(exp(pm)) 
  pm = pm/sum(pm) # pmf of posterior distribution 
  qnum = sample(1:nrow(pos.q),size=1,prob=pm) 
}



################################################################################
### Simulation Experiment ######################################################
################################################################################


# Experiment 1 : Hierarchy structure : a1  / Q-matrix : q1 
load("./Data_DINAH(st1.SH.mid.1000).RData")

## Case1 : QI = true q1-mat with 80% confidence 
qi1_1 = ifelse(q1==1, 0.8, 0.2) 
yy=response
system.time(out <- EstQSHMCMC(yy, QI=qi1_1, A=a1, niter=100)) #niter=100000

## Case2 : QI = 90% true q1-mat with 80% confidence (adding att to P7)
qi1_2 = qi1_1
qi1_2[7,] = c(0.8,0.8,0.8) 

## Case3 : QI = 90% true q1-mat with 80% confidence (deleting att to P7)
qi1_3 = qi1_1
qi1_3[7,] = c(0.8,0.4,0.2)


################################################################################

# Experiment 2 : Hierarchy structure : a2  / Q-matrix : q2 
load("./Data_DINAH(st2.SH.mid.1000).RData")

## Case1 : QI = true q2-mat with 80% confidence
qi2_1 = ifelse(q1=2, 0.8, 0.2) 
system.time(out <- EstQSHMCMC(yy, QI=qi2_1, A=a1, niter=100)) #niter=100000

## Case2 : QI = 94% true q2-mat with 80% confidence (adding att to P9)
qi2_2 = qi2_1
qi2_2[9,] = c(0.8,0.8,0.8,0.2)

## Case3 : QI = 94% true q2-mat with 80% confidence (deleting att to P9)
qi2_3 = qi2_1
qi2_3[9,] = c(0.8,0.2,0.2,0.2)












