################################################################################
##### Estimation Code ##########################################################
##### This is the estimation code with hirarchy infromation by Heesang Ann #####
##### Original Code (no hirarchy) was suggested by Chung #######################
################################################################################


#1-1 Correspondence between attribute class <-> class number 0 ~ 2^K-1
as.binary = function(x){
  ans = NULL
  while(any(x!=0)){
    ans = cbind(x%%2,ans)
    x = floor(x/2)
  }
  ans
}


#Bayesian Estimation - Main function

EstQSHMCMC = function(Y, QI=NULL, A=NULL, g.start=NULL, s.start=NULL, pi.start=NULL, lambda=1, niter){
  
  N = nrow(Y)
  J = ncol(Y)
  K = ncol(QI)
  
  Y = as.matrix(Y)
  
  if(is.null(QI) || is.null(A))
    stop('User must supply Q-matrix Info matrix and Hirarchy Adjacent matrix!\n\n')
  
  
  # 가능인지패턴 (= 접근행렬의 열 & 0벡터) 
  R=diag(K)
  S=diag(K)
  for(k in 1:K){
    S=S%*%A
    R=R+S
  }
  R=ifelse(R>0,1,0)
  pos.a = subset(t(R),rowSums(t(R))!=0)
  pos.a = rbind(c(0,0,0),pos.a)
  nclass = nrow(pos.a)
  
  
  # 시작 pi 없다면 dir(1,1,...,1)에서 임의 추출
  if(is.null(pi.start)){
    pi = rgamma(nclass, rep(1,nclass)) #dirichlet분포 함수 대신 rgamma로 구함 
    pi = pi/sum(pi)
  }
  else
    pi = pi.start
  
  
  # 시작 g, s 없다면 문항별 0.1~0.3 unif 랜덤추출 
  if(is.null(g.start))
    g = runif(J, 0.1, 0.3)   
  else
    g = g.start
  if(is.null(s.start))
    s = runif(J, 0.1, 0.3)   
  else
    s = s.start
  
  
  # QI로부터 시작 Q 설정 Qnum으로 저장 
  Q = ifelse(QI>=0.5, 1, 0) 
  Qnum = as.vector(Q%*%(2^((K-1):0)))+1 #0~2^K-1 숫자로 각 q-vec 변환 
  all.a = as.binary(0:(2^K-1)) #2^K by K 모든 인지요소 조합 
  eq.a = ifelse(all.a%*%t(R)>0,1,0) #2^K by K 각 인지요소 조합과 같은 동치류인 가능인지벡터 
  class.all = ifelse(tcrossprod(eq.a,pos.a)+tcrossprod(1-eq.a,1-pos.a)==K,1,0) #2^K by nclass 각 조합이 몇 번쨰 동치류에 해당하는지 나타내는 행렬 
  Qnum = sapply(1:J,function(j){which(class.all[Qnum[j],]==1)})-1 #각 문항 q-vec를 동치류로 표현한 벡터 
  
  
  # QI로부터 각 q벡터 시작 분포 확률 phi와 q벡터 prior인 dirichlet 분포의 parameter 생성  
  phi=matrix(NA,J,nclass-1) #J by nclass 각 문항 q-vec 동치류 예상확률 > 초기 확률로 사용 
  for(j in 1:J){
    qi_j=tcrossprod(rep(1,2^K),QI[j,])
    p_j=apply(all.a*qi_j+(1-all.a)*(1-qi_j),1,prod)
    P_j=p_j%*%class.all
    P_j=P_j/(1-P_j[1])
    phi[j,]=P_j[-1]
  }
  qdis.par=lambda*phi
    
  Yt = t(Y)
  alpha.out = NULL; pi.out = pi; g.out = g; s.out = s; Qnum.out = Qnum; phi.out = phi; 
  
  
  for(ii in 1:niter){
    #alpha Update (초기설정포함)
    etajm = tcrossprod(Q,pos.a) #J by nclass 각 인지요소 조합이 해당 문항과 겹치는 인지요소 수 
    natt = apply(Q,1,sum) #현재 Q에서 문항별 필요한 요소 개수
    etajm = (etajm == natt) #J by nclass 각인지요소 조합이 맞출 수 있으면 true, o.w. False      
    pp = g*(1-etajm) + (1-s)*etajm #J by nclass 각 인지요소 조합 별 문항 맞출 확률 
    ll = Y %*% log(pp) + (1-Y)%*%log(1-pp) #N by nclass 각 인지요소 조합 별 log응답확률
    ll = sweep(ll,2,log(pi),'+') #N by nclass 각 인지요소 조합 별 log posterior (loglikelihood+logPrior)
    pp = exp(ll) #N by nclass 각 인지요소 조합별 posterior (비례)
    pp = apply(pp,1,cumsum) #nclass by N 각 인지요소 조합 별 posterior 누적합 
    pp = sweep(pp,2,pp[nclass,],'/') #각 인지요소 조합 별 posterior 누적합을 마지막 항으로 나눠줌 (0~1구간분할)
    #각 사람마다 posterior을 pmf로 하여인지요소조합번호 추출 > 인지요소패턴을 얻는 과정 
    u = runif(N) #N개의 0~1사이 random num
    alphanum = apply(sweep(pp,2,u,'<'),2,sum)+1 #누적합상 u값을 넘지 못하는 row index == u값이 들어가게 되는 구간의 인지요소조합번호  
    alpha = pos.a[alphanum,]    # N by K 행렬 가장 likelihood가 높은 인지요소 조합 
    alpha.out = cbind(alpha.out,alpha) #iter 누적
    
    #pi Update
    cc = apply(outer(alphanum,1:nclass,'=='),2,sum) #각 class별 인원수 벡터 
    pi = rgamma(nclass, 1+cc) #dirichlet분포 함수 대신 rgamma로 구함 
    pi = pi/sum(pi) #각 class별 분포 사후확률  
    pi.out = rbind(pi.out,pi) #iter 누적 matrix  
    
    #g,s Update 
    etaim = tcrossprod(Q,alpha) # J by N 각 피험자가 해당 인지요소와 겹치는 것의 수 
    etaim = (etaim == natt) # J by N 각 피험자가 문항 맞출 수 있으면 True, o.w. False 
    ga = apply((1-etaim)*Yt,1,sum) # 각 문항별 guess상황인 피험자 수 벡터 
    gb = apply((1-etaim)*(1-Yt),1,sum) #각 문항 별 기대,실제 모두 0상황인 피험자 수 벡터
    sa = apply(etaim*(1-Yt),1,sum) #각 문항별 slip상황인 피험자 수 벡터 
    sb = apply(etaim*Yt,1,sum)     #각 문항별 기대, 실제 모두 1상황인 피험자 수 벡터
    g = qbeta(runif(J, 0,pbeta(1-s,1+ga,1+gb)),1+ga,1+gb) #guess update
    s = qbeta(runif(J,0,pbeta(1-g,1+sa,1+sb)),1+sa,1+sb)  #slip update 
    g.out = rbind(g.out,g) #iter 누적 matrix 저장
    s.out = rbind(s.out,s) #iter 누적 matrix 저장
    
    #Q Update > 한 문항씩 순서대로 q-vec num추출하여 붙히는 방식 
    Qnum = sapply(1:J,function(j){sample.Q(pos.a[-1,],g[j],s[j],Y[,j],alpha,phi[j,])}) #Sample.Q로부터 문항별 q벡터num 받아서 vec로 합침
    Qnum.out = rbind(Qnum.out,Qnum) #iter 누적 저장
    Q = t(sapply(1:J,function(j){pos.a[Qnum[j]+1,]}))
    
    #phi Update 
    par = qdis.par + diag(K)[Qnum,]
    phi = t(sapply(1:J,function(j){rgamma(nclass-1, par[j,])}))
    phi = phi/rowSums(phi)
    phi.out = cbind(phi.out, phi) #iter 누적 저장
  }
  alpha.out = array(alpha.out, c(N,K,niter)) #alpha matrix 누적 > 3차원으로  잘라 저장
  phi.out = array(phi.out, c(J,K,niter)) #p matrix 누적 > 3차원으로  잘라 저장
  out = list(alpha=alpha.out, pi=pi.out, s=s.out, g=g.out, Qnum=Qnum.out, phi=phi.out)
  class(out) = 'cdmcmc'
  out
}

#각 문항별 q-vec를 Update하는 함수 (output : 동치류의 번호)
sample.Q = function(pos.q, g_j, s_j, Y_j, alpha, phi_j){
  natt = apply(pos.q,1,sum) #각 q-vec 동치류 후보 별 필요 인지요소 개수 벡터
  cc = tcrossprod(pos.q,alpha)  # nclass-1 by N 각 피험자와 q-vec후보가 겹치는 인지요소 수 
  etaim = (cc==natt) #(nclass-1 by N) 행렬 : 각 피험자별 각 q-vec후보를 포함하면 true o.w. False 
  phi_j[phi_j<1e-8] = 1e-8 
  phi_j[phi_j>1-1e-8] = 1-1e-8
  ga = (1-etaim)%*%Y_j #문항 j 각 q-vec에 대해 guess상황 명수 
  gb = (1-etaim)%*%(1-Y_j) #문항 j 각 q-vec에 대해 기대,응답0상황 명수
  sa = etaim%*%(1-Y_j) #문항j 각 q-vec에 대해 slip상황 명수
  sb = etaim%*%Y_j     #문항j 각 q-vec에 대해 기대, 응답1 상황 명수 
  pm = ga*log(g_j) + gb*log(1-g_j) + sa*log(s_j) + sb*log(1-s_j) #문항 j 각 q-vec에 대해 loglikelihood
  pm = pm + log(phi_j) #문항 j 각 q-vec 동치류 대해 log 사후 확률 
  pm = pm - max(pm)
  pm = as.vector(exp(pm)) 
  pm = pm/sum(pm) #문항 j 각 q-vec에 대해 사후확률 
  qnum = sample(1:nrow(pos.q),size=1,prob=pm) #사후확률 pmf로 하여 q-vec 추출 (조합 번호로)
}


 
#QI matrix : true q1-mat with 80% confidence
qi1_1 = ifelse(q1==1, 0.8, 0.2) 

#QI matrix : 90% true q1-mat with 80% confidence (adding att to P7)
qi1_2 = qi1_1
qi1_2[7,] = c(0.8,0.8,0.8) 

#QI matrix : 90% true q1-mat with 80% confidence (deleting att to P7)
qi1_3 = qi1_1
qi1_3[7,] = c(0.8,0.4,0.2)



#Estimation 메인함수 시행 
yy=response
system.time(out <- EstQSHMCMC(yy, QI=qi1_3, A=a1, niter=100)) #niter=100000













#QI matrix : true q2-mat with 80% confidence
qi2_1 = ifelse(q1=2, 0.8, 0.2) 

#QI matrix : 94% true q2-mat with 80% confidence (adding att to P9)
qi2_2 = qi2_1
qi2_2[9,] = c(0.8,0.8,0.8,0.2)

#QI matrix : 94% true q2-mat with 80% confidence (deleting att to P9)
qi2_3 = qi2_1
qi2_3[9,] = c(0.8,0.2,0.2,0.2)
