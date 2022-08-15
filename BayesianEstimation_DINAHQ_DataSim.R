rm(list=ls())
###############################
#####   Data Simulation   #####
###############################

# DINAH데이터 생성함수 
### Q=Q행렬, N=피험자 수 
### 인접행렬 A는 요소 순서를 잘 설정하여 상부삼각행렬이 되도록 한다.
### Prob벡터를 설정하면 조건부숙달확률, 숙달확률을 조절할 수 있다. (default 0.8/0.7)
### 문항모수 s,g값은 단일 값을 넣어도 되고, vec를 넣어도 된다.
### 강한위계가 default이며 약한위계로 설정하기 위해서는 SH=FALSE로 설정한다.  

DINAH.SIM = function (Q, N, A, Prob=NULL, s=0.2, g=0.2, SH=TRUE){
  J=nrow(Q)
  K=ncol(Q)
  
  #인접행렬 검증
  S=0
  for(k in 2:K){
    for(kk in 1:(k-1)){S=S+A[k,kk]}
  }
  if(S>0)
    stop('Adjacency Matrix A should be an upper triangular matrix! \n\n')
  
  #근원요소집합
  root_att=NULL
  for(k in 1:K){
    if(sum(A[,k])==0){root_att=c(root_att,k)} #근원숙달확률=0.7
  }

  #숙달확률 및 조건부숙달확률 설정이 따로 없을 경우   
  if(is.null(Prob)){
    Prob=rep(0.8,K) #조건부숙달확률=0.8 
    for(k in root_att){Prob[k]=0.7} #근원숙달확률=0.7
  }
  
  nonmas_Prob=rep(0,K)
  #약한위계로 설정하였을 경우 (조건만족 안해도 0~0.1확률로 숙달)
  if(!SH){nonmas_Prob=runif(K,0,0.1)}
  for(k in root_att){nonmas_Prob[k]=NA} #근원요소이면 NA값 가지도록 
  

  #피험자 인지요소패턴 
  alpha = matrix(NA, N, K)
  for(k in 1:K){
    #근원요소일 경우 
    if(k %in% root_att){
      alpha[,k]=rbinom(N,1,Prob[k])
    }
    #선행요소가 있을 경우
    else{
      pre_att=which(A[,k]!=0)
      pre_att_mas=rep(1,N)
      for(pa in pre_att){
        pre_att_mas=pre_att_mas*alpha[,pa]
      }
      if(SH){pre_att_prob=pre_att_mas*Prob[k]} #강한 위계 
      else if(!SH){pre_att_prob=ifelse(pre_att_mas==1,Prob[k],nonmas_Prob[k])} #약한 위계 
      comp=c(runif(N, 0, 1))
      alpha[,k]=ifelse(comp<pre_att_prob,1,0)
    }
  }
    
  #피험자 기대응답 
  tm = alpha%*%t(Q) 
  at = matrix(rep(apply(Q, 1, sum), N), nrow=N, byrow=T)
  eta = ifelse(tm==at, 1, 0)   
  
  #피험자 실제응답
  if (length(s)==1){s=rep(s,J)}
  if (length(g)==1){g=rep(g,J)}
  y=eta
  for(j in 1:J){y[,j]=ifelse(eta[,j],1-s[j],g[j])}
  comp=c(runif(N*J, 0, 1))
  y=ifelse(y>comp, 1, 0)
  
  #최종 output
  out = list(alpha=alpha, eta=eta, y=y, slip=s, guess=g, Prob=Prob, nProb=nonmas_Prob)
  out
}

############################
#####   Simulation I   #####
############################
q1 = matrix(c(1,0,0,
              0,1,0,
              0,0,1,
              1,0,0,
              0,1,0,
              0,0,1,
              1,1,0,
              1,0,1,
              0,1,1,
              1,1,1), byrow=TRUE, nrow=10, ncol=3)

a1 = matrix(c(0,1,0,
              0,0,1,
              0,0,0), byrow=TRUE, nrow=3, ncol=3)

q2 = matrix(c(1,0,0,0,
             0,1,0,0,
             0,0,1,0,
             0,0,0,1,
             1,0,0,0,
             0,1,0,0,
             0,0,1,0,
             0,0,0,1,
             1,1,0,0,
             1,0,1,0,
             1,0,0,1,
             0,1,1,0,
             0,1,0,1,
             0,0,1,1,
             1,1,1,0,
             1,1,0,1,
             1,0,1,1,
             0,1,1,1), byrow=TRUE, nrow=18, ncol=4)

a2 = matrix(c(0,1,1,0,
              0,0,0,1,
              0,0,0,1,
              0,0,0,0), byrow=TRUE, nrow=4, ncol=4)


simdata = DINAH.SIM(Q=q2, N=1000, A=a2,SH=TRUE)
alpha=simdata$alpha
eta=simdata$eta
response=simdata$y
#slip=simdata$slip 
#guess=simdata$guess
#Prob=simdata$Prob
#nProb=simdata$nProb
save.image("Data_DINAH(st2.SH.mid.1000).RData")

#SIM.st1.Hard.WH.500 = DINAH.SIM(Q=q1, N=500, A=a1, s=0.3, g=0.2, SH=FALSE)
