rm(list=ls())
###############################
#####   Data Simulation   #####
###############################

# DINAH������ �����Լ� 
### Q=Q���, N=������ �� 
### ������� A�� ��� ������ �� �����Ͽ� ��λﰢ����� �ǵ��� �Ѵ�.
### Prob���͸� �����ϸ� ���Ǻμ���Ȯ��, ����Ȯ���� ������ �� �ִ�. (default 0.8/0.7)
### ���׸�� s,g���� ���� ���� �־ �ǰ�, vec�� �־ �ȴ�.
### �������谡 default�̸� ��������� �����ϱ� ���ؼ��� SH=FALSE�� �����Ѵ�.  

DINAH.SIM = function (Q, N, A, Prob=NULL, s=0.2, g=0.2, SH=TRUE){
  J=nrow(Q)
  K=ncol(Q)
  
  #������� ����
  S=0
  for(k in 2:K){
    for(kk in 1:(k-1)){S=S+A[k,kk]}
  }
  if(S>0)
    stop('Adjacency Matrix A should be an upper triangular matrix! \n\n')
  
  #�ٿ��������
  root_att=NULL
  for(k in 1:K){
    if(sum(A[,k])==0){root_att=c(root_att,k)} #�ٿ�����Ȯ��=0.7
  }

  #����Ȯ�� �� ���Ǻμ���Ȯ�� ������ ���� ���� ���   
  if(is.null(Prob)){
    Prob=rep(0.8,K) #���Ǻμ���Ȯ��=0.8 
    for(k in root_att){Prob[k]=0.7} #�ٿ�����Ȯ��=0.7
  }
  
  nonmas_Prob=rep(0,K)
  #��������� �����Ͽ��� ��� (���Ǹ��� ���ص� 0~0.1Ȯ���� ����)
  if(!SH){nonmas_Prob=runif(K,0,0.1)}
  for(k in root_att){nonmas_Prob[k]=NA} #�ٿ�����̸� NA�� �������� 
  

  #������ ����������� 
  alpha = matrix(NA, N, K)
  for(k in 1:K){
    #�ٿ������ ��� 
    if(k %in% root_att){
      alpha[,k]=rbinom(N,1,Prob[k])
    }
    #�����Ұ� ���� ���
    else{
      pre_att=which(A[,k]!=0)
      pre_att_mas=rep(1,N)
      for(pa in pre_att){
        pre_att_mas=pre_att_mas*alpha[,pa]
      }
      if(SH){pre_att_prob=pre_att_mas*Prob[k]} #���� ���� 
      else if(!SH){pre_att_prob=ifelse(pre_att_mas==1,Prob[k],nonmas_Prob[k])} #���� ���� 
      comp=c(runif(N, 0, 1))
      alpha[,k]=ifelse(comp<pre_att_prob,1,0)
    }
  }
    
  #������ ������� 
  tm = alpha%*%t(Q) 
  at = matrix(rep(apply(Q, 1, sum), N), nrow=N, byrow=T)
  eta = ifelse(tm==at, 1, 0)   
  
  #������ ��������
  if (length(s)==1){s=rep(s,J)}
  if (length(g)==1){g=rep(g,J)}
  y=eta
  for(j in 1:J){y[,j]=ifelse(eta[,j],1-s[j],g[j])}
  comp=c(runif(N*J, 0, 1))
  y=ifelse(y>comp, 1, 0)
  
  #���� output
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