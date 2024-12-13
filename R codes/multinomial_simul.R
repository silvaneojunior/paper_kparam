##############################################################################
source('R codes/base.r')

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

stan_model=stan_model(file='Stan codes/model.stan')


r=3
k=r-1
T=55
H=0.01
structure <- polynomial_block(p = 1, order = 1, H= H,R1=1) * k

#### If starting a new simulation ####
write("model, y, var, time, sample, n.obs, N, mean, icl, icu",file='data/simul_data_multinomial.csv')
write("model, n.obs, sample, N, time_spent",file='data/time_data_multinomia.csv')

#### If continuing a simulation ####
# 
# init=449 # The last completed set in the last run
# 
# read.csv('data/simul_data_multinomial.csv') %>%
#   filter(sample<init) %>%
#   write.csv(file='data/simul_data_multinomial.csv',row.names=FALSE)
# read.csv('data/time_data_multinomial.csv') %>%
#   filter(sample<init) %>%
#   write.csv(file='data/time_data_multinomial.csv',row.names=FALSE)
####################################

for(rep in init:2000){
  for(N.smp in c(1,2,5,10,20,50,100)){
    # m_t = m_t-1 + omega_t
    if(H>0){
      omega=matrix(rnorm(k*T,0,sqrt(H)),k,T)
    }else{
      omega=matrix(0,k,T)
    }
    mt[,1]=rnorm(k,0,0.1)
    for(i in 1:k){
      mt[i,]=cumsum(omega[i,])
    }
    
    # y|N, p ~ Multinom(N,p)
    # log(p_i/p_k)=mt_i
    y=matrix(NA,r,T)
    for(i in 1:T){
      alpha=exp(c(mt[,i],0))/sum(exp(c(mt[,i],0)))
      y[,i]=rmultinom(1,N.smp,alpha)
    }

    for(t in c(2,5,10,20,50)){
      cat(paste0(rep,' - ',N.smp,' - ',t,'                     \r'))

      #### kDGLM ####
      
      outcome <- Multinom(p = structure$pred.names, data = t(y)[1:t,])

      start.kDGLM=Sys.time()
      fitted.kDGLM <- fit_model(structure, outcome)
      end.kDGLM=Sys.time()
      
      #### STAN ####

      init.func=function(){list(theta_1=fitted.kDGLM$mts[,1])}

      start.stan=Sys.time()
      stan_fit=sampling(stan_model,chains=5,iter=2000,init=init.func,
                        data=list(t=fitted.kDGLM$t,n=fitted.kDGLM$n,k=fitted.kDGLM$k,
                                  y=t(y)[1:t,],
                                  F=fitted.kDGLM$FF |> aperm(c(3,2,1)),
                                  G=fitted.kDGLM$G[,,-1,drop=FALSE] |> aperm(c(3,1,2)),
                                  a1=fitted.kDGLM$a1,R1=fitted.kDGLM$R1,
                                  W=fitted.kDGLM$W[,,t]))
      end.stan=Sys.time()


      stan_smp=extract(stan_fit)
      N=length(stan_smp$theta[,1,1])

      # Testing convergence
      # coda::effectiveSize(stan_smp$theta[,1,])
      # coda::raftery.diag(stan_smp$theta[,1,])
      
      ####################################################
      
      simul.data=data.frame()
      for(var in 1:k){
        simul.data=simul.data %>% rbind(data.frame(model='kdglm',y=mt[var,1:t],var=paste0('theta.',var),
                                                   time=1:t,sample=rep, n.obs=t,
                                                   N=N.smp,
                                                   mean=fitted.kDGLM$mts[var,],
                                                   icl=fitted.kDGLM$mts[var,]-1.96*sqrt(fitted.kDGLM$Cts[var,var,]),
                                                   icu=fitted.kDGLM$mts[var,]+1.96*sqrt(fitted.kDGLM$Cts[var,var,])))
        simul.data=simul.data %>% rbind(data.frame(model='stan',y=mt[var,1:t],var=paste0('theta.',var),
                                                   time=1:t,sample=rep, n.obs=t,
                                                   N=N.smp,
                                                   mean=colMeans(stan_smp$theta[,,var]),
                                                   icl=colQuantile(stan_smp$theta[,,var],0.025),
                                                   icu=colQuantile(stan_smp$theta[,,var],0.975)))
      }

      for(t.adv in 1:5){
        theta.smp.kDGLM=rmvnorm(N,fitted.kDGLM$mts[,t],fitted.kDGLM$Cts[,,t]+t.adv*fitted.kDGLM$W[,,t])
        # Note that the prediction that m_t=m_{t-1}+omega_t
        # As such, given the sample for m_t obtained by STAN, to obtain a sample for m_{t+h}, it is enough to sample though omega_{t+1} to omega_{t+h}
        # Addionally, as we are interested in the marginal distribution of m_t+h and the omega_t's are independent, we can obtain a sample for each $h$ separately, adapting the covariance matrix to reflect the cummulative effect of omega_{t+1} to omega_{t+h}.
        theta.smp.stan=t(stan_smp$theta[,t,])+rmvnorm(N,rep(0,k),t.adv*fitted.kDGLM$W[,,t])
        eta.smp.kDGLM=exp(theta.smp.kDGLM)
        eta.smp.stan=exp(theta.smp.stan)
        y.smp.kDGLM=y.smp.stan=matrix(NA,r,N)
        for(i in 1:N){
          y.smp.kDGLM[,i]=rmultinom(1,sum(y[,t+t.adv]),c(eta.smp.kDGLM[,i],1)/sum(c(eta.smp.kDGLM[,i],1)))
          y.smp.stan[,i]=rmultinom(1,sum(y[,t+t.adv]),c(eta.smp.stan[,i],1)/sum(c(eta.smp.stan[,i],1)))
        }

        for(var in 1:k){
          simul.data=simul.data %>% rbind(data.frame(model='kdglm',y=y[var,t+t.adv],var=paste0('y.',var),
                                                     time=t+t.adv,sample=rep, n.obs=t,
                                                     N=N.smp,
                                                     mean=mean(y.smp.kDGLM[var,]),
                                                     icl=quantile(y.smp.kDGLM[var,],0.025),
                                                     icu=quantile(y.smp.kDGLM[var,],0.975)))
          simul.data=simul.data %>% rbind(data.frame(model='stan',y=y[var,t+t.adv],var=paste0('y.',var),
                                                     time=t+t.adv,sample=rep, n.obs=t,
                                                     N=N.smp,
                                                     mean=mean(y.smp.stan[var,]),
                                                     icl=quantile(y.smp.stan[var,],0.025),
                                                     icu=quantile(y.smp.stan[var,],0.975)))
        }
      }

      # Appending is faster for large simulations
      write.table(simul.data,file='data/simul_data_multinomial.csv',append =TRUE,col.names = FALSE,row.names = FALSE,sep=',',dec='.')


      time.data=data.frame(model=c('kdglm','stan'),
                           n.obs=t,sample=rep,
                           N=N.smp,
                           time_spent=c(difftime(end.kDGLM,start.kDGLM,units='secs'),
                                        difftime(end.stan,start.stan,units='secs')))

      write.table(time.data,file='data/time_data_multinomial.csv',append =TRUE,col.names = FALSE,row.names = FALSE,sep=',',dec='.')
    }
  }
}