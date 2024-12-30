source('R codes/base.R')


# Creating dataset
set.seed(1635135)
T = 365
time.cutoff=20

phi=0.95
W=0.05

ht=rep(NA,T)
ht_i=0
for(i in 1:T){
  ht_i=phi*ht_i+rnorm(1,0,sqrt(W))
  ht[i]=ht_i
}

sigma2=exp(-ht+mean(ht))
data_ori=rnorm(T,0,sqrt(sigma2))
data_sqrd=data_ori**2

#### Prior for the AR coefficient #### 
# The mean is centered in or best guess. We recommend the use of an informative prior, if possible, otherwise the model may take too much time to get close to the true value.
# For the paper we chose 0.8 as our best guess, avoiding a value too close to the true, but keeping most of the probability mass within a reasonable range.
# As described below, our prior will be such that  P(0.6<phi<1)=0.95, which seemed like a reasonable vague guess.
# The variance is chosen in such a way that the upper bound for the 95% credibility interval is 1.
# The prior for the AR coefficient is not too important, but one should be careful to not choose prior that is too vague.
# We recomend using a prior such that at least 95% of the mass in within the statiory region of the AR(1), i.e., P(-1<phi<1)>=0.95
# m0=0.8
# C0=((1-m0)*0.5)**2
m0=0.95
C0=0

level=polynomial_block(mu=1,R1=1)
volat=TF_block(tau=1,order=1,
               noise.var=W, # To facilitate the comparison with STAN, we consider the variance of the evolutional noise as known.
               a1.coef=m0,R1.coef=C0,
               a1=0,R1=1)
time.norm=Sys.time()
model_norm=fit_model(level,volat,Normal(mu='mu',Tau='tau',data = data_ori))
time.norm=Sys.time()-time.norm

time.gamma=Sys.time()
model_gamma=fit_model(level,volat,Gamma(phi=1/2,mu='tau',data = data_sqrd))
time.gamma=Sys.time()-time.gamma

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

if(FALSE){
  # Fitting the model with STAN
  model=stan_model(file='Stan codes/normal_gamma_model_sim.stan')
  
  time.stan=Sys.time()
  stan_fit=sampling(
    model,
    data=list(n=T,y=data_ori,W=W,m0=m0,C0=C0)
  )
  time.stan=Sys.time()-time.stan
  
  stan_smp=rstan::extract(stan_fit)
  # Saving the model
  save(stan_smp,file='data/volatility_stan_fit_simul.bkp')
}else{
  # Loading the model
  load(file='data/volatility_stan_fit_simul.bkp')
}


############################################################

plot.data=rbind(
  data.frame(
    Time=1:T,
    y=-model_norm$mts[2,],
    ymin=qnorm(0.025,-model_norm$mts[2,],sqrt(model_norm$Cts[2,2,])),
    ymax=qnorm(0.975,-model_norm$mts[2,],sqrt(model_norm$Cts[2,2,])),
    Model='Normal',
    plot='Volatility',
    true_val=-ht
  )
)


plot.data=rbind(plot.data,
                data.frame(
                  Time=1:T,
                  y=model_gamma$mts[2,],
                  ymin=qnorm(0.025,model_gamma$mts[2,],sqrt(model_gamma$Cts[2,2,])),
                  ymax=qnorm(0.975,model_gamma$mts[2,],sqrt(model_gamma$Cts[2,2,])),
                  Model='Gamma',
                  plot='Volatility',
                  true_val=-ht
                )
)

###########################################################

# The kDGLM package provides metrics for model comparison
# But, as the models have different scales (one model the squared log returns and the other the log returns), a correction of scale usting the jacobian is necessary
pred.norm=
  log.like.norm=rep(NA,T)
pred.gamma=coef(model_gamma,eval.pred=TRUE)
pred.gamma=pred.gamma$data$Prediction
log.like.gamma=coef(model_gamma,eval.pred=TRUE,eval.metric=TRUE,
                    lag=1,t.eval = (time.cutoff+1):T)
log.like.gamma=log.like.gamma$metrics$log.like

smp1=simulate(model_norm,50000,lag = 1)

for(i in 1:T){
  mu=smp1$param$Series.1[1,i,]
  sigma2=1/smp1$param$Series.1[2,i,]
  y.log.like=dgamma((data_ori[i]-mu)**2,0.5,0.5/sigma2,log=TRUE)
  # OR
  # y.log.like=dnorm(data_ori[i],mu,sqrt(sigma2),log=TRUE)-0.5*log((data_ori[i]-mu)**2)
  max.log.like=max(y.log.like)
  log.like.norm[i]=log(mean(exp(y.log.like-max.log.like)))+max.log.like
  
  y=rnorm(500000,mu,sqrt(sigma2))
  pred.norm[i]=mean((y-mu)**2)
}

sum(log.like.norm[(time.cutoff+1):T])
sum(log.like.gamma)

###########################################################

pdf(file = 'Figures/normal_gamma_simul_2.pdf',
    width = base.size/dpi, height = base.size/(dpi*2.5),
    family=family_font)

m <- matrix(c(1,2),nrow = 1,ncol = 2,byrow = TRUE)
graphics::layout(mat = m)


par(mgp=c(1.5,0.5,0),mar = c(3, 2.5, 1, 2.5))
ref.data.prop1=plot.data%>%filter(Time>time.cutoff & Model=="Gamma" & plot=="Volatility")
ref.data.prop2=plot.data%>%filter(Time>time.cutoff & Model=="Normal" & plot=="Volatility")
plot(-ht,type='p',pch=16,cex=0.5,xlim=c(time.cutoff,T),
     xlab='Time',ylab='Volatility',xaxs='i',
     ylim=c(-2.5,2.5),
     # main='Volatility',
     family=family_font)
mtext('(A)',side=3,at=-15)
base_ribbon2(
  x=ref.data.prop2$Time,
  ymin=ref.data.prop2$ymin,
  ymax=ref.data.prop2$ymax,
  col='#aaaaaa'
)
lines(ref.data.prop2$Time,ref.data.prop2$y,col='#aaaaaa')

base_ribbon2(
  x=ref.data.prop1$Time,
  ymin=ref.data.prop1$ymin,
  ymax=ref.data.prop1$ymax,
  col='#222222'
)
lines(ref.data.prop1$Time,ref.data.prop1$y,col='#222222')

par(mgp=c(1.5,0.5,0),mar = c(3, 2.5, 1, 2.5))
ref.data.prop1=plot.data%>%filter(Time>time.cutoff & Model=="Gamma" & plot=="Volatility")
ref.data.prop2=plot.data%>%filter(Time>time.cutoff & Model=="Normal" & plot=="Volatility")
plot(-ht,type='p',pch=16,cex=0.5,xlim=c(time.cutoff,T),
     xlab='Time',ylab='Volatility',xaxs='i',
     ylim=c(-2.5,2.5),
     # main='Volatility',
     family=family_font)
mtext('(B)',side=3,at=-15)
base_ribbon2(
  x=ref.data.prop2$Time,
  ymin=ref.data.prop2$ymin,
  ymax=ref.data.prop2$ymax,
  col='#aaaaaa'
)
lines(ref.data.prop2$Time,ref.data.prop2$y,col='#aaaaaa')

base_ribbon2(
  x=(time.cutoff+1):T,
  ymin=colQuantile(-stan_smp$theta,0.025)[(time.cutoff+1):T],
  ymax=colQuantile(-stan_smp$theta,0.975)[(time.cutoff+1):T],
  col='#222222'
)
lines((time.cutoff+1):T,colMeans(-stan_smp$theta,0.025)[(time.cutoff+1):T],col='#222222')

dev.off()
par(mfrow=c(1,1))

# Mean and quantiles of the AR coefficient
c(stan_smp$rho %>% mean, stan_smp$rho %>% quantile(c(0.025,0.975)))%>% round(3) %>% paste(collapse = ' & ')
(model_norm$mts[3,model_norm$t] + c(0,-1,1)*sqrt(model_norm$Cts[3,3,model_norm$t])) %>% round(3) %>% paste(collapse = ' & ')
(model_gamma$mts[3,model_gamma$t] + c(0,-1,1)*sqrt(model_gamma$Cts[3,3,model_gamma$t])) %>% round(3) %>% paste(collapse = ' & ')

