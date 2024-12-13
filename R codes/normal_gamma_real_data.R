source('R codes/base.R')

# Dados IBM
data=read.csv('data/m-ibmln.txt',header=FALSE)

T = dim(data)[1]
cutoff=432
media=mean(data[[1]])
desv.pad=sd(data[[1]])
# Change of location to avoid problems with observed 0.
# (optional) Change in scale to avoid possible numerical problems. If you wish to remove the change of scale, remember to adapt the prior to reflect this change.
y <- (data[[1]]-media)/desv.pad
date_label <- seq(as.Date("1926/10/1"), by = "month", length.out = T)
data_ori=y
data_sqr=(data_ori-mean(data_ori))**2

# A sensibility analysis was performed to select this value.
W1=0.05

#### Prior for the AR coefficient #### 
# The mean is centered in or best guess. We recommend the use of an informative prior, if possible, otherwise the model may take too much time to get close to the true value.
# For the paper we chose 0.8 as our best guess, avoiding a value too close to the true, but keeping most of the probability mass within a reasonable range.
# As described below, our prior will be such that  P(0.6<phi<1)=0.95, which seemed like a reasonable vague guess.
# The variance is chosen in such a way that the upper bound for the 95% credibility interval is 1.
# The prior for the AR coefficient is not too important, but one should be careful to not choose prior that is too vague.
# We recomend using a prior such that at least 95% of the mass in within the statiory region of the AR(1), i.e., P(-1<phi<1)>=0.95
m0=0.8
C0=((1-m0)*0.5)**2

level=polynomial_block(mu=1,R1=1)
volat=TF_block(tau=1, order=1,
               noise.var=W1,
               a1.coef=m0, R1.coef=C0,
               a1=0, R1=1)

model_norm=fit_model(level,volat,Normal(mu='mu',Tau='tau',data = data_ori))
model_gamma=fit_model(level,volat,Gamma(phi=1/2,mu='tau',data = data_sqr))

###########################################################

data_norm=coef(model_norm,lag=1,eval.pred = TRUE)$data %>% 
  mutate(Serie='(B)',
         Observation=Observation*desv.pad+media,
         Prediction=Prediction*desv.pad+media,
         C.I.lower=C.I.lower*desv.pad+media,
         C.I.upper=C.I.upper*desv.pad+media,
         max=0.5,min=-0.5)
data_gamma=coef(model_gamma,lag=1,eval.pred = TRUE)$data %>% 
  mutate(Serie='(A)',
         Observation=Observation*desv.pad**2,
         Prediction=Prediction*desv.pad**2,
         C.I.lower=C.I.lower*desv.pad**2,
         C.I.upper=C.I.upper*desv.pad**2,
         max=0.15,min=0)

plot.data=rbind(data_norm,data_gamma) %>% 
  mutate(Time=as.Date(date_label[Time])) 

ggplot(plot.data %>% filter(format(Time,format='%Y')>1960))+
  geom_point(aes(x=Time,y=Observation),alpha=0.3,shape=16)+
  geom_point(aes(x=Time,y=min),alpha=0)+
  geom_point(aes(x=Time,y=max),alpha=0)+
  geom_line(aes(x=Time,y=Prediction))+
  geom_line(aes(x=Time,y=C.I.lower),linetype='dashed')+
  geom_line(aes(x=Time,y=C.I.upper),linetype='dashed')+
  facet_wrap(~Serie,scale='free')+
  coord_cartesian(expand=0)+
  ylab('')+
  scale_x_date('Year',expand=c(0,0))+
  theme_bw()+
  theme(text=element_text(size=font_size,family=family_font),
        legend.text=element_text(size=font_size,family=family_font),
        legend.position = 'bottom',panel.grid=element_blank(),
        strip.background = element_blank(),
        strip.clip = 'off',
        strip.text = element_text(hjust = -0.15,vjust=0.75))

ggsave(
  'Figures/normal_gamma_real_data_1.pdf',
  device='pdf',
  units='px',
  width=base.size,
  height=base.size/2.5
)

###########################################################

library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

if(FALSE){
  stan_fit=stan(
    file='Stan codes/normal_gamma_model_sim.stan',
    data=list(n=T,y=data_ori,W=W1,m0=m0,C0=C0)
  )
  
  stan_smp=rstan::extract(stan_fit)
  save(stan_smp,file='data/volatility_stan_fit_ibm.bkp')
}else{
  load(file='data/volatility_stan_fit_ibm.bkp')
}



plot.data=rbind(
  data.frame(
    Time=1:T,
    y=-model_norm$mts[2,],
    ymin=qnorm(0.025,-model_norm$mts[2,],sqrt(model_norm$Cts[2,2,])),
    ymax=qnorm(0.975,-model_norm$mts[2,],sqrt(model_norm$Cts[2,2,])),
    Model='Normal',
    plot='Volatility'
  )
)


plot.data=rbind(plot.data,
                data.frame(
                  Time=1:T,
                  y=model_gamma$mts[2,],
                  ymin=qnorm(0.025,model_gamma$mts[2,],sqrt(model_gamma$Cts[2,2,])),
                  ymax=qnorm(0.975,model_gamma$mts[2,],sqrt(model_gamma$Cts[2,2,])),
                  Model='Gamma',
                  plot='Volatility'
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
                    lag=1,t.eval = (cutoff+1):T)
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

sum(log.like.norm[(cutoff+1):T])
sum(log.like.gamma)

###########################################################


pdf(file = 'Figures/normal_gamma_real_data_2.pdf',
    width = base.size/dpi, height = base.size/(dpi*2.5),
    family=family_font)

m <- matrix(c(1,2),nrow = 1,ncol = 2,byrow = TRUE)
graphics::layout(mat = m)


par(mgp=c(1.5,0.5,0),mar = c(3, 2.5, 1, 2.5))
ref.data.prop1=plot.data%>%filter(Time>cutoff & Model=="Gamma" & plot=="Volatility")
ref.data.prop2=plot.data%>%filter(Time>cutoff & Model=="Normal" & plot=="Volatility")
plot(date_label[ref.data.prop2$Time] %>% as.Date,rep(0,T-cutoff),type='p',pch=16,cex=0.5,
     xlab='',ylab='Volatility',
     ylim=c(-2.5,2.5)+log(desv.pad),
     # main='Volatility',
     xaxt='n',
     family=family_font,xaxs='i',yaxs='i')
mtext('(A)',side=3,at=as.Date('1960-01-01'))
axis(1,
     at=seq(as.Date('1965-01-01'),as.Date('1995-01-01'),'10 years'),
     labels=format(seq(as.Date('1965-01-01'),as.Date('1995-01-01'),'10 years'),'%Y'))

base_ribbon2(
  x=date_label[ref.data.prop1$Time] %>% as.Date,
  ymin=ref.data.prop1$ymin+log(desv.pad),
  ymax=ref.data.prop1$ymax+log(desv.pad),
  col='#222222'
)
lines(date_label[ref.data.prop1$Time] %>% as.Date,ref.data.prop1$y+log(desv.pad),col='#222222')
base_ribbon2(
  x=date_label[ref.data.prop2$Time] %>% as.Date,
  ymin=ref.data.prop2$ymin+log(desv.pad),
  ymax=ref.data.prop2$ymax+log(desv.pad),
  col='#aaaaaa'
)
lines(date_label[ref.data.prop2$Time] %>% as.Date,ref.data.prop2$y+log(desv.pad),col='#aaaaaa')

par(mgp=c(1.5,0.5,0),mar = c(3, 2.5, 1, 2.5))
plot(date_label[ref.data.prop2$Time] %>% as.Date,rep(0,T-cutoff),type='n',pch=16,cex=0.5,
     xlab='',ylab='Volatility',
     ylim=c(-2.5,2.5)+log(desv.pad),
     # main='Volatility',
     family=family_font,
     # main='Volatility',
     xaxt='n',
     family=family_font,xaxs='i',yaxs='i')
mtext('(B)',side=3,at=as.Date('1960-01-01'))
axis(1,
     at=seq(as.Date('1965-01-01'),as.Date('1995-01-01'),'10 years'),
     labels=format(seq(as.Date('1965-01-01'),as.Date('1995-01-01'),'10 years'),'%Y'))

base_ribbon2(
  x=date_label[ref.data.prop2$Time] %>% as.Date,
  ymin=colQuantile(-stan_smp$theta,0.025)[-c(1:cutoff)]+log(desv.pad),
  ymax=colQuantile(-stan_smp$theta,0.975)[-c(1:cutoff)]+log(desv.pad),
  col='#222222'
)
lines(date_label[ref.data.prop2$Time] %>% as.Date,colMeans(-stan_smp$theta)[-c(1:cutoff)]+log(desv.pad),col='#222222')
base_ribbon2(
  x=date_label[ref.data.prop2$Time] %>% as.Date,
  ymin=ref.data.prop2$ymin+log(desv.pad),
  ymax=ref.data.prop2$ymax+log(desv.pad),
  col='#aaaaaa'
)
lines(date_label[ref.data.prop2$Time] %>% as.Date,ref.data.prop2$y+log(desv.pad),col='#aaaaaa')


dev.off()
par(mfrow=c(1,1))

# Mean and quantiles of the AR coefficient
c(stan_smp$rho %>% mean, stan_smp$rho %>% quantile(c(0.025,0.975)))%>% round(3) %>% paste(collapse = ' & ')
(model_norm$mts[3,model_norm$t] + c(0,-1,1)*sqrt(model_norm$Cts[3,3,model_norm$t])) %>% round(3) %>% paste(collapse = ' & ')
(model_gamma$mts[3,model_gamma$t] + c(0,-1,1)*sqrt(model_gamma$Cts[3,3,model_gamma$t])) %>% round(3) %>% paste(collapse = ' & ')


