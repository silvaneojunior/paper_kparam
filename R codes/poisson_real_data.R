source('R codes/base.R')
# Install NGSSEML from the github repository (https://github.com/hadht/NGSSEML)
# DO NOT INSTALL IT USING install_github!!!
# Download the zip file in the repository and use install.packages:
# install.packages('<path to the zip file>\\NGSSEML_2.2.1.zip')
library(NGSSEML)

# Source:  West and Harrison (1997, Chap. 8, p. 268)
y=c(131.7,322.6,285.6,105.7, 80.4,285.1,347.8, 68.9,203.3,375.9,415.9, 65.8,
    177.0,438.3,463.2,136.0,192.2,442.8,509.6,201.2,196.0,478.6,688.6,259.8,
    352.5,508.1,701.5,325.6,305.9,422.2,771.0,329.3,384.0,472.0,852.0)

date.x=as.Date(paste0(sort(rep(1974:1982,4)),'-',rep(1:4,9)*3,'-1'))[-length(y)-1]

# Poisson case
T <- length(y)
w <- 2*pi/4
data=round(y)

level <- polynomial_block(rate = 1,order=2, D = 0.9,
                          a1=0,R1=matrix(c(2,1,1,1),2,2)/0.9)
season <- harmonic_block(rate = 1, period = 4, order=2, D = 0.95,R1=1/0.95)

A=level+season

# The sample size was chosen to provide a similar precision to NGSSEML
sample_size=1700

fitted_data <- fit_model(level, season,
                         Poisson(lambda = "rate", data = data))

outcome_lb <- Poisson(lambda = "rate", data = data)
outcome_lb$conj_distr=function(ft, Qt, parms) {
  h <- -3 + 3 * sqrt(1 + 2 * Qt / 3)
  alpha <- (1 / h)
  beta <- alpha * exp(-ft + 0.5 * Qt)
  return(list("alpha" = alpha, "beta" = beta))
}

fitted_data_lb <- fit_model(level, season, outcome_lb)

time=1:T
val1_season=cos(w * time)
val2_season=sin(w * time)
val3_season=cos(2*w * time)
dataset=data.frame(data=data,
                   1,time=log(time),val1_season,val2_season,val3_season)

# a0 and b0 were chosen to match the prior using the kDGLM package
fitbayes=ngssm.bayes(data~val1_season+val2_season+val3_season,
                     data = dataset,
                     prbetamu=rep(0,3),
                     prbetasigma=1*diag(3),
                     a0=0.3432753 ,
                     b0=0.039440378199,
                     nsamplex = 5000,
                     model = "Poisson",
                     verbose = TRUE
)

sample_bayes=fitbayes$samplepost

smoothpar = SmoothingF(data~val1_season+val2_season+val3_season,
                     data = dataset,
                     StaPar= sample_bayes,
                     a0=0.3432753 ,
                     b0=0.039440378199,
                     Type = 'Marg',
                     model = "Poisson",
                     samples=5000,
                     splot = TRUE)

# Some convergence tests
# library(coda)
# 
# plot(smoothpar[[2]][18,])
# 
# chain=as.mcmc(t(smoothpar[[2]]))
# summary(chain)
# heidel.diag(chain)
# raftery.diag(chain)
# max(raftery.diag(chain)$resmatrix[,2])
# 
# min(effectiveSize(chain))
# 
# geweke.diag(chain)

pred_vals=coef(fitted_data,lag=-1,eval.pred = TRUE)$data
smooth_pred=pred_vals$Prediction
smooth_ci_lower=pred_vals$C.I.lower
smooth_ci_upper=pred_vals$C.I.upper

pred_vals_lb=coef(fitted_data_lb,lag=-1,eval.pred = TRUE)$data
smooth_pred_lb=pred_vals_lb$Prediction
smooth_ci_lower_lb=pred_vals_lb$C.I.lower
smooth_ci_upper_lb=pred_vals_lb$C.I.upper

((
  ggplot()+
    geom_point(aes(x=date.x,y=data,color='Obs.',fill='Obs.',shape='Obs.'))+
    geom_line(aes(x=date.x,y=rowMeans(smoothpar[[2]]),
                  color='Gamerman et al. (2013)',fill='Gamerman et al. (2013)',shape='Gamerman et al. (2013)'))+
    # geom_line(aes(x=date.x,y=apply(smoothpar[[2]],1,function(x){quantile(x,0.025)}),
    #               color='Gamerman et al. (2013)',fill='Gamerman et al. (2013)',shape='Gamerman et al. (2013)'),
    #           linetype='dashed')+
    # geom_line(aes(x=date.x,y=apply(smoothpar[[2]],1,function(x){quantile(x,0.975)}),
    #               color='Gamerman et al. (2013)',fill='Gamerman et al. (2013)',shape='Gamerman et al. (2013)'),
    #           linetype='dashed')+
    geom_line(aes(x=date.x,y=smooth_pred,
                  color='Our approach',fill='Our approach',shape='Our approach'))+
    # geom_line(aes(x=date.x,y=smooth_ci_lower,
    #               color='Our approach',fill='Our approach',shape='Our approach'),
    #           linetype='dashed')+
    # geom_line(aes(x=date.x,y=smooth_ci_upper,
    #               color='Our approach',fill='Our approach',shape='Our approach'),
    #           linetype='dashed')+
    scale_shape_manual('',values=c('Gamerman et al. (2013)'=NULL,
                                      'Our approach'=NULL,
                                      'Obs.'=16,
                                      'West et al. (1985)'=NULL))+
    scale_color_manual('',values=c('Gamerman et al. (2013)'='#000000',
                                   'Our approach'='#aaaaaa',
                                   'Obs.'=NULL,
                                   'West et al. (1985)'='#555555'))+
    scale_fill_manual('',values=c('Gamerman et al. (2013)'='#000000',
                                  'Our approach'='#aaaaaa',
                                  'Obs.'=NULL,
                                  'West et al. (1985)'='#555555'))+
    scale_x_date('Year')+
    coord_cartesian(ylim=c(0,1000),xlim=as.Date(c('1974-01-01',NA)))+
    scale_y_continuous('Sales',expand=c(0,0),breaks=seq(0,1000,250),labels=function(x){formatC(x,big.mark='.')})+
    theme_bw()+
    theme(text=element_text(size=font_size,family=family_font),
          legend.text=element_text(size=font_size,family=family_font),
          legend.position = "bottom",panel.grid=element_blank())+
    # guides(color = guide_legend(override.aes = list(shape = NA)))+
    guides(color = 'none',shape='none',fill='none')
))

ggsave(
  'Figures/Poisson_plot_1.pdf',
  device = 'pdf',
  units='px',
  width=base.size,
  height=base.size/2
)

plot.data=data.frame(date.x=date.x,
                     y=rowMeans(smoothpar[[2]]/exp(as.matrix(dataset)[,4:6]%*%t(fitbayes$samplepost[,2:4]))),
                     ymin=apply(smoothpar[[2]]/exp(as.matrix(dataset)[,4:6]%*%t(fitbayes$samplepost[,2:4])),1,function(x){quantile(x,0.025)}),
                     ymax=apply(smoothpar[[2]]/exp(as.matrix(dataset)[,4:6]%*%t(fitbayes$samplepost[,2:4])),1,function(x){quantile(x,0.975)}),
                     Model='Gamerman et al. (2013)',
                     plot='(A)')

plot.data=rbind(plot.data,
                data.frame(date.x=date.x,
                           y=exp(fitted_data$mts[1,]+fitted_data$Cts[1,1,]/2),
                           ymin=qlnorm(0.025,fitted_data$mts[1,],sqrt(fitted_data$Cts[1,1,])),
                           ymax=qlnorm(0.975,fitted_data$mts[1,],sqrt(fitted_data$Cts[1,1,])),
                           Model='Our approach',
                           plot='(A)'))


plot.data=rbind(plot.data,
                data.frame(date.x=date.x,
                           y=rowMeans(as.matrix(dataset)[,4:6]%*%t(fitbayes$samplepost[,2:4])),
                           ymin=apply(as.matrix(dataset)[,4:6]%*%t(fitbayes$samplepost[,2:4]),1,function(x){quantile(x,0.025)}),
                           ymax=apply(as.matrix(dataset)[,4:6]%*%t(fitbayes$samplepost[,2:4]),1,function(x){quantile(x,0.975)}),
                           Model='Gamerman et al. (2013)',
                           plot='(B)'))


plot.data=rbind(plot.data,
                data.frame(date.x=date.x,
                           y=fitted_data$mts[3,]+fitted_data$mts[5,]+fitted_data$Cts[3,3,]/2+fitted_data$Cts[5,5,]/2+2*fitted_data$Cts[3,5,]/2,
                           ymin=qnorm(0.025,fitted_data$mts[3,]+fitted_data$mts[5,],sqrt(fitted_data$Cts[3,3,]+fitted_data$Cts[5,5,]+2*fitted_data$Cts[3,5,])),
                           ymax=qnorm(0.975,fitted_data$mts[3,]+fitted_data$mts[5,],sqrt(fitted_data$Cts[3,3,]+fitted_data$Cts[5,5,]+2*fitted_data$Cts[3,5,])),
                           Model='Our approach',
                           plot='(B)'))

((ggplot(plot.data,aes(x=date.x,color=Model,fill=Model))+
    geom_line(aes(y=y))+
    # geom_line(aes(y=ymin),linetype='dashed')+
    # geom_line(aes(y=ymax),linetype='dashed')+
    scale_shape_manual('',values=c('Gamerman et al. (2013)'=NULL,
                                   'Our approach'=NULL,
                                   'Obs.'=16,
                                   'West et al. (1985)'=NULL))+
    scale_color_manual('',values=c('Gamerman et al. (2013)'='#000000',
                                   'Our approach'='#aaaaaa',
                                   'Obs.'=NULL,
                                   'West et al. (1985)'='#555555'))+
    scale_fill_manual('',values=c('Gamerman et al. (2013)'='#000000',
                                  'Our approach'='#aaaaaa',
                                  'Obs.'=NULL,
                                  'West et al. (1985)'='#555555'))+
    scale_x_date('Year',expand=c(0.01,0))+
    coord_cartesian(xlim=as.Date(c('1974-01-01',NA)))+
    scale_y_continuous('Estimated value',expand=c(0.2,0.2))+
    facet_wrap(.~plot,scales='free')+
    theme_bw()+
    theme(text=element_text(size=font_size,family=family_font),
          legend.position = 'bottom',panel.grid=element_blank(),
          legend.text=element_text(size=font_size,family=family_font),
          strip.background = element_blank(),
          strip.clip = 'off',
          strip.text = element_text(hjust = -0.15,vjust=-0.5))))+
  guides(color='none')

ggsave(
  'Figures/Poisson_plot_2.pdf',
  device='pdf',
  units='px',
  width=base.size,
  height=base.size/2.5
)

q_sup=rep(NA,T)
q_inf=rep(NA,T)
log.like=rep(NA,T)
for(i in 1:T){
  sample_lambda=smoothpar[[2]][i,]
  sample_y=rpois(length(sample_lambda),sample_lambda)
  pre.log.like=dpois(data[i],sample_lambda,log=TRUE)
  m=max(pre.log.like)
  log.like[i]=log(mean(exp(pre.log.like-m)))+m
  q_sup[i]=quantile(sample_y,0.975)
  q_inf[i]=quantile(sample_y,0.025)
}

smooth_coef=coef(fitted_data,lag=-1,eval.pred = TRUE)

# Using the approximated preditive distribution of the observations.
log.like_KL=dnbinom(data,
                    smooth_coef$conj.param$Series.1$alpha,
                    smooth_coef$conj.param$Series.1$beta/(smooth_coef$conj.param$Series.1$beta+1),
                    log=TRUE)
                    


xlim=8
xlast=T  
                    
sum(log.like_KL[xlim:xlast])
sum(log.like[xlim:xlast])

ft=fitted_data$ft[1,]
Qt=fitted_data$Qt[1,1,]
h <- -3 + 3 * sqrt(1 + 2 * Qt / 3)
alpha <- (1 / h)
beta <- alpha * exp(-ft - 0.5 * Qt)

smooth_data=smooth_coef$data
smooth_data_lb=coef(fitted_data_lb,lag=-1,eval.pred = TRUE)$data

q_sup_KL=smooth_data$C.I.upper
q_sup_lb=smooth_data_lb$C.I.upper

q_inf_KL=smooth_data$C.I.lower
q_inf_lb=smooth_data_lb$C.I.lower

smooth_pred=smooth_data$Prediction
smooth_pred_lb=smooth_data_lb$Prediction

sum(log.like[xlim:xlast])
sum(log.like_KL[xlim:xlast])


mean(abs((rowMeans(smoothpar[[2]])-data))[xlim:xlast])
mean(abs((smooth_pred-data))[xlim:xlast])
mean(abs((smooth_pred_lb-data))[xlim:xlast])

mean(abs((rowMeans(smoothpar[[2]])-data)/data)[xlim:xlast])
mean(abs((smooth_pred-data)/data)[xlim:xlast])
mean(abs((smooth_pred_lb-data)/data)[xlim:xlast])

mean(abs((rowMeans(smoothpar[[2]])-data)/rowMeans(smoothpar[[2]]))[xlim:xlast])
mean(abs((smooth_pred-data)/mean(abs(diff(data))))[xlim:xlast])
mean(abs((smooth_pred_lb-data)/mean(abs(diff(data))))[xlim:xlast])

mean(abs((rowMeans(smoothpar[[2]])-data))[xlim:xlast]**2)
mean(abs((smooth_pred-data))[xlim:xlast]**2)
mean(abs((smooth_pred_lb-data))[xlim:xlast]**2)


mean(((q_sup-q_inf)*
    2 / (1 - 0.95) * (q_inf - data) * (data < q_inf) +
    2 / (1 - 0.95) * (data - q_sup) * (data > q_sup))[xlim:xlast])
mean(((q_sup_KL-q_inf_KL)*
        2 / (1 - 0.95) * (q_inf_KL - data) * (data < q_inf_KL) +
        2 / (1 - 0.95) * (data - q_sup_KL) * (data > q_sup_KL))[xlim:xlast])
mean(((q_sup_lb-q_inf_lb)*
        2 / (1 - 0.95) * (q_inf_lb - data) * (data < q_inf_lb) +
        2 / (1 - 0.95) * (data - q_sup_lb) * (data > q_sup_lb))[xlim:xlast])

