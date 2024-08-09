##############################################################################
source('base.r')

colors=c('Our proposal'='#aaaaaadd',
         'STAN'='#000000dd'
)

#### Time spent ####
time.data=read.csv('section 3-1 simulation/time_data_multinomial.csv') %>%
  group_by(model,n.obs,N) %>%
  summarize(time_spent=mean(time_spent)) %>%
  mutate(model=ifelse(model=='kdglm','Our proposal','STAN') %>%
           factor(levels=c('STAN','Our proposal')))

# Avarage time per sample unit
time.data %>%
  group_by(model) %>%
  summarize(time_spent=sum(time_spent), n.obs=sum(n.obs)) %>%
  mutate(time_avg=time_spent/n.obs)

plot.data=rbind(cbind(time.data,scale='Time~~to~~fit~~(seconds)'),
                cbind(time.data %>%
                        mutate(time_spent=log10(time_spent),
                               n.obs=n.obs),
                      scale='Time~~to~~fit~~(log[10]~~seconds)'))

plot.data$scale=factor(plot.data$scale,levels=unique(plot.data$scale))


ggplot(plot.data %>% filter(scale=='Time~~to~~fit~~(log[10]~~seconds)'))+
  geom_point(aes(x=n.obs,y=time_spent,color=model))+
  geom_line(aes(x=n.obs,y=time_spent,color=model))+
  scale_y_continuous('')+
  scale_x_continuous('Sample size')+
  scale_color_manual('Method',values=colors)+
  facet_grid(scale~N,
             labeller = label_parsed)+
  theme_bw()+
  theme(text=element_text(size=font_size,family=family_font),
        legend.position = 'bottom',panel.grid=element_blank(),
        strip.text = element_text(size = 3*font_size/4, margin = margin()),
        strip.background = element_blank(),
        strip.placement = "outside")

#####################################

simul.data=read.csv('section 3-1 simulation/simul_data_multinomial.csv') %>%
  mutate(model=ifelse(model=='kdglm','Our proposal','STAN'))

for(col in c('model','var')){
  simul.data[[col]]=as.factor(simul.data[[col]])
}
for(col in c('y','time','sample','n.obs','N','mean','icl','icl')){
  simul.data[[col]]=as.numeric(simul.data[[col]])
}

simul.data$model %>% unique()
names(simul.data)


#### Short versions ####
width=2*base.size
height=2*base.size

ggplot(simul.data%>%
         arrange(desc(model))%>%
         filter(((var=='y.1' | var=='y.2') & time-n.obs==1) | ((var=='theta.1' | var=='theta.2') & time==n.obs), n.obs==10) %>%
         mutate(var=factor(ifelse(var=='y.1','y[1]',ifelse(var=='y.2','y[2]',ifelse(var=='theta.1','theta[1]','theta[2]'))),
                           levels=c('y[1]','y[2]','theta[1]','theta[2]'))) %>%
         # mutate(mean=(mean-y)/y,icl=(icl-y)/y,icu=(icu-y)/y) %>%
         mutate(error=(mean-y)/N) %>%
         select(model,n.obs,var,N,time,error) %>%
         group_by(model,n.obs,time,var,N) %>%
         mutate(error=error,n.obs=as.factor(n.obs),N=as.factor(N)) %>%
         summarize(mean=mean(error),icl=quantile(error,0.025),icu=quantile(error,0.975)))+
  geom_hline(yintercept = 0,linetype='dashed')+
  geom_point(aes(x=N,y=mean,color=model,shape=model))+
  geom_errorbar(aes(x=N,ymin=icl,ymax=icu,color=model))+
  scale_x_discrete('N')+
  # scale_y_continuous('Prediction error',limits=c(-0.75,0.75),labels=~paste0(100*.,'%'))+
  scale_y_continuous('Bias')+
  scale_color_manual('Method',values=colors)+
  scale_shape_manual('Method',values=c(1,2,0,4))+
  # guides(shape='none')+
  facet_wrap(~var,labeller = label_parsed)+
  theme_bw()+
  theme(text=element_text(size=font_size,family=family_font),
        legend.position = 'bottom',panel.grid=element_blank(),
        strip.text = element_text(size = 3*font_size/4, margin = margin()))

ggsave(
  'figures/fig_paper_alves1_short.pdf',
  device='pdf',
  units='px',
  width=width,
  height=height
)

ggplot(simul.data %>%
         arrange(desc(model))%>%
         filter(((var=='y.1' | var=='y.2') & time-n.obs==1) | ((var=='theta.1' | var=='theta.2') & time==n.obs),
                N>2 | var=='theta.1' | var=='theta.2',n.obs==10) %>%
         mutate(var=factor(ifelse(var=='y.1','y[1]',ifelse(var=='y.2','y[2]',ifelse(var=='theta.1','theta[1]','theta[2]'))),
                           levels=c('y[1]','y[2]','theta[1]','theta[2]'))) %>%
         mutate(error=(icu-icl)/N) %>%
         select(model,n.obs,var,N,time,error) %>%
         group_by(model,n.obs,time,var,N) %>%
         mutate(error=error,n.obs=as.factor(n.obs),N=as.factor(N)) %>%
         summarize(mean=mean(error),icl=quantile(error,0.025),icu=quantile(error,0.975)))+
  geom_point(aes(x=N,y=mean,color=model,shape=model))+
  scale_x_discrete('N')+
  scale_y_continuous('Relative range')+
  scale_color_manual('Method',values=colors)+
  scale_shape_manual('Method',values=c(1,2,0,4))+
  # guides(shape='none')+
  facet_wrap(~var,
             labeller = label_parsed,scales='free')+
  theme_bw()+
  theme(text=element_text(size=font_size,family=family_font),
        legend.position = 'bottom',panel.grid=element_blank(),
        strip.text = element_text(size = 3*font_size/4, margin = margin()))

ggsave(
  'figures/fig_paper_alves2_short.pdf',
  device='pdf',
  units='px',
  width=width,
  height=height
)

#### Long versions ####
width=2*base.size
height=4*base.size

ggplot(simul.data%>%
         arrange(desc(model))%>%
         filter(var=='y.1' | var=='y.2' | var=='y.3',time-n.obs==1) %>%
         mutate(var=ifelse(var=='y.1','y[1]',ifelse(var=='y.2','y[2]','y[3]'))) %>%
         # mutate(mean=(mean-y)/y,icl=(icl-y)/y,icu=(icu-y)/y) %>%
         mutate(error=(mean-y)/N) %>%
         select(model,n.obs,var,N,time,error) %>%
         group_by(model,n.obs,time,var,N) %>%
         mutate(error=error,n.obs=as.factor(n.obs),N=as.factor(N)) %>%
         summarize(mean=mean(error),icl=quantile(error,0.025),icu=quantile(error,0.975)))+
  geom_hline(yintercept = 0,linetype='dashed')+
  geom_point(aes(x=N,y=mean,color=model,shape=model))+
  geom_errorbar(aes(x=N,ymin=icl,ymax=icu,color=model))+
  scale_x_discrete('N')+
  # scale_y_continuous('Prediction error',limits=c(-0.75,0.75),labels=~paste0(100*.,'%'))+
  scale_y_continuous('Prediction error')+
  scale_color_manual('Method',values=colors)+
  scale_shape_manual('Method',values=c(1,2,0,4))+
  # guides(shape='none')+
  facet_grid((paste0(n.obs,'~observations') %>% factor(.,levels=unique(.)))~var,labeller = label_parsed)+
  theme_bw()+
  theme(text=element_text(size=font_size,family=family_font),
        legend.position = 'bottom',panel.grid=element_blank(),
        strip.text = element_text(size = 3*font_size/4, margin = margin()))

ggsave(
  'figures/fig_paper_alves1.pdf',
  device='pdf',
  units='px',
  width=width,
  height=height
)

ggplot(simul.data %>%
         arrange(desc(model))%>%
         filter(var=='y.1' | var=='y.2' | var=='y.3',time-n.obs==1,N>2) %>%
         mutate(var=ifelse(var=='y.1','y[1]',ifelse(var=='y.2','y[2]','y[3]'))) %>%
         mutate(error=(icu-icl)/N) %>%
         select(model,n.obs,var,N,time,error) %>%
         group_by(model,n.obs,time,var,N) %>%
         mutate(error=error,n.obs=as.factor(n.obs),N=as.factor(N)) %>%
         summarize(mean=mean(error),icl=quantile(error,0.025),icu=quantile(error,0.975)))+
  geom_point(aes(x=N,y=mean,color=model,shape=model))+
  scale_x_discrete('N')+
  scale_y_continuous('Relative range')+
  scale_color_manual('Method',values=colors)+
  scale_shape_manual('Method',values=c(1,2,0,4))+
  # guides(shape='none')+
  facet_grid((paste0(n.obs,'~observations') %>% factor(.,levels=unique(.)))~var,
             labeller = label_parsed)+
  theme_bw()+
  theme(text=element_text(size=font_size,family=family_font),
        legend.position = 'bottom',panel.grid=element_blank(),
        strip.text = element_text(size = 3*font_size/4, margin = margin()))


ggsave(
  'figures/fig_paper_alves2.pdf',
  device='pdf',
  units='px',
  width=width,
  height=height
)

ggplot(simul.data%>%
         arrange(desc(model))%>%
         filter(var=='theta.1' | var=='theta.2' | var=='theta.3',time==n.obs) %>%
         mutate(var=ifelse(var=='theta.1','theta[1]',ifelse(var=='theta.2','theta[2]','theta[3]'))) %>%
         # mutate(mean=(mean-y)/y,icl=(icl-y)/y,icu=(icu-y)/y) %>%
         mutate(error=(mean-y)) %>%
         select(model,n.obs,var,N,time,error) %>%
         group_by(model,n.obs,time,var,N) %>%
         mutate(error=error,n.obs=as.factor(n.obs),N=as.factor(N)) %>%
         summarize(mean=mean(error),icl=quantile(error,0.025),icu=quantile(error,0.975)))+
  geom_hline(yintercept = 0,linetype='dashed')+
  geom_point(aes(x=N,y=mean,color=model,shape=model))+
  geom_errorbar(aes(x=N,ymin=icl,ymax=icu,color=model))+
  scale_x_discrete('N')+
  # scale_y_continuous('Prediction error',limits=c(-0.75,0.75),labels=~paste0(100*.,'%'))+
  scale_y_continuous('Prediction error')+
  scale_color_manual('Method',values=colors)+
  scale_shape_manual('Method',values=c(1,2,0,4))+
  # guides(shape='none')+
  facet_grid((paste0(n.obs,'~observations') %>% factor(.,levels=unique(.)))~var,labeller = label_parsed)+
  theme_bw()+
  theme(text=element_text(size=font_size,family=family_font),
        legend.position = 'bottom',panel.grid=element_blank(),
        strip.text = element_text(size = 3*font_size/4, margin = margin()))

ggsave(
  'figures/fig_paper_alves3.pdf',
  device='pdf',
  units='px',
  width=width,
  height=height
)

ggplot(simul.data %>%
         arrange(desc(model))%>%
         filter(var=='theta.1' | var=='theta.2' | var=='theta.3',time==n.obs) %>%
         mutate(var=ifelse(var=='theta.1','theta[1]',ifelse(var=='theta.2','theta[2]','theta[3]'))) %>%
         mutate(error=icu-icl) %>%
         select(model,n.obs,var,N,time,error) %>%
         group_by(model,n.obs,time,var,N) %>%
         mutate(error=error,n.obs=as.factor(n.obs),N=as.factor(N)) %>%
         summarize(mean=mean(error),icl=quantile(error,0.025),icu=quantile(error,0.975)))+
  geom_point(aes(x=N,y=mean,color=model,shape=model))+
  scale_x_discrete('N')+
  scale_y_continuous('Relative range')+
  scale_color_manual('Method',values=colors)+
  scale_shape_manual('Method',values=c(1,2,0,4))+
  # guides(shape='none')+
  facet_grid((paste0(n.obs,'~observations') %>% factor(.,levels=unique(.)))~var,
             labeller = label_parsed)+
  theme_bw()+
  theme(text=element_text(size=font_size,family=family_font),
        legend.position = 'bottom',panel.grid=element_blank(),
        strip.text = element_text(size = 3*font_size/4, margin = margin()))


ggsave(
  'figures/fig_paper_alves4.pdf',
  device='pdf',
  units='px',
  width=width,
  height=height
)

########################################################################################
