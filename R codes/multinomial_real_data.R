source('R codes/base.R')

# Age group labels
age.groups = list(
  'Menor 1 ano' = 'Less than 1 year',
  '1 a 4 anos' = '1 to 4 years',
  '5 a 9 anos' = '5 to 9 years',
  '10 a 14 anos' = '10 to 14 years',
  '15 a 19 anos' = '15 to 19 years',
  '20 a 29 anos' = '20 to 29 years',
  '30 a 39 anos' = '30 to 39 years',
  '40 a 49 anos' = '40 to 49 years',
  '50 a 59 anos' = '50 to 59 years',
  '60 a 69 anos' = '60 to 69 years',
  '70 a 79 anos' = '70 to 79 years',
  '80 anos e mais' = '80 years and over'
)

month.translation <- c(
  "Jan" = "Jan", "Fev" = "Feb", "Mar" = "Mar", "Abr" = "Apr",
  "Mai" = "May", "Jun" = "Jun", "Jul" = "Jul", "Ago" = "Aug",
  "Set" = "Sep", "Out" = "Oct", "Nov" = "Nov", "Dez" = "Dec"
)

convert.age=function(x){
  age.groups[[x]]
}
convert.month=function(x){
  month.translation[[x]]
}

# From 1998 to 2007
data.2007=read.csv('data/gastro_1998_2007.csv') %>%
  select(-Total) %>%
  filter(Faixa.Etaria.1!='Total') %>%
  pivot_longer(-Faixa.Etaria.1)

# From 2008 to present (at the time of extraction)
data.2008=read.csv('data/gastro_2008_2023.csv') %>%
  select(-Total) %>%
  filter(Faixa.Etaria.1!='Total') %>%
  pivot_longer(-Faixa.Etaria.1)

data.raw= data.2007 %>%
  full_join(data.2008,by=c('Faixa.Etaria.1','name')) %>%
  mutate(value=if.na(value.x,0)+if.na(value.y,0),
         value.x=NULL,
         value.y=NULL)

# The proportion of children less than 4 years old that are less then 1 year old (see bellow).
pre.prorp=rowSums(read.csv('data/populacao_1980_2012.csv')[c(2,3),-1])
prorp.inft=pre.prorp[1]/(pre.prorp[1]+pre.prorp[2])

data=data.raw %>%
  mutate(Faixa.Etaria.1=sapply(Faixa.Etaria.1,convert.age),
         name=paste0(substr(name,2,5),'-',sapply(substr(name,7,9),convert.month),'-01') %>% 
           as.Date(format='%Y-%b-%d')) %>%
  rename(Age.group=Faixa.Etaria.1,Date=name,Hosp.admt=value)

# The data for the population size groups does not contain the age groups "Less than 1 year old" and "1 to 4 years old"
# Instead, those age groups are combined into a single age group: "Less than 4 year old".
# Using an older dataset that provides information only up to 2012, we observed that the proportion of children
# under 4 years old who are also younger than 1 year old has remained stable over the years.
# Based on this observation, we decided to use this proportion to split the "Less than 4 years old" age group 
# into "Less than 1 year old" and "1 to 4 years old".
data.pop=read.csv('data/populacao_2000_2023.csv')[-12,]

less.than.4=data.pop[1,-1]
less.than.1=cbind(Faixa.Etaria='Menor 1 ano',less.than.4*prorp.inft)
one.2.four=cbind(Faixa.Etaria='1 a 4 anos',less.than.4*(1-prorp.inft))

data.pop=rbind(less.than.1,one.2.four,data.pop[-1,])

data.pop=data.pop %>%
  pivot_longer(-1) %>%
  mutate(name=substr(name,2,5)) %>%
  rename(Year=name, Exp=value,Age.group=Faixa.Etaria) %>%
  mutate(Age.group=sapply(Age.group,convert.age)) %>%
  group_by(Age.group,Year) %>%
  summarize_at('Exp',sum) %>%
  arrange(Age.group,Year) %>%
  ungroup()

data=data %>% 
  mutate(Year=Date %>% as.Date(format='%Y-%b-%d') %>% format('%Y')) %>% 
  inner_join(data.pop,by=c('Age.group','Year'))

y=data %>%
  filter(Year<2023) %>%
  select(Age.group,Date,Hosp.admt) %>%
  pivot_wider(names_from=Age.group,values_from = Hosp.admt) %>% 
  select(-Date)
exp=data %>%
  filter(Year<2023) %>%
  select(Age.group,Date,Exp) %>%
  pivot_wider(names_from=Age.group,values_from = Exp) %>% 
  select(-Date)



clusters=(y/exp)%>%  as.matrix %>% log  %>% diff(lag=12) %>% {1-cor(.)} %>% as.dist %>% 
  hclust() %>% {plot(.);cutree(.,k=4)}

# From the dendogram we notice 4 major age groups: less than 4 years old, 5 to 9 years old, 15 to 49 years old and 50 years or older.
# We chose to keep the age group of Less than 1 year separate from the age group of 1 to 4 years old, because of the introduction of the vaccine (see the paper).

clusters[1]='Less than 1 year'
clusters[2]='1 to 4 years'
clusters[clusters==2]='5 to 14 years'
clusters[clusters==3]='Reference group'
clusters[clusters==4]='60+ years'

convert.age=function(x){
  clusters[[x]]
}

data=data %>% mutate(Age.group=sapply(Age.group,convert.age) %>% factor(.,unique(.))) %>% 
  group_by(Age.group,Date,Year) %>% 
  summarize_all(sum) %>% 
  ungroup()

y=data %>%
  filter(Year<2023) %>%
  select(Age.group,Date,Hosp.admt) %>%
  pivot_wider(names_from=Age.group,values_from = Hosp.admt) %>% 
  select(-Date)
exp=data %>%
  filter(Year<2023) %>%
  select(Age.group,Date,Exp) %>%
  pivot_wider(names_from=Age.group,values_from = Exp) %>% 
  select(-Date)

k <- dim(y)[2]
T <- dim(y)[1]

ggplot(data,aes(x=Date,y=Hosp.admt/Exp,color=Age.group))+
  scale_y_log10()+
  geom_point()+
  geom_vline(xintercept=as.Date('2008-01-01'),linetype='dashed')+
  theme_bw()


level=(polynomial_block(p=1,order=2,D=0.95,name='Level')+noise_block(p=1,H=0.01,name='Noise')) %>% 
  block_mult(k-1) %>% 
  block_rename(unique(data$Age.group) %>% {.[.!='Reference group']} %>% as.character)
season=harmonic_block(p=1,period=12,D=0.975,name='Season') %>% 
          block_mult(k-1) %>% 
          block_rename(level$pred.names)
outcome=Multinom(p=level$pred.names,data = y,offset=exp,base.class='Reference group')

start=Sys.time()
fitted_data=fit_model(level,season,outcome)
difftime(Sys.time(),start)
# fitted_data
# microbenchmark::microbenchmark(fit_model(level,season,outcome))

# plot(fitted_data)
# plot(fitted_data,latent.states = '.Level')

predictions=coef(fitted_data,lag=1,eval.pred = TRUE)

obs.pred=predictions$data %>% 
  mutate(Serie=Serie %>% substr(10,80) %>% factor(.,levels=unique(.)),
         Time=unique(data$Date)[Time]) %>% 
  rename(Date=Time,Age.group=Serie)

names(clusters)=paste0('alpha_',1:5)
convert.age=function(x){
  clusters[[x]]
}

prorp.pred=obs.pred %>% 
  group_by(Date) %>% 
  mutate(Prediction =Prediction /sum(Observation),
         C.I.lower =C.I.lower /sum(Observation),
         C.I.upper =C.I.upper /sum(Observation),
         Variance  =Variance  /(sum(Observation)**2),
         Observation =Observation /sum(Observation))

tx.pred=obs.pred %>% 
  left_join(data %>% select(Date,Age.group,Exp),by=c('Date','Age.group')) %>% 
  mutate(Prediction =Prediction /Exp,
         C.I.lower =C.I.lower /Exp,
         C.I.upper =C.I.upper /Exp,
         Variance  =Variance  /(Exp**2),
         Observation =Observation /Exp) %>% 
  select(-Exp)

plot.data=rbind(cbind(obs.pred,plot='(A)'), #Total Number
                cbind(prorp.pred,plot='(B)'),# P(Age group | Admission)
                cbind(tx.pred,plot='(C)')) %>%  # P(Admission | Age group)
  mutate(plot=plot %>% factor(.))

colors=c('Less than 1 year'='#555555',
         '1 to 4 years'='black',
         '5 to 14 years'='#555555',
         '60+ years'='#aaaaaa')
shape=c('Less than 1 year'='dashed',
         '1 to 4 years'='solid',
         '5 to 14 years'='solid',
         '60+ years'='solid')

ggplot(plot.data %>% filter(format(Date,'%Y')>2002,plot!='(C)',Age.group%in% c('Less than 1 year','1 to 4 years')),aes(x=Date,color=Age.group,fill=Age.group,linetype=Age.group))+
    geom_line(aes(y=Prediction))+
    geom_point(aes(y=Observation),alpha=0.5,size=0.5)+
    geom_ribbon(aes(ymin=C.I.lower,
                    ymax=ifelse(C.I.upper>10000,10000,C.I.upper)),alpha=0.25,color=NA)+
    geom_vline(xintercept=as.Date('2006-03-01'),linetype='dashed')+
  scale_color_manual('',values=colors)+
  scale_linetype_manual('',values=shape)+
    scale_fill_manual('',values=colors)+
    scale_shape('Observation')+
    scale_y_continuous('',n.breaks=7,labels=function(x){
      ifelse(x<=0.05 & x!=0,
             round(100*x,2)%>%paste0('%'),
             ifelse(x<=1 & x!=0,
                    round(100*x)%>%paste0('%'),
                    formatC(x,big.mark=',',format='d')
                    )
             )
    })+
    scale_x_date('Year',expand=c(0,0))+
    facet_wrap(.~plot,scales='free')+
    theme_bw()+
    theme(text=element_text(size=font_size,family=family_font),
          legend.position = 'bottom',panel.grid=element_blank(),
          legend.text=element_text(size=font_size,family=family_font),
          strip.background = element_blank(),
          strip.clip = 'off',
          strip.text = element_text(hjust = -0.15,vjust=-0.5))
  
ggsave(
  'Figures/multinomial_real_data.pdf',
  device='pdf',
  units='px',
  width=base.size,
  height=base.size/2
)
p=plotly::ggplotly(ggplot(plot.data %>% filter(format(Date,'%Y')>2002,plot=='(A)'),aes(x=Date,color=Age.group,fill=Age.group))+
  geom_line(aes(y=Prediction))+
  geom_point(aes(y=Observation),alpha=0.5,size=0.5)+
  geom_ribbon(aes(ymin=C.I.lower,
                  ymax=ifelse(C.I.upper>10000,10000,C.I.upper)),alpha=0.25,linetype=0)+
  geom_vline(xintercept=75,linetype='dashed')+
  scale_color_hue('Prediction')+
  scale_fill_hue('Prediction')+
  scale_y_continuous('',expand=c(0,0,0,0),n.breaks=7,labels=function(x){
    ifelse(x<=0.05,
           round(100*x,2)%>%paste0('%'),
           ifelse(x<=1,
                  round(100*x)%>%paste0('%'),
                  formatC(x,big.mark=',',format='d')
           )
    )
  })+
  scale_x_date('Year',expand=c(0,0))+
  theme_bw()+
  theme(text=element_text(size=font_size,family=family_font),legend.position = 'bottom',panel.grid=element_blank(),
        strip.background = element_blank(),
        strip.clip = 'off',
        strip.text = element_text(hjust = -0.25,vjust=-0.5)))

path=getwd()
setwd('Multinomial_ Interative_graphs/')
htmlwidgets::saveWidget(plotly::as_widget(p), "./total.html",selfcontained = TRUE)
p=plotly::ggplotly(ggplot(plot.data %>% filter(format(Date,'%Y')>2002,plot=='(B)'),aes(x=Date,color=Age.group,fill=Age.group))+
                     geom_line(aes(y=Prediction))+
                     geom_point(aes(y=Observation),alpha=0.5,size=0.5)+
                     geom_ribbon(aes(ymin=C.I.lower,
                                     ymax=ifelse(C.I.upper>10000,10000,C.I.upper)),alpha=0.25,linetype=0)+
                     geom_vline(xintercept=75,linetype='dashed')+
                     scale_color_hue('Prediction')+
                     scale_fill_hue('Prediction')+
                     scale_y_continuous('',expand=c(0,0,0,0),n.breaks=7,labels=function(x){
                       ifelse(x<=0.05,
                              round(100*x,2)%>%paste0('%'),
                              ifelse(x<=1,
                                     round(100*x)%>%paste0('%'),
                                     formatC(x,big.mark=',',format='d')
                              )
                       )
                     })+
                     scale_x_date('Year',expand=c(0,0))+
                     theme_bw()+
                     theme(text=element_text(size=font_size,family=family_font),legend.position = 'bottom',panel.grid=element_blank(),
                           strip.background = element_blank(),
                           strip.clip = 'off',
                           strip.text = element_text(hjust = -0.25,vjust=-0.5)))

htmlwidgets::saveWidget(plotly::as_widget(p), "./composition.html",selfcontained = TRUE)
p=plotly::ggplotly(ggplot(plot.data %>% filter(format(Date,'%Y')>2002,plot=='(C)'),aes(x=Date,color=Age.group,fill=Age.group))+
                     geom_line(aes(y=Prediction))+
                     geom_point(aes(y=Observation),alpha=0.5,size=0.5)+
                     geom_ribbon(aes(ymin=C.I.lower,
                                     ymax=ifelse(C.I.upper>10000,10000,C.I.upper)),alpha=0.25,linetype=0)+
                     geom_vline(xintercept=75,linetype='dashed')+
                     scale_color_hue('Prediction')+
                     scale_fill_hue('Prediction')+
                     scale_y_continuous('',expand=c(0,0,0,0),n.breaks=7,labels=function(x){
                       ifelse(x<=0.05,
                              round(100*x,2)%>%paste0('%'),
                              ifelse(x<=1,
                                     round(100*x)%>%paste0('%'),
                                     formatC(x,big.mark=',',format='d')
                              )
                       )
                     })+
                     scale_x_date('Year',expand=c(0,0))+
                     theme_bw()+
                     theme(text=element_text(size=font_size,family=family_font),legend.position = 'bottom',panel.grid=element_blank(),
                           strip.background = element_blank(),
                           strip.clip = 'off',
                           strip.text = element_text(hjust = -0.25,vjust=-0.5)))
htmlwidgets::saveWidget(plotly::as_widget(p), "./probability.html",selfcontained = TRUE)
setwd(path)
