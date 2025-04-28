


OUD = readRDS('/Users/chongliang/Dropbox/PDA-git/UNMC_workshop/data/OUD.rds')
site.name = c('MedicalCenter', 'Lincoln', 'Omaha', 'Kearney')
head(OUD)  
table1(~time +factor(status) +age +factor(gender) +factor(RACE_NHW) +
         factor(smoking) +CCI+factor(depression) +factor(pain)|site, data=OUD)

fit.pooled = coxph(Surv(time, status)~age+gender+RACE_NHW+smoking+CCI+depression+pain+strata(site), data=OUD)
round(summary(fit.pooled)$coef, 3)
bpool = round(summary(fit.pooled)$coef[,1], 4)
sepool = round(summary(fit.pooled)$coef[,3], 4)

px = length(bpool)
K = length(site.name)
bi = sei = c() # matrix(NA, px, K)
for(sid in site.name){
  fit.i = coxph(Surv(time, status)~age+gender+RACE_NHW+smoking+CCI+depression+pain, data=OUD[site==sid,])
  cat(sid, '\n')
  print(round(summary(fit.i)$coef, 3))
  bi = cbind(bi, summary(fit.i)$coef[,1])
  sei = cbind(sei, summary(fit.i)$coef[,3])
} 
bmeta = round(rowSums(bi/sei^2,na.rm=T)/rowSums(1/sei^2,na.rm=T), 4)
semeta = round(sqrt(1/rowSums(1/sei^2,na.rm=T)), 4)
cbind(bpool, bmeta) 

control <- list(project_name = 'Opioid use disorder (ODACH)',
                sites = c('MedicalCenter', 'Lincoln', 'Omaha', 'Kearney'),
                lead_site = 'MedicalCenter', 
                step = 'initialize',
                model = 'ODAC',
                family = 'cox',
                heterogeneity = T,
                outcome = 'Surv(time, status)',
                variables = c("age","gender","RACE_NHW","smoking","CCI","depression","pain"),
                # variables_lev = list(Race=c('White', 'Black', 'Asian', 'Other')),
                init_method = "meta",   
                optim_maxit = 100,
                optim_method = 'BFGS',
                upload_date = as.character(Sys.time()) )
OUD_split = split(OUD, by='site')
res = run_ODACH_with_pda(control, mydir=getwd(), OUD_split, upload_without_confirm=T, silent_message=T)
res$btilde

Xname = c("Age 65+", "Gender male", "Race NHW", "Smoking", "CCI", "Depression", "Pain")
px = length(Xname)
methods <- c('pooled', 'meta', 'ODACH')   
nm <- length(methods)
ci.df <- data.frame(method = rep(methods, each=px),       
                    risk.factor= rep(Xname, nm),
                    beta=as.numeric(c(bpool, bmeta, res$btilde)),  
                    sd=as.numeric(c(sepool, semeta, res$setilde)),          
                    goldstandard=as.numeric(rep(bpool, nm)))
ci.df$method <- factor(ci.df$method, levels = rev(methods)) 
ci.df$risk.factor <- factor(ci.df$risk.factor, levels = Xname)  
case = 'Opioid use disorder study'

ggplot(ci.df, aes(x=method, y=beta, shape=method,color=method, group=risk.factor)) +
  geom_errorbar(aes(ymin=beta-1.96*sd, ymax=beta+1.96*sd), width=0.1) +
  geom_line() +
  geom_hline(aes(yintercept=goldstandard), linetype = "dashed") +
  geom_point(size=1.5) + 
  facet_wrap(. ~ risk.factor
             , scales= "free"  #  "fixed" # 
             , ncol=3)+
  labs(title= case, # case, # paste0(case,', 95% confidence intervals of effect size estimates'),
       x =  '', y = "Effect size (log hazard ratio)" ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)
        , axis.ticks.x = element_blank()
        , strip.text.x = element_text(size = 14)
        , panel.grid.major = element_blank() 
        , axis.title.y = element_text(size=9)
        , axis.text.x = element_text(size=9)
        , axis.text.y = element_text(size=9) 
        , legend.position = c(0.7, 0.13)
        , legend.title = element_blank() 
        , legend.text=element_text(size=15) 
        , legend.key.height = unit(1.5, "line")
  ) + coord_flip()
 
