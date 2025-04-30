################################################################################
################   UNMC workshop: PDA-OTA practice #3 ODACH  ###################
################################################################################




############################## Setup ##############################
# rm(list=ls()) # empty R， remove all objects

## install packages
require(data.table)
require(table1)
require(survival)
require(ggplot2)
# require(lme4)
require(remotes)
install_github('https://github.com/Penncil/pda')
require(pda)

## local directory
setwd('/Users/chongliang/Dropbox/PDA-git/UNMC_workshop/run_examples/ota_cloud') # my working dir
mydir = 'OUD_ODACH'       # my project working dir 
dir.create(mydir)         # create project working dir
mysite = 'Kearney'        # my site name 'MedicalCenter', 'Lincoln', 'Omaha', 'Kearney'

## read in my csv data  
url_mydata = paste0('https://github.com/chongliang-luo/UNMC_workshop/raw/main/data/OUD_', mysite, '.csv') # github repo
OUD_mydata = fread( url_mydata )
########################### END: Setup ##############################


############################ Data summary ###########################
## descriptive
head(OUD_mydata)   
table1(~time +factor(status) +factor(age_65) +factor(gender_M) +factor(race_NHW) +
         factor(smoking) +CCI+factor(depression) +factor(pain)|site, data=OUD_mydata)
site.name = c('MedicalCenter', 'Lincoln', 'Omaha', 'Kearney')
K = length(site.name)

## let's fit a Cox regression model using my data
fit.i = coxph(Surv(time, status)~age_65+gender_M+race_NHW+smoking+CCI+depression+pain, data=OUD_mydata)
round(summary(fit.i)$coef, 3)

## Discussion...
######################### END: Data summary ###########################
 

####################### ODACH workflow ##################################
## setup control
site.name = c('MedicalCenter', 'Lincoln', 'Omaha', 'Kearney')
control <- list(project_name = 'Opioid use disorder (ODACH)',
                sites = c('MedicalCenter', 'Lincoln', 'Omaha', 'Kearney'),
                lead_site = 'MedicalCenter', 
                step = 'initialize',
                model = 'ODAC',
                family = 'cox',
                heterogeneity = T,
                outcome = 'Surv(time, status)',
                variables = c("age_65","gender_M","race_NHW","smoking","CCI","depression","pain"), 
                init_method = "meta",   
                optim_maxit = 100,
                optim_method = 'BFGS',
                upload_date = as.character(Sys.time()) )

## interactive section:
# [all sites]: remove any json files if exist
file.remove(list.files(mydir,full.names = T)[grepl('.json', list.files(mydir))])
K = length(control$sites)

# [lead site]: create control.json, 
pda(site_id = control$lead_site, control=control, dir=mydir)
# and upload it to pda-ota 

# [all sites]: initialize step, 
# download the control file, calculate individual estimates,   
pda(site_id = mysite, ipdata = OUD_mydata, dir=mydir)
# and upload it to pda-ota

# [all sites]: derive step,  
# download the control file, calculate derivatives,  
pda(site_id = mysite, ipdata = OUD_mydata, dir=mydir)
# and upload it to pda-ota

# [lead site]: estimate step, calculate ODACH estimate, 
pda(site_id = mysite, dir=mydir, ipdata = OUD_mydata)
## END interactive section 


## OK pda-ODACH is completed now, let's check the results
# read in the ODACH results
config <- getCloudConfig(site_id =control$lead_site, dir=mydir)
fit.odach <- pdaGet(paste0(control$lead_site,'_estimate'), config = config)
fit.odach

##################### END: ODACH workflow #################################



########################## Cox reg Pooled data ############################
## read in pooled data 
url_mydata = 'https://github.com/chongliang-luo/UNMC_workshop/raw/main/data/OUD.csv' # github repo
OUD = fread(url_mydata)

## Cox regression, stratified by site
fit.pooled = coxph(Surv(time, status)~age_65+gender_M+race_NHW+smoking+CCI+depression+pain+strata(site), data=OUD)
round(summary(fit.pooled)$coef, 3)
bpool = round(summary(fit.pooled)$coef[,1], 4)
sepool = round(summary(fit.pooled)$coef[,3], 4)


## Cox regression at each site, then average them with inverse-variance as weights. This is called "meta-estimator" 
bi = sei = c()  
for(sid in site.name){
  fit.i = coxph(Surv(time, status)~age_65+gender_M+race_NHW+smoking+CCI+depression+pain, data=OUD[site==sid,])
  cat(sid, '\n')
  print(round(summary(fit.i)$coef, 3))
  bi = cbind(bi, summary(fit.i)$coef[,1])
  sei = cbind(sei, summary(fit.i)$coef[,3])
} 
bmeta = round(rowSums(bi/sei^2,na.rm=T)/rowSums(1/sei^2,na.rm=T), 4)
semeta = round(sqrt(1/rowSums(1/sei^2,na.rm=T)), 4)



## make CI plot: compare pooled vs meta vs ODACH
Xname = c("Age 65+", "Gender male", "Race NHW", "Smoking", "CCI", "Depression", "Pain")
px = length(Xname)
methods <- c('pooled', 'meta', 'ODACH')   
nm <- length(methods)
ci.df <- data.frame(method = rep(methods, each=px),       
                    risk.factor= rep(Xname, nm),
                    beta=as.numeric(c(bpool, bmeta, fit.odach$btilde)),  
                    sd=as.numeric(c(sepool, semeta, fit.odach$setilde)),          
                    goldstandard=as.numeric(rep(bpool, nm)))
ci.df$method <- factor(ci.df$method, levels = rev(methods)) 
ci.df$risk.factor <- factor(ci.df$risk.factor, levels = Xname)  
case = 'Opioid use disorder study'

CI_plot_OUD_ODACH = 
ggplot(ci.df, aes(x=method, y=beta, shape=method,color=method, group=risk.factor)) +
  geom_errorbar(aes(ymin=beta-1.96*sd, ymax=beta+1.96*sd), width=0.1) +
  geom_line() +
  geom_hline(aes(yintercept=goldstandard), linetype = "dashed") +
  geom_point(size=1.5) + 
  facet_wrap(. ~ risk.factor
             , scales= "free"  #  "fixed" # 
             , ncol=3)+
  labs(title= case, x =  '', y = "Effect size (log hazard ratio)" ) +
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
 
ggsave("CI_plot_OUD_ODACH.pdf",  plot = CI_plot_OUD_ODACH, width =10, height =8)

####################### END: Cox reg Pooled data ###########################
 



