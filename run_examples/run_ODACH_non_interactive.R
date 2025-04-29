
# rm(list=ls())

require(data.table)
require(table1)
require(survival)
require(ggplot2)
# require(lme4)
require(remotes)
install_github('https://github.com/Penncil/pda')
require(pda)


############################## Setup ##############################
## local directory
setwd('/Users/chongliang/Dropbox/PDA-git/UNMC_workshop/run_examples/cloud') # my working dir
mydir = 'OUD_ODACH'       # my DLMM working dir 
dir.create(mydir)   # create DLMM working dir
mysite = 'Kearney'        # my site name

## read in my csv data 
# url_mydata = paste0('https://github.com/chongliang-luo/UNMC_workshop/raw/main/data/OUD_', mysite, '.rds') # github repo
# COVID_LOS_mydata = readRDS(gzcon(url(url_mydata)))    # read in my data from github repo

url_mydata = 'https://github.com/chongliang-luo/UNMC_workshop/raw/main/data/OUD.csv' # github repo
OUD = fread( url_mydata )
########################### END: Setup ##############################


############################ Data summary ###########################
## descriptive
head(OUD)   
table1(~time +factor(status) +factor(age_65) +factor(gender_M) +factor(race_NHW) +
         factor(smoking) +CCI+factor(depression) +factor(pain)|site, data=OUD)

site.name = c('MedicalCenter', 'Lincoln', 'Omaha', 'Kearney')
## let's fit a Cox regression model at each site
bi = sei = c()  
for(sid in site.name){
  fit.i = coxph(Surv(time, status)~age_65+gender_M+race_NHW+smoking+CCI+depression+pain, data=OUD[site==sid,])
  cat(sid, '\n')
  print(round(summary(fit.i)$coef, 3))
  bi = cbind(bi, summary(fit.i)$coef[,1])
  sei = cbind(sei, summary(fit.i)$coef[,3])
}
colnames(bi) = colnames(sei) = site.name
bmeta = round(rowSums(bi/sei^2,na.rm=T)/rowSums(1/sei^2,na.rm=T), 3)
semeta = round(sqrt(1/rowSums(1/sei^2,na.rm=T)), 3)

fit.pooled = coxph(Surv(time, status)~age_65+gender_M+race_NHW+smoking+CCI+depression+pain+strata(site), data=OUD)
round(summary(fit.pooled)$coef, 3)
bpool = round(summary(fit.pooled)$coef[,1], 4)
sepool = round(summary(fit.pooled)$coef[,3], 4)
 
## Discussion...
######################### END: Data summary ###########################

   

  
####################### ODACH workflow ##################################
## setup control
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
OUD_split = split(OUD, by='site')
K = length(site.name)

## ODACH non-interactive section: data communication at local computer
# remove any json files if exist
file.remove(list.files(mydir,full.names = T)[grepl('.json', list.files(mydir))])

# [lead site]: create control.json,
pda(site_id = control$lead_site, control = control, dir = mydir,
    upload_without_confirm =T, silent_message=T)

# [all sites]: initialize step, calculate individual estimates, 
for(sid in K:1) pda(site_id = control$sites[sid], ipdata = OUD_split[[sid]], dir=mydir, 
                    upload_without_confirm =T, silent_message=T)

# [all sites]: derive step, calculate derivatives, 
for(sid in K:1) pda(site_id = control$sites[sid], ipdata = OUD_split[[sid]], dir=mydir, 
                    upload_without_confirm =T, silent_message=T)

# [lead site]: estimate step, calculate ODACH estimate,  
pda(site_id = control$lead_site, dir=mydir, ipdata = OUD_split[[control$lead_site]],
    upload_without_confirm =T, silent_message=T)
config <- getCloudConfig(site_id =control$lead_site, dir=mydir)
fit.odach <- pdaGet(paste0(control$lead_site,'_estimate'), config = config)


# # each site run ODACH estimate and then synthesize them, optional
# pda(site_id = site.name[2], dir=mydir, ipdata = OUD_split[[2]],
#     upload_without_confirm =T, silent_message=T)
# 
# pda(site_id = site.name[3], dir=mydir, ipdata = OUD_split[[3]],
#     upload_without_confirm =T, silent_message=T)
# 
# pda(site_id = site.name[4], dir=mydir, ipdata = OUD_split[[4]],
#     upload_without_confirm =T, silent_message=T)
# 
# pda(site_id = control$lead_site, dir=mydir, ipdata = OUD_split[[control$lead_site]],
#     upload_without_confirm =T, silent_message=T)
# 
# fit.odach <- pdaGet(paste0(site.name[4],'_estimate'), config = config)
# fit.odach <- pdaGet(paste0(control$lead_site,'_synthesize'), config = config)
  
  
  
  
  
## make CI plot
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

ggplot(ci.df, aes(x=method, y=beta, shape=method,color=method, group=risk.factor)) +
  geom_errorbar(aes(ymin=beta-1.96*sd, ymax=beta+1.96*sd), width=0.1) +
  geom_line() +
  geom_hline(aes(yintercept=goldstandard), linetype = "dashed") +
  geom_point(size=1.5) + 
  facet_wrap(. ~ risk.factor
             , scales= "free"  #  "fixed" # 
             , ncol=3)+
  labs(title= case,  
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
 
