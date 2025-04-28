
# rm(list=ls())

require(data.table)
require(table1)
# require(survival)
require(ggplot2)
require(lme4)
require(remotes)
install_github('https://github.com/Penncil/pda')
require(pda)

############################## Setup ##############################
## local directory
setwd('UNMC_workshop/run_examples/test/') # my working dir
mydir = 'test_DLMM'       # my DLMM working dir 
dir.create('test_DLMM')   # create DLMM working dir
mysite = 'Kearney'        # my site name

## read in pooled csv data 
# url_mydata = 'https://github.com/chongliang-luo/UNMC_workshop/raw/main/data/COVID_LOS.rds' # github repo
# COVID_LOS = readRDS(gzcon(url(url_mydata)))    # read in my data from github repo

url_mydata = 'https://github.com/chongliang-luo/UNMC_workshop/raw/main/data/COVID_LOS.csv'  # github repo
COVID_LOS = fread( url_mydata )

########################### END: Setup ##############################


############################ Data summary ###########################
## descriptive
head(COVID_LOS)   
table1(~LOS+Age+CCI+Gender+Admission+factor(Cancer)+factor(COPD)+factor(Hyperlipidemia)+factor(Hypertension)+
         factor(KidneyDisease)+factor(Obesity)+factor(HeartDisease)+factor(Diabetes)|site, data=COVID_LOS)

site.name = c('MedicalCenter', 'Lincoln', 'Omaha', 'Kearney')
## let's fit a linear regression model at each site
bi = sei = c()  
for(sid in site.name){
  fit.i = lm(LOS~Age+CCI+Gender+Admission+Cancer+COPD+Hyperlipidemia+Hypertension+
               KidneyDisease + Obesity + HeartDisease + Diabetes, data=COVID_LOS[site==sid,])
  cat(sid, '\n')
  print(round(summary(fit.i)$coef, 3))
  bi = cbind(bi, summary(fit.i)$coef[,1])
  sei = cbind(sei, summary(fit.i)$coef[,2])
}
bmeta = round(rowSums(bi/sei^2,na.rm=T)/rowSums(1/sei^2,na.rm=T), 3)
semeta = round(sqrt(1/rowSums(1/sei^2,na.rm=T)), 3)

## Discussion...
######################### END: Data summary ###########################

 
####################### DLMM workflow ################################## 
## split pooled data by sites
COVID_LOS_split = split(COVID_LOS, by='site') 
## setup control
control <- list(project_name = 'COVID hospitalization length of stay (DLMM)',
                sites = c('MedicalCenter', 'Lincoln', 'Omaha', 'Kearney'),
                lead_site = 'MedicalCenter', 
                step = 'initialize',
                model = 'DLM',
                family = 'gaussian',
                heterogeneity = T,
                heterogeneity_effect = 'random', 
                outcome = 'LOS',
                variables = c("Age","CCI","Gender","Admission","Cancer","COPD","Hyperlipidemia","Hypertension", 
                              "KidneyDisease",  "Obesity","HeartDisease","Diabetes"),
                variables_lev = list(Age=c('18-64','65-80','80+'), 
                                     CCI=c('0-1','2-5','5+'), Admission=c('Q1','Q2','Q3')),
                variables_heterogeneity = c('Intercept'),
                optim_maxit = 100,
                upload_date = as.character(Sys.time()) )
K = length(control$sites)
px = length(bpool) 

# remove any json files if exist
file.remove(list.files(mydir,full.names = T)[grepl('.json', list.files(mydir))])

# [lead site]: create control.json, 
pda(site_id = control$lead_site, control=control, dir=mydir)

# [all sites]: calculate AD, 
for(sid in K:1) pda(site_id = control$sites[sid], ipdata = COVID_LOS_split[[sid]], dir=mydir, 
                    upload_without_confirm =T, silent_message=T)

# [any site]: use AD files to fit DLMM
pda(site_id = control$lead_site, dir=mydir, ipdata = COVID_LOS_split[[control$lead_site]], 
    upload_without_confirm =T, silent_message=T)


## OK pda-DLMM is completed now, let's check the results
# read in the DLMM results
config <- getCloudConfig(site_id =control$lead_site, dir=mydir)
fit.dlmm <- pdaGet(name = paste0(control$lead_site,'_estimate'), config = config)


# fixed effects
data.frame(risk_factor=fit.dlmm$risk_factor, b_dlmm=fit.dlmm$bhat, se_dlmm=fit.dlmm$sebhat)

# var component (of random intercepts)
data.frame(RandomEffect=c(fit.dlmm$risk_factor_heterogeneity,  'random error'),
           VarComp=c(fit.dlmm$Vhat, fit.dlmm$sigmahat^2) )

# random intercepts (BLUP)
data.frame(site = control$sites, BLUP.dlmm = c(fit.dlmm$uhat))


## Next, we test if the Gender effects on LOS are the same across sites
# We don't need to re-do the AD, 

# [any site]: setup the control and re-estimate  
control$variables_heterogeneity = c('Intercept', 'Gender')
control$step = 'estimate'
pdaPut(obj=control, name='control', config=config)
pda(site_id = mysite, dir=mydir, ipdata = COVID_LOS_split[[control$lead_site]])


## Let's check the results 
fit.dlmm1 <- pdaGet(name = paste0(mysite,'_estimate'), config = config)
fit.dlmm1

# Chi-sq test for Gender random effects
1 - pchisq((fit.dlmm1$lik - fit.dlmm$lik)*2, df=1) # p-value

# fixed effects
data.frame(risk_factor=fit.dlmm1$risk_factor, b_dlmm=fit.dlmm1$bhat, se_dlmm=fit.dlmm1$sebhat)

# var components (of random intercepts and Gender effects)
data.frame(RandomEffect=c(fit.dlmm1$risk_factor_heterogeneity,  'random error'),
           VarComp=c(diag(fit.dlmm1$Vhat), fit.dlmm$sigmahat^2) )

# random intercepts and Gender effects (BLUP)
data.frame(site = control$sites, 
           BLUP.dlmm.Intercept = fit.dlmm1$uhat[,1],
           BLUP.dlmm.Gender = fit.dlmm1$uhat[,2])
##################### END: DLMM workflow #################################








########################## LMM Pooled data ############################
## LMM, random intercepts
fit.pooled <- lme4::lmer(LOS~Age+CCI+Gender+Admission+Cancer+COPD+Hyperlipidemia+Hypertension+
                           KidneyDisease + Obesity + HeartDisease + Diabetes +
                           (1|site), REML = F, data=COVID_LOS)
b_pool = round(summary(fit.pooled)$coef[,1], 4)
se_pool = round(summary(fit.pooled)$coef[,2], 4)

# fixed effects
data.frame(risk_factor=fit.dlmm$risk_factor, b_pool=b_pool, b_meta=bmeta, b_dlmm=fit.dlmm$bhat)

# var component
data.frame(RandomEffect=c(fit.dlmm$risk_factor_heterogeneity,  'random error'),
           VarComp.dlmm=data.frame(summary(fit.pooled)$varcor)$vcov, 
           VarComp.pool=c(fit.dlmm$Vhat, fit.dlmm$sigmahat^2) )

# random intercepts (BLUP) 
data.frame(site = control$sites, 
           BLUP.pool.Intercept = round(ranef(fit.pooled)$site$`(Intercept)`, 4),
           BLUP.dlmm.Intercept = fit.dlmm$uhat[,1])



## LMM, random intercepts and Gender effects
COVID_LOS[,Gender:=ifelse(Gender=='F',0,1)] # create dummy var 
fit.pooled1 <- lme4::lmer(LOS~Age+CCI+Gender+Admission+Cancer+COPD+Hyperlipidemia+Hypertension+
                            KidneyDisease + Obesity + HeartDisease + Diabetes +
                            (1|site)+(0+Gender|site), REML = F, data=COVID_LOS)
# compare p-value:
# Chi-sq test for Gender random effect (likelihood ratio test, LRT)
anova(fit.pooled, fit.pooled1) 
# Chi-sq test for Gender random effect, via DLMM 
1 - pchisq((fit.dlmm1$lik - fit.dlmm$lik)*2, df=1) # p-value

# fixed effects
b_pool1 = round(summary(fit.pooled1)$coef[,1], 4)
se_pool1 = round(summary(fit.pooled1)$coef[,2], 4)
data.frame(risk_factor=fit.dlmm1$risk_factor, b_pool=b_pool1, b_meta=bmeta, b_dlmm=fit.dlmm1$bhat)

# var component (of random intercepts and Gender effects)
data.frame(RandomEffect=c(fit.dlmm1$risk_factor_heterogeneity,  'random error'),
           VarComp.dlmm=data.frame(summary(fit.pooled1)$varcor)$vcov, 
           VarComp.pool=c(diag(fit.dlmm1$Vhat), fit.dlmm1$sigmahat^2) )

# random intercepts and Gender effects (BLUP)
data.frame(site = control$sites, 
           BLUP.pool.Intercept = round(ranef(fit.pooled1)$site$`(Intercept)`, 4),
           BLUP.dlmm.Intercept = fit.dlmm1$uhat[,1], 
           BLUP.pool.GenderM = round(ranef(fit.pooled1)$site$Gender, 4),
           BLUP.dlmm.GenderM = fit.dlmm1$uhat[,2])

####################### END: LMM Pooled data ############################


