################################################################################
################   UNMC workshop: PDA-OTA practice #1 DLMM  ####################
################################################################################



############################## Setup ##############################
# rm(list=ls()) # empty R， remove all objects

## install packages 
install.packages(c('data.table', 'table1', 'lme4', 'ggplot2', 'jsonlite', 'minqa', 'numDeriv')) 
require(data.table)
require(table1) 
require(lme4)
require(ggplot2)
require(jsonlite)
require(minqa)
require(numDeriv)


## let's directly source pda functions from Github, please don't try to install from CRAN or Github
source('https://github.com/Penncil/pda/raw/master/R/pda.R')
source('https://github.com/Penncil/pda/raw/master/R/DLM.R')
source('https://github.com/Penncil/pda/raw/master/R/dlmm.R')
source('https://github.com/Penncil/pda/raw/master/R/ODAL.R')


## local working directory
# ***CHANGE THIS TO YOUR DIR***
setwd('/Users/chongliang/Dropbox/PDA-git/UNMC_workshop/run_examples/ota_cloud') 
# ***CHANGE THIS TO YOUR SITE***
mysite = 'MedicalCenter'      # 'MedicalCenter', 'Lincoln', 'Omaha', 'Kearney'



mydir = 'COVID_LOS_DLMM'      # my project working dir 
dir.create(mydir)             # create project working dir
site.name = c('MedicalCenter', 'Lincoln', 'Omaha', 'Kearney')
K = length(site.name)



## read in RDS data file, pooled data
url_mydata = 'https://github.com/chongliang-luo/UNMC_workshop/raw/main/data/COVID_LOS.rds' # github repo
COVID_LOS = readRDS(gzcon(url(url_mydata)))    # read in my data from github repo
 
########################### END: Setup ##############################


############################ Data summary ###########################
## descriptive
COVID_LOS   
table1(~LOS+Age+CCI+Gender+Admission+factor(Cancer)+factor(COPD)+factor(Hyperlipidemia)+factor(Hypertension)+
         factor(KidneyDisease)+factor(Obesity)+factor(HeartDisease)+factor(Diabetes)|site, data=COVID_LOS)
 
## Whiteboard Discussion...


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

## Whiteboard Discussion...

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
 

## DLMM non-interactive section: data communication at local computer
# STEP 0: remove any json files if exist
file.remove(list.files(mydir,full.names = T)[grepl('.json', list.files(mydir))])



# STEP 1: create control.json, read the R output
pda(site_id = control$lead_site, control=control, dir=mydir)



# STEP 2: calculate AD, read the R output  
pda(site_id ='Kearney', ipdata = COVID_LOS_split$Kearney, dir=mydir)
pda(site_id ='Omaha', ipdata = COVID_LOS_split$Omaha, dir=mydir)
pda(site_id ='Lincoln', ipdata = COVID_LOS_split$Lincoln, dir=mydir)
pda(site_id ='MedicalCenter', ipdata = COVID_LOS_split$MedicalCenter, dir=mydir)



# STEP 3: use AD files to fit DLMM, read the R output
pda(site_id = control$lead_site, dir=mydir, ipdata = COVID_LOS_split[[control$lead_site]] )

## END DLMM non-interactive section

##################### END: DLMM workflow #################################



 
 
############################   Conclusion    ###############################

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


## Let's also compare DLMM with other methods... 
## Let's fit LMM using the pooled data
fit.pooled <- lme4::lmer(LOS~Age+CCI+Gender+Admission+Cancer+COPD+Hyperlipidemia+Hypertension+
                           KidneyDisease + Obesity + HeartDisease + Diabetes +
                           (1|site), REML = F, data=COVID_LOS)
b_pool = round(summary(fit.pooled)$coef[,1], 4)
se_pool = round(summary(fit.pooled)$coef[,2], 4)

# fixed effects
data.frame(risk_factor=fit.dlmm$risk_factor, b_pool=b_pool, b_dlmm=fit.dlmm$bhat)

# var component
data.frame(RandomEffect=c(fit.dlmm$risk_factor_heterogeneity,  'random error'),
           VarComp.dlmm=data.frame(summary(fit.pooled)$varcor)$vcov, 
           VarComp.pool=c(fit.dlmm$Vhat, fit.dlmm$sigmahat^2) )

# random intercepts (BLUP) 
data.frame(site = control$sites, 
           BLUP.pool.Intercept = round(ranef(fit.pooled)$site$`(Intercept)`, 4),
           BLUP.dlmm.Intercept = fit.dlmm$uhat[,1])




## Let's try another naive federated analysis method: 
# By averaging individual estimates, this is called "meta-estimator" 

# site "Kearney" has no Q1, let's manually correct it (fill it with NA)
colnames(bi) = colnames(sei) = site.name
bi[,4] = c(bi[1:6,4],NA,bi[7:15,4]) 
sei[,4] = c(sei[1:6,4],NA,sei[7:15,4])

## "meta-estimator" (inverse-variance weighted average),
bmeta = round(rowSums(bi/sei^2,na.rm=T)/rowSums(1/sei^2,na.rm=T), 4)
semeta = round(sqrt(1/rowSums(1/sei^2,na.rm=T)), 4)

 


## Let's compare pooled vs meta vs DLMM, 
## with confidence interval (CI) plots of the fixed effects: 
Xname = names(b_pool)
px = length(Xname)
methods <- c('pooled', 'meta', 'DLMM')   
nm <- length(methods)
ci.df <- data.frame(method = rep(methods, each=px),       
                    risk.factor= rep(Xname, nm),
                    beta=as.numeric(c(b_pool, bmeta, fit.dlmm$bhat)),  
                    sd=as.numeric(c(se_pool, semeta, fit.dlmm$sebhat)),          
                    goldstandard=as.numeric(rep(b_pool, nm)))
ci.df$method <- factor(ci.df$method, levels = rev(methods)) 
ci.df$risk.factor <- factor(ci.df$risk.factor, levels = Xname)  
case = 'COVID LOS study'

CI_plot_COVID_LOS_DLMM <- 
  ggplot(ci.df, aes(x=method, y=beta, shape=method,color=method, group=risk.factor)) +
  geom_errorbar(aes(ymin=beta-1.96*sd, ymax=beta+1.96*sd), width=0.1) +
  geom_line() +
  geom_hline(aes(yintercept=goldstandard), linetype = "dashed") +
  geom_point(size=1.5) + 
  facet_wrap(. ~ risk.factor
             , scales= "free"  #  "fixed" # 
             , ncol=5)+
  labs(title= case, x =  '', y = "Effect estimate" ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)
        , axis.ticks.x = element_blank()
        , strip.text.x = element_text(size = 14)
        , panel.grid.major = element_blank() 
        , axis.title.y = element_text(size=9)
        , axis.text.x = element_text(size=9)
        , axis.text.y = element_text(size=9) 
        , legend.position = "inside",
        , legend.position.inside = c(0.7, 0.13)
        , legend.title = element_blank() 
        , legend.text=element_text(size=15) 
        , legend.key.height = unit(1.5, "line")
  ) + coord_flip()

CI_plot_COVID_LOS_DLMM

## save this plot as pdf file
ggsave("CI_plot_COVID_LOS_DLMM.pdf",  plot = CI_plot_COVID_LOS_DLMM, width =12, height =8)

## Lead site can upload this plot to OTA 

####################### END: Conclusion ############################


