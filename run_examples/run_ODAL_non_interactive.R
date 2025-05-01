################################################################################
################   UNMC workshop: PDA-OTA practice #2 ODAL  ####################
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



mydir = 'fetal_loss_ODAL'       # my project working dir 
dir.create(mydir)               # create project working dir
site.name = c('MedicalCenter', 'Lincoln', 'Omaha', 'Kearney')
K = length(site.name)




## read in RDS data file 
url_mydata = 'https://github.com/chongliang-luo/UNMC_workshop/raw/main/data/fetal_loss.rds' # github repo
fetal_loss = readRDS(gzcon(url(url_mydata)))    # read in my data from github repo

########################### END: Setup ##############################



############################ Data summary ###########################
## descriptive
fetal_loss   
table1(~factor(FetalLoss)+factor(MedX)+Race+Age+Weight+BMI|site, fetal_loss )

## Whiteboard Discussion...


## Logistic regression at each site,  
bi = sei = c()  
for(sid in site.name){
  fit.i = glm(FetalLoss~MedX+Race+Age+Weight+BMI, family='binomial', data=fetal_loss[site==sid,])
  cat(sid,'\n')
  print(round(summary(fit.i)$coef, 3))
  bi = cbind(bi, summary(fit.i)$coef[,1])
  sei = cbind(sei, summary(fit.i)$coef[,2])
}

## Whiteboard Discussion...
 
######################### END: Data summary ###########################



####################### ODAL workflow ##################################
## split pooled data by sites
fetal_loss_split = split(fetal_loss, by='site')

## setup control
control <- list(project_name = 'Medication use and fetal loss (ODAL)',
                sites = c('MedicalCenter', 'Lincoln', 'Omaha', 'Kearney'),
                lead_site = 'MedicalCenter', 
                step = 'initialize',
                model = 'ODAL',
                family = 'binomial',
                heterogeneity = FALSE,
                outcome = 'FetalLoss',
                variables = c('MedX', 'Race', 'Age', 'Weight', 'BMI'),
                variables_lev = list(Race=c('White', 'Black', 'Asian', 'Other')),
                init_method = "meta",   
                optim_maxit = 100,
                optim_method = 'BFGS',
                upload_date = as.character(Sys.time()) )
 


## ODAL non-interactive section: data communication at local computer
# STEP 0: remove any json files if exist
file.remove(list.files(mydir,full.names = T)[grepl('.json', list.files(mydir))])


# STEP 1: create control.json, read the R output
pda(site_id = control$lead_site, control=control, dir=mydir)


# STEP 2: calculate AD (individual estimates), read the R output  
pda(site_id ='Kearney', ipdata = fetal_loss_split$Kearney, dir=mydir)
pda(site_id ='Omaha', ipdata = fetal_loss_split$Omaha, dir=mydir)
pda(site_id ='Lincoln', ipdata = fetal_loss_split$Lincoln, dir=mydir)
pda(site_id ='MedicalCenter', ipdata = fetal_loss_split$MedicalCenter, dir=mydir)
# pda(site_id = control$lead_site, ipdata = fetal_loss_split[[control$lead_site]], dir=mydir)
# STEP 3: control file is automatically updated... 



# STEP 4: calculate AD (likelihood derivatives), read the R output 
pda(site_id ='Kearney', ipdata = fetal_loss_split$Kearney, dir=mydir)
pda(site_id ='Omaha', ipdata = fetal_loss_split$Omaha, dir=mydir)
pda(site_id ='Lincoln', ipdata = fetal_loss_split$Lincoln, dir=mydir)
pda(site_id ='MedicalCenter', ipdata = fetal_loss_split$MedicalCenter, dir=mydir)
# pda(site_id = control$lead_site, ipdata = fetal_loss_split[[control$lead_site]], dir=mydir)



# STEP 5: calculate ODAL estimate, read the R output  
pda(site_id = control$lead_site, ipdata = fetal_loss_split[[control$lead_site]], dir=mydir)


## END ODAL non-interactive section 

##################### END: ODAL workflow #################################





###########################    Conclusion    #############################

## OK pda-ODAL is completed now, let's check the results
# read in the ODAL results
config <- getCloudConfig(site_id=mysite, dir=mydir)
fit.odal <- pdaGet(name = paste0(control$lead_site,'_estimate'), config = config)
fit.odal
 

## Let's also compare ODAL with other methods...
## Let's fit Logistic regression using pooled data 
fit.pooled = glm(FetalLoss~MedX+Race+Age+Weight+BMI, family='binomial', data=fetal_loss)
round(summary(fit.pooled)$coef, 3)
bpool = round(summary(fit.pooled)$coef[,1], 4)
sepool = round(summary(fit.pooled)$coef[,2], 4)


## Let's try another naive federated analysis method: 
# By averaging individual estimates, this is called "meta-estimator" 
# first fit a logistic regression model at each site

# site "Kearney" has no coef for Asian, let's manually correct it (fill it with NA)
bi[,4] = c(bi[1:3,4],NA,bi[4:7,4]) 
sei[,4] = c(sei[1:3,4],NA,sei[4:7,4]) 

## "meta-estimator" (inverse-variance weighted average),
bmeta = round(rowSums(bi/sei^2,na.rm=T)/rowSums(1/sei^2,na.rm=T), 4)
semeta = round(sqrt(1/rowSums(1/sei^2,na.rm=T)), 4)




## Let's compare pooled vs meta vs DLMM, 
## with confidence interval (CI) plots of the effect estimates: 
Xname = c("MedX", "RaceBlack", 'RaceAsian', 'RaceOther', "Age", "Weight", "BMI")
px = length(Xname)
methods <- c('pooled', 'meta', 'ODAL')   
nm <- length(methods)
ci.df <- data.frame(method = rep(methods, each=px),       
                    risk.factor= rep(Xname, nm),
                    beta=as.numeric(c(bpool[-1], bmeta[-1], fit.odal$btilde[-1])),  
                    sd=as.numeric(c(sepool[-1], semeta[-1], fit.odal$setilde[-1])),          
                    goldstandard=as.numeric(rep(bpool[-1], nm)))
ci.df$method <- factor(ci.df$method, levels = rev(methods)) 
ci.df$risk.factor <- factor(ci.df$risk.factor, levels = Xname)  
case = 'Fetal loss study'

CI_plot_fetal_loss_ODAL <- 
  ggplot(ci.df, aes(x=method, y=beta, shape=method,color=method, group=risk.factor)) +
  geom_errorbar(aes(ymin=beta-1.96*sd, ymax=beta+1.96*sd), width=0.1) +
  geom_line() +
  geom_hline(aes(yintercept=goldstandard), linetype = "dashed") +
  geom_point(size=1.5) + 
  facet_wrap(. ~ risk.factor
             , scales= "free"  #  "fixed" # 
             , ncol=3)+
  labs(title= case, x =  '', y = "Effect estimate (log odds ratio)" ) +
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

CI_plot_fetal_loss_ODAL

## save this plot as pdf file
ggsave("CI_plot_fetal_loss_ODAL.pdf",  plot = CI_plot_fetal_loss_ODAL, width =10, height =8)

## Lead site can upload this plot to OTA 

####################### END: Logistic reg Pooled data ##########################



