################################################################################
################   UNMC workshop: PDA-OTA practice #2 ODAL  ####################
################################################################################
 

############################## Setup ##############################
# rm(list=ls()) # empty R， remove all objects

## install packages
require(data.table)
require(table1) 
require(ggplot2)
require(remotes)
install_github('https://github.com/Penncil/pda')
require(pda)


## local directory
setwd('/Users/chongliang/Dropbox/PDA-git/UNMC_workshop/run_examples/ota_cloud') # my working dir
mydir = 'fetal_loss_ODAL'     # my project working dir 
dir.create(mydir)             # create project working dir
mysite = 'MedicalCenter'      # my site name 'MedicalCenter', 'Lincoln', 'Omaha', 'Kearney'

## read in  csv data file 
# url_mydata = paste0('https://github.com/chongliang-luo/UNMC_workshop/raw/main/data/fetal_loss_', mysite, '.rds') # github repo
# COVID_LOS_mydata = readRDS(gzcon(url(url_mydata)))    # read in my data from github repo

url_mydata = 'https://github.com/chongliang-luo/UNMC_workshop/raw/main/data/fetal_loss.csv' # github repo
fetal_loss = fread( url_mydata )
fetal_loss[,Race:=factor(Race, levels = c('White', 'Black', 'Asian', 'Other'))]
########################### END: Setup ##############################



############################ Data summary ###########################
## descriptive
head(fetal_loss)   
table1(~factor(FetalLoss)+factor(MedX)+Race+Age+Weight+BMI|site, fetal_loss )
site.name = c('MedicalCenter', 'Lincoln', 'Omaha', 'Kearney')
K = length(site.name)

## Logistic regression at each site,  
bi = sei = c()  
for(sid in site.name){
  fit.i = glm(FetalLoss~MedX+Race+Age+Weight+BMI, family='binomial', data=fetal_loss[site==sid,])
  cat(sid,'\n')
  print(round(summary(fit.i)$coef, 3))
  bi = cbind(bi, summary(fit.i)$coef[,1])
  sei = cbind(sei, summary(fit.i)$coef[,2])
}
# site "Kearney" has no coef for Asian, let's manually correct it (fill it with NA)
bi[,4] = c(bi[1:3,4],NA,bi[4:7,4]) 
sei[,4] = c(sei[1:3,4],NA,sei[4:7,4]) 

## Average them with inverse-variance as weights. This is called "meta-estimator" 
bmeta = round(rowSums(bi/sei^2,na.rm=T)/rowSums(1/sei^2,na.rm=T), 4)
semeta = round(sqrt(1/rowSums(1/sei^2,na.rm=T)), 4)

## Discussion...
######################### END: Data summary ###########################


####################### ODAL workflow ##################################
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
fetal_loss_split = split(fetal_loss, by='site')
 
## ODAL interactive section:  data communication at local computer
# [all sites]: remove any json files if exist
file.remove(list.files(mydir,full.names = T)[grepl('.json', list.files(mydir))])


# [lead site]: create control file, 
pda(site_id = control$lead_site, control=control, dir=mydir)


# [all sites]: initialize step, calculate individual estimates,   
for(sid in site.name) pda(site_id = sid, ipdata = fetal_loss_split[[sid]], dir=mydir, 
                    upload_without_confirm =T, silent_message=T)

  
# [all sites]: derive step, calculate derivatives, 
for(sid in site.name) pda(site_id = sid, ipdata = fetal_loss_split[[sid]], dir=mydir, 
                    upload_without_confirm =T, silent_message=T)
 

# [lead site]: estimate step, calculate ODAL estimate, 
pda(site_id = mysite, dir=mydir, ipdata = fetal_loss_split[[mysite]])
## END interactive section 


## OK pda-ODAL is completed now, let's check the results
# read in the ODAL results
config <- getCloudConfig(site_id=mysite, dir=mydir)
fit.odal <- pdaGet(name = paste0(control$lead_site,'_estimate'), config = config)
fit.odal

##################### END: ODAL workflow #################################






########################## Logistic reg Pooled data ############################

## Logistic regression using pooled data 
fit.pooled = glm(FetalLoss~MedX+Race+Age+Weight+BMI, family='binomial', data=fetal_loss)
round(summary(fit.pooled)$coef, 3)
bpool = round(summary(fit.pooled)$coef[,1], 4)
sepool = round(summary(fit.pooled)$coef[,2], 4)


## make CI plots: compare pooled vs meta vs ODAL
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
  labs(title= case, x =  '', y = "Effect size (log odds ratio)" ) +
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

ggsave("CI_plot_fetal_loss_ODAL.pdf",  plot = CI_plot_fetal_loss_ODAL, width =10, height =8)


####################### END: Logistic reg Pooled data ##########################



