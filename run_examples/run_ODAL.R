

setwd('/Users/chongliang/Dropbox/PDA-git/UNMC_workshop/test')
rm(list=ls())

require(data.table)
require(table1) 
require(ggplot2)
 

fetal_loss = readRDS('/Users/chongliang/Dropbox/PDA-git/UNMC_workshop/data/fetal_loss.rds')
site.name = c('MedicalCenter', 'Lincoln', 'Omaha', 'Kearney')
head(fetal_loss) 
table1(~factor(FetalLoss)+factor(MedX)+Race+Age+Weight+BMI|site, fetal_loss)

fit.pooled = glm(FetalLoss~MedX+Race+Age+Weight+BMI, family='binomial', data=fetal_loss)
round(summary(fit.pooled)$coef, 3)
bpool = round(summary(fit.pooled)$coef[,1], 4)
sepool = round(summary(fit.pooled)$coef[,2], 4)

px = length(bpool)
K = length(site.name)
bi = sei = c()  
for(sid in 1:K){
  fit.i = glm(FetalLoss~MedX+Race+Age+Weight+BMI, family='binomial', data=fetal_loss[site==site.name[sid],])
  print(round(summary(fit.i)$coef, 3))
  bi = cbind(bi, summary(fit.i)$coef[,1])
  sei = cbind(sei, summary(fit.i)$coef[,2])
}
bi[,4] = c(bi[1:3,4],NA,bi[4:7,4]) 
sei[,4] = c(sei[1:3,4],NA,sei[4:7,4]) 
bmeta = round(rowSums(bi/sei^2,na.rm=T)/rowSums(1/sei^2,na.rm=T), 4)
semeta = round(sqrt(1/rowSums(1/sei^2,na.rm=T)), 4)
cbind(bpool, bmeta) 
 
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

fetal_loss_split = split(fetal_loss,by='site')

res = run_ODAL_with_pda(control, mydir=getwd(), fetal_loss_split, upload_without_confirm=T, silent_message=T)




## CI plots
Xname = c("MedX", "RaceBlack", 'RaceAsian', 'RaceOther', "Age", "Weight", "BMI")
px = length(Xname)
methods <- c('pooled', 'meta', 'ODAL')   
nm <- length(methods)
ci.df <- data.frame(method = rep(methods, each=px),       
                    risk.factor= rep(Xname, nm),
                    beta=as.numeric(c(bpool[-1], bmeta[-1], res$btilde[-1])),  
                    sd=as.numeric(c(sepool[-1], semeta[-1], res$setilde[-1])),          
                    goldstandard=as.numeric(rep(bpool[-1], nm)))
ci.df$method <- factor(ci.df$method, levels = rev(methods)) 
ci.df$risk.factor <- factor(ci.df$risk.factor, levels = Xname)  
case = 'Fetal loss study'

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
 