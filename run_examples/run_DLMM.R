

COVID_LOS = readRDS('/Users/chongliang/Dropbox/PDA-git/UNMC_workshop/data/COVID_LOS.rds')
site.name = c('MedicalCenter', 'Lincoln', 'Omaha', 'Kearney')
head(COVID_LOS)  
table1(~LOS+Age+CCI+Gender+Admission+factor(Cancer)+factor(COPD)+factor(Hyperlipidemia)+factor(Hypertension)+
         factor(KidneyDisease)+factor(Obesity)+factor(HeartDisease)+factor(Diabetes)|site, data=COVID_LOS)

## Linear mixed-effects model, random intercepts
fit.pooled <- lme4::lmer(LOS~Age+CCI+Gender+Admission+Cancer+COPD+Hyperlipidemia+Hypertension+
                          KidneyDisease + Obesity + HeartDisease + Diabetes +
                          (1|site), REML = F, data=COVID_LOS)
bpool = round(summary(fit.pooled)$coef[,1], 4)
sepool = round(summary(fit.pooled)$coef[,2], 4)

# px = length(bpool)
# K = length(site.name)
# bi = sei = c() # matrix(NA, px, K)
# for(sid in site.name){
#   fit.i = coxph(Surv(time, status)~age+gender+RACE_NHW+smoking+CCI+depression+pain, data=OUD[site==sid,])
#   cat(sid, '\n')
#   print(round(summary(fit.i)$coef, 3))
#   bi = cbind(bi, summary(fit.i)$coef[,1])
#   sei = cbind(sei, summary(fit.i)$coef[,3])
# } 
# bmeta = round(rowSums(bi/sei^2,na.rm=T)/rowSums(1/sei^2,na.rm=T), 3)
# semeta = round(sqrt(1/rowSums(1/sei^2,na.rm=T)), 3)
# cbind(bpool, bmeta) 

COVID_LOS_split = split(COVID_LOS, by='site' )
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
res = run_DLMM_with_pda(control, mydir=getwd(), mydata=COVID_LOS_split, upload_without_confirm=T, silent_message=T)


# fixed effects
cbind(bpool, bdlmm=res$bhat, sepool, sedlmm=res$sebhat)

# var component
cbind(data.frame(summary(fit.pooled)$varcor)$vcov, 
      c(res$Vhat, res$sigmahat^2) )

# random effects (BLUP)
cbind(u.pool = ranef(fit.pooled)$site,
      u.dlmm = c(res$uhat))


## Chi-sq test for Gender random effect (likelihood ratio test, LRT)
COVID_LOS[,Gender:=ifelse(Gender=='F',0,1)] # create dummy var 
fit.pooled1 <- lme4::lmer(LOS~Age+CCI+Gender+Admission+Cancer+COPD+Hyperlipidemia+Hypertension+
                            KidneyDisease + Obesity + HeartDisease + Diabetes +
                            (1|site)+(0+Gender|site), REML = F, data=COVID_LOS)
anova(fit.pooled, fit.pooled1) 

# Chi-sq test for Gender random effect, via GLMM
control$variables_heterogeneity = c('Intercept', "Gender")
res1 = run_DLMM_with_pda(control, mydir=getwd(), mydata=COVID_LOS_split, upload_without_confirm=T, silent_message=T)
1 - pchisq((res1$lik - res0$lik)*2, df=1) # p-value

# var component
cbind(data.frame(summary(fit.pooled1)$varcor)$vcov, 
      c(diag(res1$Vhat), res1$sigmahat^2) )

# random effects (BLUP)
cbind(u.pool1 = ranef(fit.pooled1)$site,
      u.dlmm1.Intercept = res1$uhat[,1],
      u.dlmm1.Gender = res1$uhat[,2])
