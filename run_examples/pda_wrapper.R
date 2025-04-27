
# pda_wrapper: ODAL

source("~/Dropbox/PDA-git/pda/R/ODAL.R")

Rcpp::sourceCpp("/Users/chongliang/Dropbox/PDA-git/pda/src/rcpp_coxph.cpp")
source("/Users/chongliang/Dropbox/PDA-git/pda/R/ODAC.R")

source("/Users/chongliang/Dropbox/PDA-git/pda/R/dlmm.R")
source("/Users/chongliang/Dropbox/PDA-git/pda/R/DLM.R")

source("~/Dropbox/PDA-git/pda/R/pda.R")

## run ODACH or ODACH_CC with pda silently at local 
run_ODAL_with_pda <- function(control, mydir, mydata, upload_without_confirm=T, silent_message=T){
  file.remove(list.files(mydir)[grepl('.json', list.files(mydir))])
  
  K = length(control$sites)
  # px = length(control$variables) 
  
  # ############################  STEP 1: initialize  ###############################
  ## assume lead site1: enter "1" to allow transferring the control file
  pda(site_id = control$lead_site, control = control, dir = mydir, 
      upload_without_confirm= upload_without_confirm, silent_message=silent_message)
  for(sid in K:1) pda(site_id = control$sites[sid], ipdata = mydata[[sid]], dir=mydir, 
                      upload_without_confirm =upload_without_confirm, silent_message=silent_message) 
  pda(site_id = control$lead_site, ipdata = mydata[[control$lead_site]], dir=mydir, 
      upload_without_confirm =upload_without_confirm, silent_message=silent_message) 
  
  # ############################  STEP 2: derivative  ###############################
  for(sid in K:1) pda(site_id = control$sites[sid], ipdata = mydata[[sid]], dir=mydir, 
                      upload_without_confirm = upload_without_confirm, silent_message=silent_message) 
  
  # ############################  STEP 3: estimate  ###############################
  tt=tryCatch(pda(site_id = control$lead_site, ipdata = mydata[[control$lead_site]], dir=mydir,
                  upload_without_confirm = upload_without_confirm, silent_message=silent_message),error=function(e) NA)
  
  config <- getCloudConfig(site_id =control$lead_site, dir=mydir) 
  
  if(!is.na(tt[1])){
    fit.pda <- pdaGet(paste0(control$lead_site,'_estimate'), config = config)
    return(list(btilde = fit.pda$btilde, setilde=fit.pda$setilde)) # sqrt(diag(solve(fit.pda$Htilde))/nrow(dd))
  }else{
    return(list(btilde = NA, setilde=NA))
  }
  
}
 

## run ODACH or ODACH_CC with pda silently at local 
run_ODACH_with_pda <- function(control, mydir, mydata, upload_without_confirm=T, silent_message=T){
  file.remove(list.files(mydir)[grepl('.json', list.files(mydir))])
  
  K = length(control$sites)
  # px = length(control$variables)
  # line #180 in ODAC.R 
  # for(sid in 1:length(mydata)) mydata[[sid]][,1]=round(mydata[[sid]][,1],4)
  
  # ############################  STEP 1: initialize  ###############################
  ## assume lead site1: enter "1" to allow transferring the control file
  pda(site_id = control$lead_site, control = control, dir = mydir, 
      upload_without_confirm= upload_without_confirm, silent_message=silent_message)
  for(sid in K:1) pda(site_id = control$sites[sid], ipdata = mydata[[sid]], dir=mydir, 
                      upload_without_confirm =upload_without_confirm, silent_message=silent_message)
  pda(site_id = control$lead_site, ipdata = mydata[[control$lead_site]], dir=mydir, 
      upload_without_confirm =upload_without_confirm, silent_message=silent_message)
  
  # ############################  STEP 3: derivative  ###############################
  for(sid in K:1) pda(site_id = control$sites[sid], ipdata = mydata[[sid]], dir=mydir, 
                      upload_without_confirm = upload_without_confirm, silent_message=silent_message)
  
  # ############################  STEP 4: estimate  ###############################
  tt=tryCatch(pda(site_id = control$lead_site, ipdata = mydata[[1]], dir=mydir, 
                  upload_without_confirm = upload_without_confirm, silent_message=silent_message),error=function(e) NA)
  
  config <- getCloudConfig(site_id =control$lead_site, dir=mydir)
  if(!is.na(tt[1])){
    fit.pda <- pdaGet(paste0(control$lead_site,'_estimate'), config = config)
    return(list(btilde = fit.pda$btilde, setilde=fit.pda$setilde)) # sqrt(diag(solve(fit.pda$Htilde))[31]/nrow(dd))
  }else{
    return(list(btilde = NA, setilde=NA))
  }
  
}






## run DLMM with pda silently at local 
run_DLMM_with_pda <- function(control, mydir, mydata, upload_without_confirm=T, silent_message=T){
  file.remove(list.files(mydir)[grepl('.json', list.files(mydir))])
  K = length(control$sites)
  
  pda(site_id = control$lead_site, control=control, dir=mydir,
      upload_without_confirm =upload_without_confirm, silent_message=silent_message)
  
  for(sid in K:1) pda(site_id = control$sites[sid], ipdata = mydata[[sid]], dir=mydir, 
                      upload_without_confirm =upload_without_confirm, silent_message=silent_message)

  pda(site_id = control$lead_site, dir=mydir, ipdata = mydata[[sid]], 
      upload_without_confirm =upload_without_confirm, silent_message=silent_message)
  
  config <- getCloudConfig(site_id =control$lead_site, dir=mydir)
  fit.dlmm <- pdaGet(name = paste0(control$lead_site,'_estimate'), config = config)
  return(fit.dlmm) 
}
