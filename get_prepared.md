UNMC workshop: PDA-OTA practice software setup
==============================================
   
  
1. Install R (v4.4.3 or later) and Rstudio:
  - https://posit.co/download/rstudio-desktop/
  
2. Install R packages in Rstudio using the script:
```r
install.packages(c('data.table', 'table1', 'lme4', 'ggplot2')) 
require(data.table)
require(table1) 
require(lme4)
require(ggplot2)
```

3. Download practice R code to your folder (Code -> Download ZIP): 
- https://github.com/chongliang-luo/UNMC_workshop/
  
4. Source pda code. 
For this workshop we will use the pda package on Github. Please directly source the code from Github, no need to install the package.
```r
source('https://github.com/Penncil/pda/raw/master/R/pda.R')
source('https://github.com/Penncil/pda/raw/master/R/DLM.R')
source('https://github.com/Penncil/pda/raw/master/R/dlmm.R')
source('https://github.com/Penncil/pda/raw/master/R/ODAL.R')
```
 
