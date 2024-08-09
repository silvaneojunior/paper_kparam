library(tidyverse)
library(mvtnorm)
library(latex2exp)
library(extrafont)

# Latest version on GitHub. Also available on CRAN
# devtools::install_github('silvaneojunior/kDGLM')
library(kDGLM)

#### Loading fonts ####
# Run only once
# font_import()

loadfonts(device = "win")
font_size=12
family_font='serif'
par(family=family_font)
#######################

Sys.setenv(R_GSCMD = "C:/Program Files/gs/gs10.01.1/bin/gswin64c.exe")
options(scipen = 999)

# Creates ribbons using base R functions
base_ribbon <- function(x,ymin,ymax,...){
  l=length(x)
  polygon(x=c(x[1],x,x[l:1]),
          y=c(ymax[1],ymin,ymax[l:1]),
          ...)
}
# Ribbon without fill
base_ribbon2 <- function(x,ymin,ymax,...){
  lines(x,ymax,lty=2,...)
  lines(x,ymin,lty=2,...)
}


# Visual options
prop.title=0.15
base.size=800
dpi=300