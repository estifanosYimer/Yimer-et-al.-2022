rm(list = ls())

# to get where this script is saved in the local disk
list.of.packages <- c("rstudioapi","reticulate")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(rstudioapi)
library(reticulate)

working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

setwd(paste0(working_dir,'/','Precipitation and PET data'))

source("Chirps data and bias correction.R")

setwd(working_dir)
source("PET estimation.R")

working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste0(working_dir,'/','GOF'))

source("GOF fitdist AD.R")

working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste0(working_dir,'/','GOF'))
source("GOF fitdist KS.R")

working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste0(working_dir,'/','SPEI'))

source("SPEI.R")

# now we can run the SW test as it is on the estimated SPEI values contrary to the classical GOF test using AD and KS where they assess the fit of distribution
# to a water balance values
working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste0(working_dir,'/','GOF'))

source("SW test.R")

working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste0(working_dir,'/','Acceptance frequency each month'))

source("GOF fitdist AD each month.R")

working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste0(working_dir,'/','SPAEI'))

source("GLEAM AET.R")

working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste0(working_dir,'/','SPAEI'))

source("GOF fitdist AD SPAEI.R")

working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(paste0(working_dir,'/','SPAEI'))

source("SPAEI.R")





