####################### Load the data #############################################################
# to get where this script is saved in the local disk
rm(list = ls())

list.of.packages <- c("rstudioapi","SPEI","SCI","")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(rstudioapi)

working_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

setwd(working_dir)

library(SPEI)
library(SCI)
library(lubridate)
library(lattice)
library(RColorBrewer)
library(ncdf4.helpers)
library(PCICt)
library(ncdf4)
library(sf)
library(terra)
require("evd")

################## trial of reading a netcdf file and converting it to an array ###################

setwd('change your directory to your local one') # change your directory

name <- "study_area1.nc"

ncin <- nc_open(name)
print(ncin)

lon <- ncvar_get(ncin,"lon")
nlon <- dim(lon)
head(lon)

lat <- ncvar_get(ncin,"lat")
nlat <- dim(lat)
head(lat)

print(c(nlon,nlat))

##################### getting the time as timeseries #########

variables = names(ncin[['var']])

#tas_time <- nc.get.time.series(ncin, v = "Band1",
#                               time.dim.name = "time")

tas_time <- format(seq(as.Date("1981-01-13"), as.Date("2021-12-13"), by="month"),format="%m-%Y") #define your time

dname <- "Band1"

tmp_array <- ncvar_get(ncin,dname)
dlname <- ncatt_get(ncin,dname,"long_name")
dunits <- ncatt_get(ncin,dname,"units")
fillvalue <- ncatt_get(ncin,dname,"_FillValue")
dim(tmp_array)

tmp_array[tmp_array==fillvalue$value] <- NA # replacing missing values

dimnames(tmp_array)<- list(lon,lat,tas_time)

save(tmp_array,file='study area1.RData')

#####################################################################################
############### Part 2: Drought analysis using SCI or SPEI package ##################
#####################################################################################


acc.lst= c("1","6","9")  #accumulation period

amt.acc <- length(acc.lst)

date.lst <- dimnames(tmp_array)[[3]]

amt.date <- length(date.lst)

#yr.lst <- year(date.lst)

lon.lst <- dimnames(tmp_array)[[1]]

amt.lon <- length(lon.lst)

lat.lst <- dimnames(tmp_array)[[2]]

amt.lat <- length(lat.lst)

SPI.value <- array(NA, 
                           dim = c(lon = amt.lon, lat = amt.lat,date = amt.date,accumulation=amt.acc), 
                           dimnames = list(lon = lon.lst,lat = lat.lst, date = date.lst,accumulation=acc.lst)) 


S.value <- array(NA, 
                  dim = c(lon = amt.lon, lat = amt.lat,accumulation=amt.acc,S=20), 
                  dimnames = list(lon = lon.lst,lat = lat.lst,accumulation=acc.lst,S=1:20)) 

D.value <- array(NA, 
                 dim = c(lon = amt.lon, lat = amt.lat,accumulation=amt.acc,D=20), 
                 dimnames = list(lon = lon.lst,lat = lat.lst,accumulation=acc.lst,D=1:20)) 

for(i.lon in seq(1, amt.lon)){
  for(i.lat in seq(1, amt.lat)){
    
    #print(i.lon)
    
    time.arr.to.use.SM <- tmp_array[i.lon,i.lat,]
    
    for (spi.month.to.use in acc.lst){
      
      tryCatch({
        
        time.arr.to.use.SM <- as.vector(t(time.arr.to.use.SM))
        S.value <- spi(time.arr.to.use.SM, spi.month.to.use)$fitted
        
        # if you want to use the SCI package, use the next two lines
        
        #spi.para <- fitSCI(time.arr.to.use,first.mon=1,time.scale=spi.month.to.use,distr="gamma",p0=TRUE)
        
        #S.value <- transformSCI(time.arr.to.use,first.mon=1,obj=spi.para)
        
        SPI.value[i.lon,i.lat,,spi.month.to.use]<-S.value
        
        #Drought duration and severity
        tryCatch({
          thres<-1  ### Drought defined as SPI<= -1
          CL<-clusters(-S.value,u=thres,r=2)  ## Cluster of SPI-values below -1
          n<-length(CL) ## number of clusters
        
          ##################################
          SD<-matrix(0,n,3)  
          colnames(SD)<-c("PEAK","SEVERITY","DURATION")
          rownames(SD)<-1:n
          ####################################################
          ### Loop over clusters
          ### (Year: year where the cluster maximum is reached)
          ####################################################
          
          # year list 
          
          yr.lst.na.rm <- seq(1981,2021,1) 
          
          year <- c()
          for (y in yr.lst.na.rm){
            year <- append(year,as.numeric(rep(y,12)))
          }
          
          
          i<-1
          for (cl in CL) {
            rownames(SD)[i]<-year[as.integer(names(cl)[max(cl)==cl])]
            SD[i,]<-c(max(cl),sum(cl),length(cl))
            i<-i+1
          }
          
          # severity
          S.value[i.lon,i.lat,spi.month.to.use,1:n]<- SD[1:n,2]
          
          # duration
          D.value[i.lon,i.lat,spi.month.to.use,1:n]<- SD[1:n,3]
        
        },error=function(e){
          
          print(" ")
        })
        
      },error=function(e){
        
        print(" ")
      })
      
      
      
      
      
    }
    
  }
  
}

save(SPI.value,file='SPI_values.RData')

save(S.value,file='SPI_severity.RData')

save(D.value,file='SPI_duration.RData')

