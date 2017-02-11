## ----eval = FALSE--------------------------------------------------------
## # Check to make sure the required packages are installed on your machine
## # If not, they will be installed
## 
## reqPackages <- c("devtools","raster","sp","maptools","geoshpere","rgeos","MASS")
## get.packages <- reqPackages[!(reqPackages %in% installed.packages()[,"Package"])]
## if(length(get.packages)>0) install.packages(get.packages, dependencies=TRUE)
## 
## # Install necessary packages from Github using the devtools library #
## library(devtools)
## install_github("SWotherspoon/SGAT")
## install_github("SLisovski/TwGeos")


.libPaths("Packages")
## ---- warning = FALSE, message = FALSE-----------------------------------
library(raster)
library(sp)
library(rgeos)
library(geosphere)
library(SGAT)
library(TwGeos)
library(MASS)
library(maptools)

## ------------------------------------------------------------------------
# read in a simple world map from the maptools package #
Americas<-raster::shapefile("Spatial_Layers/Americas.shp")
PETEdist<-raster::shapefile("Spatial_Layers/Sternula_lorata.shp")

WGS84<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

laea <- "+proj=laea +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

PETEdistB <- rgeos::gBuffer(spTransform(PETEdist,CRS(laea)),byid = TRUE, width = 1500000)
PETEdistBuffer <- spTransform(PETEdistB, CRS(WGS84))

## ------------------------------------------------------------------------
TernFiles <- list.files("Data",pattern = ".lig",full.names=TRUE)

# Read just the file names for Bird ID
BirdId <- list.files("Data", pattern=".lig")

# Determine the number of birds
nBirds<-length(BirdId)

## ------------------------------------------------------------------------

# Read all files in as Lig files #

Terndata <- lapply(TernFiles,readLig,skip = 1) 

# Keep only the relevant data 

for(i in 1:nBirds){
Terndata[[i]] <- Terndata[[i]][,c(2,4)]
}


## Subset the data to only include when on the birds ##
tern008 <- as.POSIXct(c("2015-09-11","2016-09-02"),format = "%Y-%m-%d",tz = "GMT")
tern009 <- as.POSIXct(c("2015-09-11","2016-09-03"),format = "%Y-%m-%d",tz = "GMT")
tern010 <- as.POSIXct(c("2015-09-12","2016-09-06"),format = "%Y-%m-%d",tz = "GMT")

Terndata[[1]] <- subset(Terndata[[1]],Date > tern008[1] & Date < tern008[2])
Terndata[[2]] <- subset(Terndata[[2]],Date > tern009[1] & Date < tern009[2])
Terndata[[3]] <- subset(Terndata[[3]],Date > tern010[1] & Date < tern010[2])

## ------------------------------------------------------------------------
# Set the capture coordinates for each bird #

CapLocs<-array(NA,c(nBirds,2))

CapLocs[1,] <- c(-70.357,-23.0544)    # 008
CapLocs[2,] <- c(-70.357,-23.0543)    # 009
CapLocs[3,] <- c(-70.359,-23.0474)    # 010

### ----------------------------------------------------------------------------------------- ### 

par(mfrow = c(3,1))
for(i in 1:nBirds){
lightImage(Terndata[[i]],
           offset = 19,
           zlim = c(0,5))
}

# Make empty list to store twilight events
twl <- vector('list',nBirds)

seed <- as.POSIXct(rep("2016-01-01 04:00:00",3),format = "%Y-%m-%d %H:%M:%S", tz = "GMT")

for(i in 1:nBirds){
twl[[i]] <- findTwilights(tagdata = Terndata[[i]],
                        threshold = 1,
                        include = seed[i],
                        dark.min = 240) # 0 hours minimum dark period
}


### ---- Edit twl ---- ###
twlEdit <- vector('list',nBirds)

for(i in 1:nBirds){
twlEdit[[i]] <- twilightEdit(twilights = twl[[i]], 
                    window = 4,           # two days before and two days after
                    outlier.mins = 35,    # difference in mins
                    stationary.mins = 25, # are the other surrounding twilights within 25 mins of one another
                    plot = TRUE)
}

## ----eval = FALSE--------------------------------------------------------

 for(i in 1:nBirds){
   twlEdit[[i]]<-twilightAdjust(twilights=twlEdit[[i]], interval=120)
 }


twlEdit <- lapply(twlEdit,subset,Deleted == FALSE)


## ----eval = FALSE--------------------------------------------------------
## for(i in 1:nBirds){
##   saveRDS(twlEdit[[i]],paste0("Twilights/",BirdId[[i]],".rds")) # Saves rds file with the BirdId.rds as the name
##     write.csv(twlEdit[[i]],paste0("Twilights/",BirdId[[i]],".csv"),row.names = FALSE)
## }


## ------------------------------------------------------------------------
# Create a vector with the dates known to be at deployment #
calibration.dates <- vector('list',nBirds)

for(i in 1:nBirds){
calibration.dates[[i]] <- c(Terndata[[i]][1,1],as.POSIXct("2015-09-18",tz="GMT"))
}

## ------------------------------------------------------------------------
# Extract twilight data during calibration period
calibration.data<-vector('list',nBirds)

for(i in 1:nBirds){
  calibration.data[[i]]<-subset(twlEdit[[i]],twlEdit[[i]]$Twilight>=calibration.dates[[i]][1] & 
                                             twlEdit[[i]]$Twilight<=calibration.dates[[i]][2])
}

## ------------------------------------------------------------------------
# create empty vectors to store data #
sun<-z<-zenith0<-zenith1<-twl_t<-twl_deviation<-alpha<-fitml<-vector("list",nBirds)

# loop through each of the two individuals #
for(i in 1:nBirds){
  
  # Calculate solar time from calibration data 
  sun[[i]]  <- solar(calibration.data[[i]][,1])
  
  # Adjust the solar zenith angle for atmospheric refraction
  z[[i]]<- refracted(zenith(sun = sun[[i]],
                            lon = CapLocs[i,1],
                            lat = CapLocs[i,2]))
  
  twl_t[[i]]   <- twilight(tm = calibration.data[[i]][,1], 
                           lon = CapLocs[i,1],
                           lat = CapLocs[i,2],
                           rise = calibration.data[[i]][,2],
                           zenith = quantile(z[[i]],probs=0.5))
  
  # Determine the difference in minutes from when the sun rose and the geolocator said it rose 
  twl_deviation[[i]] <- ifelse(calibration.data[[i]]$Rise, as.numeric(difftime(calibration.data[[i]][,1], twl_t[[i]], units = "mins")),
                         as.numeric(difftime(twl_t[[i]], calibration.data[[i]][,1], units = "mins")))
  
  twl_deviation[[i]]<-subset(twl_deviation[[i]],twl_deviation[[i]]>=0)
  
  # Describe the distribution of the error 
  fitml[[i]] <- fitdistr(twl_deviation[[i]], "log-Normal")
  # save the Twilight model parameters
  alpha[[i]] <- c(fitml[[i]]$estimate[1], fitml[[i]]$estimate[2]) 
}

## ----echo=FALSE----------------------------------------------------------
meanAlpha1<-mean(c(alpha[[1]][1],alpha[[2]][1],mean(alpha[[3]][1])))
meanAlpha2<-mean(c(alpha[[1]][2],alpha[[2]][2],mean(alpha[[3]][2])))

ALPHA<-alpha[[1]]
ALPHA[1]<-meanAlpha1
ALPHA[2]<-meanAlpha2

## ---- echo = FALSE, fig.cap = ("**Figure 4** *Left* The deviation in twilights from the true twilight at the capture site. *Right* The mean (point estimate) and 95% CI for the zenith angle at the capture location.")----

b<-unlist(twl_deviation)
cols<-c("red","blue","green","yellow","orange","purple","brown","gray","black","pink")
seq <- seq(0,60, length = 100)
par(mfrow=c(1,2),mar=c(4,4,0,0))
hist(b, freq = F,
     yaxt="n",
     ylim = c(0, 0.15),
     xlim = c(0, 60),
     breaks=15,
     col="gray",
     main = "",
     xlab = "Twilight error (mins)")
axis(2,las=2)
for(i in 1:nBirds){
lines(seq, dlnorm(seq, alpha[[i]][1], alpha[[i]][2]), col = cols[i], lwd = 3, lty = 2)
}

#Zenith angle plot
par(bty="l")
plot(median(z[[1]],na.rm=TRUE),xlim=c(1,10),ylim=c(85,100),pch=19,ylab="Zenith Angle",xlab="PETU",col=cols[1])
segments(1,quantile(z[[1]],probs=0.025),1,quantile(z[[1]],probs=0.975),col=cols[1])
for(i in 2:nBirds){
  par(new = TRUE)
  plot(median(z[[i]],na.rm=TRUE)~i,xlim=c(1,10),ylim=c(85,100),pch=19,yaxt="n",xaxt="n",ylab="",xlab="",col=cols[i])
  segments(i,quantile(z[[i]],probs=0.025),i,quantile(z[[i]],probs=0.975),col=cols[i])
}

## ------------------------------------------------------------------------
# Create empty vectors to store objects #
d.twl<-path<-vector('list',nBirds)

zenith0<-zenith1<-rep(NA,nBirds)

# loop through the birds #
for(i in 1:nBirds){
  # Store the zenith (sun-elevation angle)
  zenith0[i] <-quantile(z[[i]],prob=0.5)
  zenith1[i]<-quantile(z[[i]],prob=0.95)
}

## ----echo = FALSE--------------------------------------------------------
zenith1

tolvalues <- array(NA,c(nBirds,3))
tolvalues[,1]<-0.1
tolvalues[,2]<-0
tolvalues[,3]<-0.1


#Manual adjustments
#tolvalues[2,1]<-0.1
#tolvalues[3,1]<-0.2

#tolvalues[1,2]<-0.25
#tolvalues[2,2]<-0.25
#tolvalues[3,2]<-0.25

## ------------------------------------------------------------------------
path <- vector('list',nBirds)
for(i in 1:nBirds){  
   path[[i]] <- thresholdPath(twilight = twlEdit[[i]]$Twilight,
                             rise = twlEdit[[i]]$Rise,
                             zenith = zenith1[i],
                             tol = tolvalues[i,])

## ----echo = FALSE, fig.cap="**Figure 5** The initial annual cycle path of Olive-sided Flycatchers captured breeding in Alaska - *blue* = Fall, *green* = Spring, *red vertical lines* spring and fall equniox"----
  layout(matrix(c(1,3,
                  2,3), 2, 2, byrow = TRUE))
  par(mar=c(2,4,2,0))
  plot(path[[i]]$time, path[[i]]$x[, 2], type = "b", pch = 16, cex = 0.5, ylab = "Lat", xlab = '',xaxt="n")
  abline(h = CapLocs[i,2])
  abline(v = as.POSIXct("2014-09-23"),col="red",lty=2,lwd=1.5)
  abline(v = as.POSIXct("2015-03-20"),col="red",lty=2,lwd=1.5)
  par(mar=c(2,4,2,0))
  plot(path[[i]]$time, path[[i]]$x[, 1], type = "b", pch = 16, cex = 0.5, ylab = "Long", xlab = '')
  abline(h = CapLocs[i,1])
  abline(v = as.POSIXct("2014-09-23"),col="red",lty=2,lwd=1.5)
  abline(v = as.POSIXct("2015-03-20"),col="red",lty=2,lwd=1.5)
  
  
  plot(Americas, col = "grey95",xlim = c(-170,-60),ylim=c(-40,5))
  box()
  lines(path[[i]]$x, col = "blue")
  points(path[[i]]$x, pch = 16, cex = 0.5, col = "blue")
Sys.sleep(5)
}

## ------------------------------------------------------------------------
x0 <- z0 <- vector('list',nBirds)

for(i in 1:nBirds){
  # Take the location estimates created above
x0[[i]]<- path[[i]]$x

  # the model also needs the mid-points - generate those here
z0[[i]]<- trackMidpts(x0[[i]])
}

## ------------------------------------------------------------------------
beta <- c(0.7, 0.08)

## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
# set xlim and ylim values need to span the range of your dataset
xlim <- c(-100, -60)
ylim <- c(-90, 90)

## ------------------------------------------------------------------------
## Function to construct a land/sea mask
distribution.mask <- function(xlim, ylim, res = c(0.25,0.25), land = TRUE, shape) {
    r <- raster(res = res,
                xmn = xlim[1], 
                xmx = xlim[2], 
                ymn = ylim[1], 
                ymx = ylim[2], 
                crs = proj4string(shape))

    r <- cover(rasterize(shape, shift = c(-360, 0), r, 1, silent = TRUE), 
        rasterize(shape, r, 1, silent = TRUE), rasterize(elide(shape, 
           shift = c(360, 0)), r, 1, silent = TRUE))
    r <- as.matrix(is.na(r))[nrow(r):1, ]
   if (land) 
        r <- !r
   xbin <- seq(xlim[1], xlim[2], length = ncol(r) + 1)
    ybin <- seq(ylim[1], ylim[2], length = nrow(r) + 1)

    function(p) {
        r[cbind(.bincode(p[, 2], ybin), .bincode(p[, 1], xbin))]
    }
}

## ------------------------------------------------------------------------

crs(Americas)<-WGS84
crs(PETEdist)<-WGS84
crs(PETEdistBuffer) <- WGS84
## ------------------------------------------------------------------------
## Define mask for Ovenbird distribution
 is.dist <- distribution.mask(shape=PETEdistBuffer,
                             xlim = xlim,
                             ylim = ylim,
                             res = c(0.25,0.25), # 0.25 x 0.25 degree resolution
                             land = TRUE)

## ------------------------------------------------------------------------
# Define the log prior for x and z
log.prior <- function(p) {
    f <- is.dist(p)
    ifelse(f | is.na(f), 0, -10)
}

## ------------------------------------------------------------------------
# set time when birds are stationary 
fixedx <- vector('list',nBirds)
for(i in 1:nBirds){
fixedx[[i]] <- rep(FALSE, length(twlEdit[[i]]$Twilight))
fixedx[[i]][1:60]<-TRUE
fixedx[[i]][(length(fixedx[[i]])-60):length(fixedx[[i]])]<-TRUE
x0[[i]][fixedx[[i]],1] <- CapLocs[i,1]
x0[[i]][fixedx[[i]],2] <- CapLocs[i,2]
}
## ------------------------------------------------------------------------
# Define the threshold model - slimilar to above #
model <-  vector('list', nBirds)

for(i in 1:nBirds){
model[[i]]<- thresholdModel(twilight = twlEdit[[i]]$Twilight,
                            rise = twlEdit[[i]]$Rise,
                            twilight.model = "ModifiedLogNormal",
                            alpha = alpha[[i]],
                            beta = beta,
                            # Here is where we set the constraints for land
                            logp.x = log.prior, 
                            logp.z = log.prior, 
                            x0 = x0[[i]],
                            z0 = z0[[i]],
                            fixedx = fixedx[[i]],
                            zenith = zenith1[i])
}

## ------------------------------------------------------------------------
# This defines the error distribution around each location #
proposal.x <- proposal.z <- vector('list',nBirds)

for(i in 1:nBirds){
proposal.x[[i]] <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(x0[[i]]))
proposal.z[[i]] <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(z0[[i]]))
}

## ----eval = FALSE--------------------------------------------------------
 fit <- xsum <- zsum <- vector('list', nBirds)
 
 for(i in 1:nBirds){
 fit[[i]] <- estelleMetropolis(model = model[[i]],
                               proposal.x = proposal.x[[i]],
                               proposal.z = proposal.z[[i]],
                               iters = 2000, # This value sets the number of iterations to run
                               thin = 10,
                               chains = 3)

 xsum[[i]] <- locationSummary(fit[[i]]$x,collapse = TRUE)
 zsum[[i]] <- locationSummary(fit[[i]]$z,collapse = TRUE)

proposal.x[[i]] <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(cbind(xsum[[i]]$'Lon.50%',xsum[[i]]$'Lat.50%')))
proposal.z[[i]] <- mvnorm(S=diag(c(0.0025,0.0025)),n=nlocation(cbind(zsum[[i]]$'Lon.50%',zsum[[i]]$'Lat.50%')))

 fit[[i]] <- estelleMetropolis(model = model[[i]],
                              proposal.x = proposal.x[[i]],
                               proposal.z = proposal.z[[i]],
                               x0 = cbind(xsum[[i]]$'Lon.50%',xsum[[i]]$'Lat.50%'),
                               z0 = cbind(zsum[[i]]$'Lon.50%',zsum[[i]]$'Lat.50%'),
                               iters=2000, # This value sets the number of iterations to run
                               thin=10,
                               chains=3)
 
# Final Run
  xsum[[i]] <- locationSummary(fit[[i]]$x,collapse = TRUE)
  zsum[[i]] <- locationSummary(fit[[i]]$z,collapse = TRUE)

 proposal.x[[i]] <- mvnorm(chainCov(fit[[i]]$x),s=0.1)
 proposal.z[[i]] <- mvnorm(chainCov(fit[[i]]$z),s=0.1)
 
 # Note the increase in number of interations - this takes a bit longer to run
fit[[i]] <- estelleMetropolis(model = model[[i]],
                               proposal.x = proposal.x[[i]],
                               proposal.z = proposal.z[[i]],
                               x0=cbind(xsum[[i]]$'Lon.50%',xsum[[i]]$'Lat.50%'),
                               z0=cbind(zsum[[i]]$'Lon.50%',zsum[[i]]$'Lat.50%'),
                               iters=5000,  # This value sets the number of iterations to run
                               thin=10,
                               chains=3)
}
#
#


## ------------------------------------------------------------------------
# This step makes an empty raster #
r <- raster(res = c(0.25,0.25),
            xmn=xlim[1],
            xmx=xlim[2],
            ymn=-90,
            ymx=90)

## ------------------------------------------------------------------------

S <- Sp <- vector('list',nBirds)
for(i in 1:nBirds){
S[[i]] <- slices(type="intermediate",
                 weights = rep(0.5,length(fit[[i]][[1]]$time)),
                 breaks="day",
                 mcmc=fit[[i]],
                 grid=r,
                 include.lowest = FALSE, 
                 right = FALSE)
Sp[[i]] <- slices(type="primary",
                 weights = rep(0.5,length(fit[[i]][[1]]$time)),
                 breaks="day",
                 mcmc=fit[[i]],
                 grid=r,
                 include.lowest = FALSE, 
                 right = FALSE)
}


source("Functions/MigSchedules_copy.R")
Schedule <- vector('list',nBirds)
for(i in 1:nBirds){
Schedule[[i]] <- MigSchedule(MCMC = S[[i]], prob = 0.95,
                                    known.breed = c("2015-09-01","2015-10-01"),
                                    known.winter = c("2016-01-01","2016-02-01"),
						rm.lat.equinox = TRUE,
						days.omit = 10,
						progress = TRUE,
						plot = TRUE,
						plot.legend = FALSE,
						latAllow = c(0.5,2))
}


birds <- cols <- vector('list',nBirds)
cols[[1]] <- colorRampPalette(c("darkorchid1","darkorchid4"),alpha = T)
cols[[2]] <- colorRampPalette(c("tomato1","tomato3"),alpha = T)
cols[[3]] <- colorRampPalette(c("springgreen2","springgreen4"),alpha = T)
for(b in 1:nBirds){
KML(Schedule[[b]]$movements, subfolder.name = paste0(substr(BirdId[b],start = 1, stop = 16),".kml"), col = cols[[b]](100),overwrite = TRUE)
}
