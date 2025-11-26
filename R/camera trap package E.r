# Estimate population density using camera trapping data
#============================================================================================================================
# library(cameratrapR)


# plot camera trapping results
#############################################################################################################################
#' plot camera trapping results, showing the locations of the cameras and the number of pictures
#'
#' @description This function shows the locations of each camera in a camera trap grid, the number of pictures
#'  taken by each camera for one species.
#'
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param x A data.frame with column names "Lon", "Lat", "Group_size", "Date", "Time"
#'
#' @param circle.size A value controls the size the circles. The defualt value is 0.2.
#'
#' @param point.scatter A value controls the distance between every point (representing a picture)
#'  at the same camera. The default value is 5.
#'
#' @examples
#'
#'  attach(trapresult)
#'  plotCamtrap(trapresult)
#'  plotCamtrap(trapresult, circle.size = .4, point.scatter = 5)
#'
#' @export

plotCamtrap = function(x, circle.size=.2, point.scatter=5){
  list.spe = c('Cape hare','Red deer','Roe deer','Moose','Chinese goral','Wild boar','Red fox',
               'Corsac fox','Raccoon dog','Eurasian lynx','Leopard cat',
               'Yellow-throated marten','Sable','Siberian weasel')
  if(!x$Species[1] %in% list.spe) print(paste(x$Species[1], "is not in the list. Movement parameters are not available", sep=""))

  sum = aggregate(x$Group_size, by = list(x$Lat, x$Lon), sum)
  colnames(sum) = c('Lat','Lon','Count')
  plot(x$Lon, x$Lat, pch=1, cex = 1, xlab='Longitude', ylab='Latitude')
  # points(sum$Lon, sum$Lat, pch = 16,
  #       col = colorRampPalette(c("grey90", "grey50"))(length(sum$Count))[round(rank(sum$Count))], cex=sum$Count*circle.size)
  points(sum$Lon, sum$Lat, pch = 16,
         col = adjustcolor("grey10", alpha.f = 0.3), cex=sum$Count*circle.size)
  X = x[x$Group_size > 0, ]
  dates = as.numeric(as.Date(X$Date, origin = "1900-01-01"))
  points(jitter(X$Lon, factor=point.scatter), jitter(X$Lat, factor=point.scatter),
         pch = 1, cex=1, col = colorRampPalette(c("red", "yellow", "green"))(length(dates))[round(rank(dates))])
}
#############################################################################################################################





# AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
#' plot the daily rhythm of species' activity
#'
#' @description This function plots the curves of probability density of animal activity (density of picture time) in 24 hours
#'
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param x A data.frame with column names "Lon", "Lat", "Group_size", "Date", "Time"
#'
#' @examples
#'
#'  attach(trapresult)
#'  dailyRhythm(trapresult)
#'
#' @export

dailyRhythm = function(x){
  times <- x$Time
  baseline = as.numeric(strptime('0:0:0', "%H:%M:%S"))
  times <- (as.numeric(strptime(times, "%H:%M:%S")) - baseline) / 3600
  plot(density(times, bw=0.3), xlim=c(0,24), xlab='Hour', ylab = 'Activity', main='')
  lines(density(times, bw=2),lwd = 2)
}
# AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA




# CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#' Convert longitude & latitude from geographic to projected coordinate system
#'
#' @description This function converts longitude & latitude from geographic to projected coordinate system.
#'  Users can specify their coordinate reference systems, or leave them as default. The default setting will
#'  using "WGS83" as geographic coordinates system and "UTM" as the projection coordinates system. The "zone"
#'  for UTM projection system will be determined using the centriod of points (longitude & latitude).
#'
#' @author Huidong Tian (tienhuitung@gmail.com)
#'
#' @param data data.frame with column names "Lon" & "Lat"
#'
#' @param geo geographic coordinate system, e.g. "+proj=longlat +datum=WGS84"
#'
#' @param proj projection coordinate system, e.g. "+proj=utm +zone=30 ellps=WGS84"
#'
#' @return Return a data.frame with column names "Lon" & "Lat". The unit of "Lon" & "Lat" is meter now.
#'
#' @examples
#'
#'  camera.geo <- data.frame(Lon = rnorm(10, 120, 10), Lat = rnorm(10, 30, 5))
#'  par(mfrow = c(1, 2))
#'  with(camera.geo, plot(Lon, Lat, xlab = "Longitude", ylab = "Latitude", main = "Geographic coordinate system"))
#'  camera.proj <- geo2proj(camera.geo)
#'  with(camera.proj, plot(Lon, Lat, xlab = "Longitude", ylab = "Latitude", main = "Universal Transverse Mercator (UTM) \n Coordinate System"))
#'
#' @import sp
#' @export
#'

geo2proj <- function(data, geo = NULL, proj = NULL) {

  if (!requireNamespace("sp", quietly = TRUE)) {
    stop("Package 'sp' needed for this function to work. Please install it.",
         call. = FALSE)
  }

  if (max(abs(data$Lat)) > 90 | max(abs(data$Lon)) >180) {
    stop("The longitude and/or latitude is not in geographic coordinate system!")
  }


  coordinates(data) <- c("Lon", "Lat")
  ## Assign geographic CRS
  if (!is.null(geo)) {
    proj4string(data) <- CRS(geo)
  } else {
    proj4string(data) <- CRS("+proj=longlat +datum=WGS84")
  }
  ## Convert to projection CRS
  if (!is.null(proj)) {
    res <- spTransform(data, CRS(proj))
  } else {
    zone <- floor((mean(data$Lon) + 180)/6)%%60 + 1
    res <- spTransform(data, CRS(sprintf("+proj=utm +zone=%s ellps=WGS84", zone)))
    res <- as.data.frame(res)
    names(res)[c(ncol(res)-1, ncol(res))] <- c("Lon","Lat")
  }
  return(as.data.frame(res))
}
# CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

# library(sp)
# attach(trapresult)
# trapresult = geo2proj(trapresult) # converts longitude & latitude from geographic coordinate system to equal area UTM system









# simulate camera trapping processes
#SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
#' Simulate animal movement within the range of the camera trap grid and obtain the pseudo camera trap result
#'
#' @description Correlated random walk of the target animal is simulated within the range of camera trap grid,
#'  using the distributions of step length, turning angles, and size of home range from footprint chain data.
#'  The simulated movement of the default 1-10 individuals generate pseudo camera trap data, which are matched
#'  with the real data using the random forest algorithm, in order to find the best fit of animal abundance
#'  among the abundance from 1 to 10 taken by each camera for one species. Such simulation can be repeated for
#'  several times defined by number of iteration.
#'
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param x A data.frame with column names "Lon", "Lat", "Group_size", "Date", "Time"
#' @param detect The detection radias (m) of a camera.
#' @param bearing The bearing direction of a camera.
#' @param step.N The number of steps the animal walks during the camera trapping.
#' @param step.L The mean step length (m) of the animal.
#' @param step.V The standard diviation of the step length
#' @param bias The standard diviation of the changing angle (degree) between two steps
#' @param range Maximum distance (m) the animal moves from the original site.
#' @param ind The number of individuals that are simulated.
#' @param iteration The number of simulations.
#'
#' @return A dataframe with the first column to be the number of individuals, and the rest columns
#'  are number of pictures (simulated) for each camera
#'
#' @examples
#'
#' par(mfrow = c(1, 2))
#' # maximum number of individuals in the camera grid is 10 (ind=10)
#' sim.out = simuCamtrap(trapresult, ind = 10, iteration = 2) # more iterations are expected for higher prediction accuracy
#'
#' @export
#'

simuCamtrap = function(x, detect = 50, bearing = runif(camera.N, 0, 2*pi), # detect = 0.200 / (40000/360) for latlon
                        step.N = 5000, step.V = 2,
                        step.L = 10, bias = 30/360*2*pi, range = 4000, ind = 10, # parameters were given based on real data
                        iteration = 3){

  # Movement parameters for species based on footprint chain data
  known <- data.frame(Species = c('Cape hare','Red deer','Roe deer','Moose','Chinese goral','Wild boar','Red fox',
                                  'Corsac fox','Raccoon dog','Eurasian lynx','Leopard cat','Yellow-throated marten','Sable','Siberian weasel'),
                       Step.L  = c(13.49,15.46,13.82,55.97,17.2,16.13,14.5,12.92,12.12,21.08,16.97,13.98,11.27,13.06),
                       Bias    = c(37.09,40.43,43.08,38.27,32.94,43.21,21.89,26.75,30.28,50.15,31.93,52.27,53.2,50.89), #SD of Angular deflection
                       Range   = c(0.25,1.95,0.66,4,9.2,1.15,2.2,0.9,2,8.81,4,2,2,1.2), # home range size (square km)
                       Ind     = c(20,30,40,20,20,50,10,15,20,5,10,20,20,20))  # maximum number of individuals in the camera grid

  if (!exists("known")) {
    stop("Please load movement parameters from dataset: known")
  }

  if (!x$Species[1] %in% as.character(known$Species)) {
    stop(paste(x$Species[1], "is not in the list. Movement parameters are not available", sep=""))
  }

  step.L = known[as.character(known$Species) == as.character(x$Species[1]), 'Step.L']
  bias   = known[as.character(known$Species) == as.character(x$Species[1]), 'Bias'] / 360*2*pi
  range  = known[as.character(known$Species) == as.character(x$Species[1]), 'Range'] * 1000
  # ind    = known[as.character(known$Species) == as.character(x$Species[1]), 'Ind'] # this parameter is open for adjusting

  camera = unique(x[,c('Lon','Lat')]) # x=trap.out
  camera = cbind(camera, Count=0)
  camera.N <- nrow(camera)
  out = as.data.frame(matrix(data = 0, nrow = ind*iteration, ncol = camera.N))


  for (ite in 1:iteration){#  iterations
    plotCamtrap(x)
    for (k in 1:ind){
      footchain = data.frame(ID=1:step.N, X=NA, Y=NA) # for plotting footchain
      loc.x = runif(1, min(x$Lon) + 0.25*(max(x$Lon) - min(x$Lon)), max(x$Lon)-0.25*(max(x$Lon) - min(x$Lon))) # random initial location
      loc.y = runif(1, min(x$Lat) + 0.25*(max(x$Lat) - min(x$Lat)), max(x$Lat)-0.25*(max(x$Lat) - min(x$Lat)))
      loc.x.0 = loc.x; loc.y.0 = loc.y

      ### MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM  simulating movement  MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
      theta  <- runif(1, 0, 2*pi) # moving bearing of the first step
      for (i in 1:step.N){   #walking steps
        for (j in 1:camera.N){ # number of cameras
          d = ((loc.x - camera$Lon[j])^2 + (loc.y - camera$Lat[j])^2)^.5
          bear = atan((loc.y - camera$Lat[j])/(loc.x - camera$Lon[j]))
          if(d < detect & abs(bear-bearing[j]) < 50/360*2*pi)   camera$Count[j] <- camera$Count[j] + 1 # 50 degree detetion region
        }
        move = step.L * rnorm(1, 1, step.V) #
        loc.x = loc.x + move*cos(theta) #
        loc.y = loc.y + move*sin(theta)
        theta.f = theta + rnorm(1,0,bias) # move forward
        if (loc.x > loc.x.0)    theta.b <- pi + atan((loc.y - loc.y.0)/(loc.x - loc.x.0))+ rnorm(1,0,bias*2)#return to origin
        if (loc.x < loc.x.0)    theta.b <-      atan((loc.y - loc.y.0)/(loc.x - loc.x.0))+ rnorm(1,0,bias*2)#return to origin
        dist = ((loc.x - loc.x.0)^2 + (loc.y - loc.y.0)^2)^.5
        theta = sample(c(theta.f, theta.b), 1, prob=c(abs((1-dist/range)), dist/range)) #abs() to make sure the p is positive
        footchain$X[i] = loc.x;   footchain$Y[i] = loc.y # for plotting footchain
      }
      # plot footprint chain
      lines(footchain$X, footchain$Y, col=rainbow(100)[sample(1:100, 1)])
      points(loc.x.0,loc.y.0, col='blue',pch=17,cex=.8)
      points(footchain$X[step.N], footchain$Y[step.N], col='red',pch=15,cex=.8)
      # points(camera$Lon, camera$Lat, col='black',pch=16)
      ### MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM  simulating movement  MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM

      camera$Count = sort(camera$Count) # NOT spatial explicit
      out[k+(ite-1)*ind,] = camera$Count            # generate out
      print(paste('Iteration ', ite, ';  Ind. ', k, sep=''))
    }
    camera$Count = camera$Count * 0
  }

  Ind = rep(1:ind, iteration)
  out = cbind(Ind = Ind, out)
  assign('sim.out', out, envir = .GlobalEnv)
  return(out)
}

library(compiler)
simuCamtrap <- cmpfun(simuCamtrap) # system.time(simuCamtrap(trap.out, ind = 3, iteration = 1))
#SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS

# x=trapresult





#PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
#' Predict number of individuals within the range of camera trap grid
#'
#' @description This function matches the real camera trap result with pseudo camera trap result simulated for a series number of individuals,
#'  to find the best fit of animal abundance with real camera trap result.
#'
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param simu The result of function simuCamtrap().
#' @param x A data.frame with column names "Lon", "Lat", "Group_size", "Date", "Time"
#' @param plot a boolean variable, if TRUE, plot the probability density of estimated animal abundance
#'
#' @return The mean value, 95% confidence intervals of the predicted animal abundance, based on
#'  a vector of the predicted animal abundance from 1000 random forest trees.
#'
#' @examples
#'
#'  attach(trapresult)
#'  # sim.out = simuCamtrap(trapresult, ind = 10, iteration = 2) # need a few minutes
#'  predictCamtrap(sim.out, trapresult, plot=T)
#'
#' @import randomForest
#' @export
#'

predictCamtrap = function(simu, x, plot=F){
  library(randomForest)
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop("randomForest needed for this function to work. Please install it.", call. = FALSE)
  }

  RF = randomForest(simu[,c(2:ncol(simu))], simu[,1], prox=TRUE, importance=F, ntree=1000)

  sum = aggregate(x$Group_size, by = list(x$Lat, x$Lon), sum)
  colnames(sum) = c('Lat','Lon','Count')
  obs = sort(sum$Count)
  pred <- predict(RF, obs, type="response", predict.all=TRUE)
  pred.rf.int <- apply( pred$individual, 1, function(x) { quantile(x, c(0.025, 0.5, 0.975) )})
  if(plot){
    plot(density(pred$individual),xlab='Number of individuals', ylab='Frequency',col='darkgrey',xlim=c(0,max(simu$Ind)),lwd=2,
         main=paste('Number of individuals:',round(pred.rf.int[2,1],1), sep=' '))
    abline(v = pred.rf.int[2,1], lwd=2) # median   # pred$aggregate # mean value of 1000 trees
    abline(v = c(pred.rf.int[1,1], pred.rf.int[3,1]), lwd=1, col='black',lty=2)
  }
  predicted = data.frame(No_ind=round(pred.rf.int[2,1],1),CI_L=round(pred.rf.int[1,1],1),CI_U=round(pred.rf.int[3,1],1))
  # write.csv(predicted, paste(x$Species[1], ".csv", sep=""), row.names = F)
  return(predicted)
}
predictCamtrap <- cmpfun(predictCamtrap)
#PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP








# Footprint chain analysis
## FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
#' Plot the trajectory of footprint chain
#'
#' @description Plot the trajectory of footprint chain, highlighting the starting and ending points.
#'
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param chain A data.frame with column names "Lon", "Lat", and "Date", in the order of recording time
#'  at the interval of one second
#'
#' @examples
#'
#'  attach(footprintchain)
#'  plotFootprint(footprintchain)
#'
#' @export
#'


plotFootprint = function(chain){
  plot(chain$Lon, chain$Lat, type='l', xlab='Longitude', ylab='Latitude')
  points(chain$Lon[1], chain$Lat[1], col = 'blue', cex=2, pch=17)
  points(chain$Lon[nrow(chain)], chain$Lat[nrow(chain)], col = 'red', cex=2, pch=15)
}
## FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF








## BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
#' Show how likely the animal would turn during moving
#'
#' @description Plot the probability density of step length, and probability density of turning angle between two steps.
#'  The step length can be defined by "scale" as the distance between every one second, two seconds, etc.
#'
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param chain A data.frame with column names "Lon", "Lat", and "Date", in the order of recording time
#'  at the interval of one second
#'
#' @param scale A value defining the step length
#'
#' @examples
#'
#'  attach(footprintchain)
#'  moveBias(footprintchain, scale=2, plot=T)
#'
#' @export
#'

moveBias = function(chain, scale, plot){
  chain = cbind(chain, dist=NA, theta=NA, delta.th=NA)
  chain = chain[seq(1, nrow(chain), by=scale),]
  N = nrow(chain)

  for (i in 1:(N-1)){
    chain$dist[i+1] = (((chain$Lat[i+1]-chain$Lat[i])*39946.79/360)^2+
                         ((chain$Lon[i+1]-chain$Lon[i])*pi*12756.32/360*cos(chain$Lat[i]*pi*2/360))^2)^0.5
    chain$theta[i+1] = asin((chain$Lat[i+1]-chain$Lat[i])*39946.79/360 / chain$dist[i+1]) *360/2/pi # bearing
  }

  chain = chain[!is.na(chain$theta),] #remove records with 0 distance
  chain$dist = chain$dist*1000
  chain = chain[!chain$dist>400,] # remove records across survey regions

  # Angular deflection
  N = nrow(chain);N
  for (j in 1:(N-1)){
    chain$delta.th[j+1] = chain$theta[j+1] - chain$theta[j]
  }
  chain = chain[-1,]
  if(plot){
     par(mfrow=c(1,2))
     plot(density(chain$dist),xlab="Step length (m)", main=paste('Mean step length: ', round(mean(chain$dist),2), "m"), cex=1.5)
     plot(density(chain$delta.th), xlab="Angular deflection (degree)", main=paste('SD of Angular deflection = ', round(sd(chain$delta.th),2)), cex=1.5)
  }

  out = c(round(mean(chain$dist),2), round(sd(chain$dist),2), round(sd(chain$delta.th),2))
  names(out) = c("Step length (m)", "SD of Step length (m)","SD of angular deflection")
  return(out)
}
## BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
#




## RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
#' Estimate number of individuals using Rowcliffe's equation (Rowcliffe 2008)
#'
#' @description Calculate animal density using Rowcliffe's random encounter models.
#'
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @param x A data.frame with column names "Lon", "Lat", "Group_size", "Date", "Time".
#' @param r Detection range (km) of the camera.
#' @param theta Detection angle of the camera
#' @param v Mean velocity (km/h) of the animal.
#' @param duration Days for camera trapping.
#'
#' @examples
#'
#' attach(trapresult)
#' Rowcliffe(trapresult, r=0.02, theta = 40, v = 2, duration = 40) # unit: km
#'
#' @export
#'

## population density (Rowcliffe 2008)
Rowcliffe = function(x, r, theta, v, duration){
  sum = aggregate(x$Group_size, by = list(x$Lat, x$Lon), sum)
  colnames(sum) = c('Lat','Lon','Count')
  D = sum$Count * pi / (duration*v*r*(2+theta*2*pi/360)) # population density (Rowcliffe 2008)
  density = c(round(mean(D),2), round(sd(D),2))
  names(density) = c('Density', 'SD')
  return(density)
}
## RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR



#' HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
#'
#' @description Calculate home ranges based on footprint chains
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#' @examples
#' attach(footprintchain)
#' HOME = homeRange(footprintchain); HOME
#'
#' @export
#'
homeRange = function(chain) {
  chain <- geo2proj(chain)
  chain$Lon = chain$Lon/1000; chain$Lat = chain$Lat/1000
  home = (max(chain$Lon)-min(chain$Lon)) * (max(chain$Lat)-min(chain$Lat))
  home = round(home, 3)
  names(home) = "Home range (square km)"
  return(home)
}
## HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH



#' DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
#' @description Calculate footprint chain length and number of steps
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#' @examples
#' attach(footprintchain)
#' footprintlength = chainLength(footprintchain); footprintlength
#'
#' @export
chainLength = function(chain) {
  chain = cbind(chain, Distance=NA)
  chain.U <- geo2proj(chain)
  chain.U$Lon = chain.U$Lon/1000; chain.U$Lat = chain.U$Lat/1000

  for (i in 2:nrow(chain)){
    chain$Distance[i] = ((chain.U$Lon[i]-chain.U$Lon[i-1])^2 + (chain.U$Lat[i]-chain.U$Lat[i-1])^2)^0.5
  }
  chain$Distance[chain$Distance>100] <- 0
  chain$Distance[1] = 0
  DIST = round(sum(chain$Distance), 3)
  names(DIST) = "Length of the footprint chain (km)"
  return(DIST)
}
## DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD









#' HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
#'
#' @description Calculate home ranges based on footprint chains
#'
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#'
#' @examples
#'
#' attach(chainDB)
#' HOME = homeRange(chainDB);HOME
#'
#' @export
#'
homeRanges = function(chain) {
  chain <- geo2proj(chain)
  chain$Lon = chain$Lon/1000; chain$Lat = chain$Lat/1000
  list = names(table(chain$Species))
  home = data.frame(Species=list, Area=NA, SD=NA, N=NA)
  for (i in 1:length(list)){
    spe = chain[chain$Species == list[i],]
    Lines = as.numeric(names(table(spe$Line)))
    Table = data.frame(Line = Lines, A=NA)
    for (j in 1:length(Lines)){
      spe_L = spe[spe$Line==Lines[j],]
      Table$A[j] = (max(spe_L$Lon)-min(spe_L$Lon)) * (max(spe_L$Lat)-min(spe_L$Lat))
    }
    home$Area[i] = round(mean(Table$A),3)
    home$SD[i] = round(sd(Table$A),3)
    home$N[i] = length(Table$A)
  }
  return(home)
}
## HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH



#' DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
#' @description Calculate footprint chain length and number of steps
#' @author Xinhai Li (Xinhai_li_edu@126.com)
#' @examples
#' attach(chainDB)
#' footprintlength = chainLength(chainDB); footprintlength
#'
#' @export
chainLengths = function(chain) {
  chain = cbind(chain, Distance=NA)
  chain.U <- geo2proj(chain)
  chain.U$Lon = chain.U$Lon/1000; chain.U$Lat = chain.U$Lat/1000

  for (i in 2:nrow(chain)){
    chain$Distance[i] = ((chain.U$Lon[i]-chain.U$Lon[i-1])^2 + (chain.U$Lat[i]-chain.U$Lat[i-1])^2)^0.5
  }
  chain$Distance[chain$Distance>100] <- 0
  chain$Distance[1] = 0

  list = names(table(chain$Species))
  DIST = data.frame(Species=list, Distance=NA, Dist_SD=NA, Step_L=NA, Step_No=NA, N=NA)
  for (i in 1:length(list)){
    spe = chain[chain$Species == list[i],]
    Lines = as.numeric(names(table(spe$Line)))
    Table = data.frame(Line = Lines, Dist=NA, Step=NA)
    for (j in 1:length(Lines)){
      spe_L = spe[spe$Line==Lines[j],]
      Table$Dist[j] = sum(spe_L$Distance)
      Table$Step[j] = moveBias(spe_L, scale=3, plot=F)[1] # geographic coordination system for lat/lon
    }
    DIST$Distance[i] = round(mean(Table$Dist),3)
    DIST$Dist_SD[i]  = round(sd(Table$Dist),3)
    DIST$Step_L[i]   = round(mean(Table$Step),1)
    DIST$Step_No[i]  = floor(DIST$Distance[i]*1000/DIST$Step_L[i]*30) # 30 days
    DIST$N[i] = length(Table$Dist)
  }
  return(DIST)
}
## DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
