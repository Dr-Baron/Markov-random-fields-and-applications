# DATA PRE-PROCESSING: data from https://catalog.data.gov/?q=911+calls&sort=relevance

setwd("C:\\Users\\baron\\Documents\\Research\\MRF\\Data and R")
 BALT = read.csv("Baltimore_911_Calls_for_Service.csv")
# BALT = read.csv("BaltimoreShort.csv")

# BALT = BALT[1:999999,]            # Only for coding speed

# Keep  "callDateTime" "priority" "location"       
# Erase "district" "description" "callNumber" "incidentLocation"
           
BALT = BALT[ , (names(BALT) %in% c("callDateTime", "priority", "location"))]

Ntotal = dim(BALT)[1]           # Ntotal = 2831903 calls

B = na.omit(BALT)

#B$month = as.numeric(substr(B$callDateTime,start=1,stop=2))
#B$day   = as.numeric(substr(B$callDateTime,start=4,stop=5))
#B$year  = as.numeric(substr(B$callDateTime,start=7,stop=10))
#B$hour  = as.numeric(substr(B$callDateTime,start=12,stop=13))
#B$min   = as.numeric(substr(B$callDateTime,start=15,stop=16))
#B$sec   = as.numeric(substr(B$callDateTime,start=18,stop=19))			# But sec = 0 always
#B$ampm  = substr(B$callDateTime,start=21,stop=22)
#B$hour = B$hour*(B$hour < 12) + 12*(B$hour < 12 & B$ampm=="PM") + 12*(B$hour == 12 & B$ampm=="PM")
#B$time = 60*B$min + 3600*B$hour

B$Day = as.Date(B$callDateTime,format="%m/%d/%Y")
B$t = as.numeric(B$Day)

# Check for any gaps in timeline
for (day in min(B$t):max(B$t)){ 
   if (sum(B$t==day)<10){ print(paste("Gap in timeline!!! Only",
         sum(B$t==day),"entries on day",day)) } 
                              }

# LATITUDE and LONGITUDE - from B$location of type " (lat,long) "
position1 = regexpr(",",B$location);    # Character position of ","
position2 = regexpr(")",B$location);    # Character position of ")"

B$Y = as.numeric(substr(B$location, 2, position1-1))              # Y = Latitude
B$X = as.numeric(substr(B$location, position1+1, position2-1))    # X = Longitude

# Erase "callDateTime" "location"     "hour"         "min"          "sec"          "ampm"         
# keep = c("year","month","day","time","Y","X","priority")

keep = c("t","X","Y","priority")
B = B[ , (names(B) %in% keep)]
B = na.omit(B)

# Delete location outliers, those outside of h IQRs from Q1 and Q3

h = 3    # h IQRs
y1 = quantile(B$Y,.25) - h*IQR(B$Y)
y2 = quantile(B$Y,.75) + h*IQR(B$Y)
x1 = quantile(B$X,.25) - h*IQR(B$X)
x2 = quantile(B$X,.75) + h*IQR(B$X)

Z = as.numeric(B$Y > y1 & B$Y < y2 & B$X > x1 & B$X < x2)

B = B[Z==1,]

B$level = 1*( B$priority=="Non-Emergency") + 3*( B$priority=="Low") + 4*( B$priority=="Medium") + 2*( B$priority=="High") + 6*( B$priority=="Emergency") 
# 1 = black, 2 = red, 3 = green, 4 = blue, 6 = maroon

plot(B$X, B$Y,xlim=c(x1,x2),ylim=c(y1,y2),col=B$level, pch=20)

# Baltimore center
xx = c(-76.67,-76.57); yy = c(39.27,39.32);

plot(B$X, B$Y,xlim=xx,ylim=yy)
points(B$X, B$Y,xlim=xx,ylim=yy,col=B$level, pch=20)
lines(xx,c(1,1)*yy[1] ,lwd=3); lines(xx,c(1,1)*yy[2] ,lwd=3); lines(xx[1]*c(1,1),yy ,lwd=3); lines(xx[2]*c(1,1),yy ,lwd=3);

# Grid
Nbins = 30; DeltaX = (xx[2]-xx[1])/(Nbins +1); DeltaY = (yy[2]-yy[1])/(Nbins +1);

for (k in 0 : (Nbins )){ lines(xx,c(1,1)*(yy[1]+k*DeltaY)); lines(c(1,1)*(xx[1]+k*DeltaX),yy); }
lines(xx,c(1,1)*yy[1] ,lwd=3); lines(xx,c(1,1)*yy[2] ,lwd=3); lines(xx[1]*c(1,1),yy ,lwd=3); lines(xx[2]*c(1,1),yy ,lwd=3);


#######################################################################
### Estimate MRF parameters - 2D king size neighborhood ###############
#######################################################################

# Find standardized locations Xbin and Ybin 
Xbin = ceiling((B$X - xx[1])/DeltaX);
Ybin = ceiling((B$Y - yy[1])/DeltaY);
Tbin = B$t - min(B$t) + 1;            # standardized time
Ndays   = max(Tbin);

## Create a 3-dimensional binary tensor in (Time, X, Y)

Events = array(0, c(Ndays, Nbins, Nbins) )

IndexEvents = which( Xbin >= 1 & Xbin <= Nbins & Ybin >= 1 & Ybin <= Nbins );


Nlocations = Nbins*Nbins;    # The 1st and the last rows and columns are
                             # used only to complete neigborhoods

# The next two functions transform a matrix into two vectors, X and Y
 LocationX = rep(0,Nlocations); LocationY = LocationX;
 for (i in 1:Nbins){ LocationX[(Nbins*(i-1)+1) : (Nbins*i)] = i;
                     LocationY[(Nbins*(i-1)+1) : (Nbins*i)] = seq(1,Nbins);
                  } 

# Neighboring locations
 i0 = seq(2,(Nbins-1));   # The point itself
 i1 = seq(1,(Nbins-2));   # Its left or bottom neighbor
 i2 = seq(3,Nbins);       # Its right or top neighbor


# Define the 3D tensor:
for (i in IndexEvents){ Events[ Tbin[i], Xbin[i], Ybin[i] ] = 1 }  # This is our Z


# Consider each day, save the logistic regression coefficients

Parameters = array(0, dim=c(Ndays,10))

# Check the number of events every day
Nevents = numeric(Ndays); for (t in 1:Ndays){ Nevents[t]=sum(Events[t,,]) }


for (t in 1:Ndays){
  DayEvents = Events[t,,];       # Nbins*Nbins matrix of 0s and 1s 
  # Define 10 Delta-statistics for each location, 
  # Each of them is a (Nbins-2)*(Nbins-2) matrix of 0s and 1s

# Using each location as is...
  Nalpha  = DayEvents[i0,i0]     
  Nlambda = DayEvents[i0,i0]*(DayEvents[i2,i0] + DayEvents[i1,i0])
  Nmu     = DayEvents[i0,i0]*(DayEvents[i0,i2] + DayEvents[i0,i1])
  Npi     = DayEvents[i0,i0]*(DayEvents[i2,i2] + DayEvents[i1,i1])
  Nrho    = DayEvents[i0,i0]*(DayEvents[i2,i1] + DayEvents[i1,i2])
  Neta    = DayEvents[i0,i0]*(DayEvents[i0,i2]*DayEvents[i1,i0] + DayEvents[i2,i0]*DayEvents[i2,i2] + DayEvents[i0,i1]*DayEvents[i1,i1])
  Nzeta   = DayEvents[i0,i0]*(DayEvents[i0,i1]*DayEvents[i2,i0] + DayEvents[i0,i2]*DayEvents[i2,i2] + DayEvents[i1,i0]*DayEvents[i1,i1])
  Ntheta  = DayEvents[i0,i0]*(DayEvents[i0,i1]*DayEvents[i1,i0] + DayEvents[i0,i2]*DayEvents[i2,i1] + DayEvents[i0,i2]*DayEvents[i1,i2])
  Niota   = DayEvents[i0,i0]*(DayEvents[i0,i2]*DayEvents[i2,i0] + DayEvents[i0,i1]*DayEvents[i2,i1] + DayEvents[i1,i0]*DayEvents[i1,i2])
  Nkappa  = DayEvents[i0,i0]*(DayEvents[i0,i2]*DayEvents[i2,i0]*DayEvents[i2,i2] + DayEvents[i0,i1]*DayEvents[i2,i0]*DayEvents[i2,i1] + 
                              DayEvents[i0,i1]*DayEvents[i1,i0]*DayEvents[i1,i1] + DayEvents[i0,i2]*DayEvents[i1,i0]*DayEvents[i1,i2])

# Replacing each location with the opposite...
  nalpha  = (1-DayEvents[i0,i0])     
  nlambda = (1-DayEvents[i0,i0])*(DayEvents[i2,i0] + DayEvents[i1,i0])
  nmu     = (1-DayEvents[i0,i0])*(DayEvents[i0,i2] + DayEvents[i0,i1])
  npi     = (1-DayEvents[i0,i0])*(DayEvents[i2,i2] + DayEvents[i1,i1])
  nrho    = (1-DayEvents[i0,i0])*(DayEvents[i2,i1] + DayEvents[i1,i2])
  neta    = (1-DayEvents[i0,i0])*(DayEvents[i0,i2]*DayEvents[i1,i0] + DayEvents[i2,i0]*DayEvents[i2,i2] + DayEvents[i0,i1]*DayEvents[i1,i1])
  nzeta   = (1-DayEvents[i0,i0])*(DayEvents[i0,i1]*DayEvents[i2,i0] + DayEvents[i0,i2]*DayEvents[i2,i2] + DayEvents[i1,i0]*DayEvents[i1,i1])
  ntheta  = (1-DayEvents[i0,i0])*(DayEvents[i0,i1]*DayEvents[i1,i0] + DayEvents[i0,i2]*DayEvents[i2,i1] + DayEvents[i0,i2]*DayEvents[i1,i2])
  niota   = (1-DayEvents[i0,i0])*(DayEvents[i0,i2]*DayEvents[i2,i0] + DayEvents[i0,i1]*DayEvents[i2,i1] + DayEvents[i1,i0]*DayEvents[i1,i2])
  nkappa  = (1-DayEvents[i0,i0])*(DayEvents[i0,i2]*DayEvents[i2,i0]*DayEvents[i2,i2] + DayEvents[i0,i1]*DayEvents[i2,i0]*DayEvents[i2,i1] + 
                              DayEvents[i0,i1]*DayEvents[i1,i0]*DayEvents[i1,i1] + DayEvents[i0,i2]*DayEvents[i1,i0]*DayEvents[i1,i2])
    
# And taking differences = Delta-statistics...
  Dalpha  = -as.numeric(as.vector(Nalpha - nalpha))
  Dlambda = -as.numeric(as.vector(Nlambda - nlambda))
  Dmu     = -as.numeric(as.vector(Nmu - nmu))
  Dpi     = -as.numeric(as.vector(Npi - npi))
  Drho    = -as.numeric(as.vector(Nrho - nrho))
  Deta    = -as.numeric(as.vector(Neta - neta))
  Dzeta   = -as.numeric(as.vector(Nzeta - nzeta))
  Dtheta  = -as.numeric(as.vector(Ntheta - ntheta))
  Diota   = -as.numeric(as.vector(Niota - niota))
  Dkappa  = -as.numeric(as.vector(Nkappa - nkappa))

  EventResponse = as.factor(Nalpha);

  logreg = glm( EventResponse ~ Dalpha + Dlambda + Dmu
    + Dpi + Drho + Deta + Dzeta + Dtheta + Diota + Dkappa, 
    family="binomial", 
    control=list(maxit = 500) ) 

  Parameters[t,] = coef(logreg)[2:11]
}

Parameters = data.frame(Parameters)
names(Parameters)=c("Alpha","Lambda","Mu","Pi","Rho","Eta","Zeta","Theta","Iota","Kappa")
write.csv(Parameters,"MRF param estim from log regr.csv")

# par(mfrow=c(5,2))
# plot(Alpha,main="Alpha")
# plot(Lambda,main="Lambda")
# plot(Mu,main="Mu")
# plot(Pi,main="Pi")
# plot(Rho,main="Rho")
# plot(Eta,main="Eta")
# plot(Zeta,main="Zeta")
# plot(Theta,main="Theta")
# plot(Iota,main="Iota")
# plot(Kappa,main="Kappa")

write.csv(Parameters,"MRF param estim from log regr.csv")

########################################################################
### Listing all (T,X,Y) combinations and relating them to the line number
########################################################################
#
# Ndays = 5; Nbins = 3; Nbins = 3;
# N = Ndays*Nbins*Nbins;
# T = rep(0,N); X=T; Y=T;
# for (t in 1:Ndays){ T[((t-1)*Nbins*Nbins + 1) : (t*Nbins*Nbins)] = t; 
#   for (x in 1:Nbins){ X[((t-1)*Nbins*Nbins + (x-1)*Nbins + 1) : ((t-1)*Nbins*Nbins + x*Nbins)] = x;
#   Y[((t-1)*Nbins*Nbins + (x-1)*Nbins + 1) : ((t-1)*Nbins*Nbins + x*Nbins)] = 1:Nbins }}
# Z = rep(0,N)
# data.frame(T,X,Y,Z)
# t=2; x=2; y=2;
# (t-1)*Nbins^2 + (x-1)*Nbins + y
# t=5; x=3; y=1;
# (t-1)*Nbins^2 + (x-1)*Nbins + y

