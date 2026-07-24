#Redoing NYimore Data using "tidyverse" package

library(tidyverse)
setwd("C:\\Users\\korin\\OneDrive\\Desktop\\Summer Research - 2026 - Dr. Baron")
   NY = read.csv("NY_Week_Sample.csv")

#Keep Day, Latitude, Longitude

NY = NY[ , (names(NY) %in% c("Day", "Latitude", "Longitude"))]
N = na.omit(NY)     #Ntotal = 148558

N$Day = as.POSIXct(N$Day, format = "%Y-%m-%d")
N$t = as.numeric(N$Day)

N$X = as.numeric(N$Longitude)
N$Y = as.numeric(N$Latitude)

#Keep "X", "Y", "t"
N = N[ , (names(N) %in% c("X", "Y", "t"))]
N = na.omit(N)

# Delete Longitude outliers, those outside of h IQRs from Q1 and Q3

h = 1.5    # h IQRs
y1 = quantile(N$Y,.25) - h*IQR(N$Y)
y2 = quantile(N$Y,.75) + h*IQR(N$Y)
x1 = quantile(N$X,.25) - h*IQR(N$X)
x2 = quantile(N$X,.75) + h*IQR(N$X)

Z = as.numeric(N$Y > y1 & N$Y < y2 & N$X > x1 & N$X < x2)
N = N[Z==1,] #Ntotal = 142340

# NY center
xx = c(-74.26,-73.70); yy = c(40.49,40.92);

# Grid
Nbins = 30; DeltaX = (xx[2]-xx[1])/(Nbins); DeltaY = (yy[2]-yy[1])/(Nbins);

#######################################################################
### Estimate MRF parameters - 2D king size neighborhood ###############
#######################################################################

# Find standardized Longitudes Xbin and Ybin 
Xbin = ceiling((N$X - xx[1])/DeltaX)
Ybin = ceiling((N$Y - yy[1])/DeltaY)
Tbin = floor((N$t - min(N$t))/86400) + 1 #standardize days
Ndays   = max(Tbin)

NLongitudes = Nbins * Nbins
Events = rep(0, Ndays * NLongitudes)
dim(Events) = c(Ndays, Nbins, Nbins)
IndexEvents = which( Xbin >= 1 & Xbin <= Nbins & Ybin >= 1 & Ybin <= Nbins )

LongitudeX = rep(0, NLongitudes)
LongitudeY = LongitudeX
for (i in 1:Nbins){
  LongitudeX[(Nbins*(i-1)+1) :  (Nbins*i)]= i
  LongitudeY[(Nbins*(i-1)+1) : (Nbins*i)] = seq(1,Nbins)}

#Neighborhood
i0 = seq(2, (Nbins-1)) #The cell itself
i1 = seq(1,(Nbins-2)) #The cell's left or bottom neighbor
i2 = seq(3,Nbins) #The cells's right or top neighbor

#Define the 3D array
for (i in IndexEvents){Events[Tbin[i], Xbin[i], Ybin[i]] = 1}

Parameters = matrix(rep(0, Ndays*10), Ndays, 10)

for (t in 1:Ndays){
  DayEvents = Events[t,,]

#Using each Longitude as is
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

#Replacing each Longitude with the opposite
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

#And taking differences = Delta-statistics
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

EventResponse = as.factor(Nalpha)

if (length(unique(EventResponse)) < 2){
  Parameters[t,]= NA
  next()
}

logreg = glm( EventResponse ~ Dalpha + Dlambda + Dmu
              + Dpi + Drho + Deta + Dzeta + Dtheta + Diota + Dkappa, 
              family="binomial", 
              control=list(maxit = 500) ) 

Parameters[t,] = coef(logreg)[2:11]
}

Alpha  = Parameters[,1]
Lambda = Parameters[,2]
Mu     = Parameters[,3]
Pi     = Parameters[,4]
Rho    = Parameters[,5]
Eta    = Parameters[,6]
Zeta   = Parameters[,7]
Theta  = Parameters[,8]
Iota   = Parameters[,9]
Kappa  = Parameters[,10]
