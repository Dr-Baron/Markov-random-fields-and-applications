# Get full Seattle Data

setwd("C:\\Users\\korin\\OneDrive\\Desktop\\Summer Research - 2026 - Dr. Baron")
  SEATTLE = read.csv("Seattle_911_calls_for_service.csv")

# Keep Datetime, Latitude, Longitude
  
S = SEATTLE[ , (names(SEATTLE) %in% c("Datetime", "Latitude", "Longitude"))]
S = na.omit(S)

S$Day = as.POSIXct(S$Datetime, format = "%Y %b %d %I:%M:%S %p")
S$t = floor(as.numeric(S$Day) / 86400)

# Location
h = 3    # h IQRs
y1 = quantile(S$Latitude,.25) - h*IQR(S$Latitude)
y2 = quantile(S$Latitude,.75) + h*IQR(S$Latitude)
x1 = quantile(S$Longitude,.25) - h*IQR(S$Longitude)
x2 = quantile(S$Longitude,.75) + h*IQR(S$Longitude)

Z = as.numeric(S$Latitude > y1 & S$Latitude < y2 & S$Longitude > x1 & S$Longitude < x2)
S = S[Z==1,]


# Seattle center
xx = c(-122.355, -122.314); yy = c(47.589, 47.629)

plot(S$Longitude, S$Latitude, xlim = xx, ylim = yy)
points(S$Longitude, S$Latitude, xlim = xx, ylim = yy)
lines(xx,c(1,1)*yy[1] ,lwd=3); lines(xx,c(1,1)*yy[2] ,lwd=3); lines(xx[1]*c(1,1),yy ,lwd=3); lines(xx[2]*c(1,1),yy ,lwd=3);

# Grid
Nbins = 30; DeltaX = (xx[2]-xx[1])/(Nbins); DeltaY = (yy[2]-yy[1])/(Nbins);

for (k in 0 : (Nbins )){ lines(xx,c(1,1)*(yy[1]+k*DeltaY)); lines(c(1,1)*(xx[1]+k*DeltaX),yy); }
lines(xx,c(1,1)*yy[1] ,lwd=3); lines(xx,c(1,1)*yy[2] ,lwd=3); lines(xx[1]*c(1,1),yy ,lwd=3); lines(xx[2]*c(1,1),yy ,lwd=3);