dum_dtadata<-read.csv("HBFO_ts52x02_260616/dta260616_01.csv",skip = 8,header = T)
dtime_sam<-dum_dtadata$Relative.Time[which(dum_dtadata$Channel==3)]
dtemp_sam<-dum_dtadata$Reading[which(dum_dtadata$Channel==3)]
dtime_ref<-dum_dtadata$Relative.Time[which(dum_dtadata$Channel==4)]
dtemp_ref<-dum_dtadata$Reading[which(dum_dtadata$Channel==4)]

plot(dtime_sam,dtemp_sam,type = "l")
seq(0,max(dtime_sam),by=0.2)
dum_intpl<-approx(dtime_sam,dtemp_sam,xout = seq(0,max(dtime_sam),by=0.2))
lines(dum_intpl$x,dum_intpl$y,col="red")
lines(dum_intpl$x,dum_sg,col="blue")
install.packages("signal")
library(signal)
library(pracma)

dum_intplsig<-signal::interp1(dtime_sam,dtemp_sam,xi = seq(0,max(dtime_sam),by=0.2),method = "linear",extrap = T)
dum_sg<-savgol(dum_intplsig,31,2)

dum_intplref<-signal::interp1(dtime_ref,dtemp_ref,xi = seq(0,max(dtime_sam),by=0.2),method = "linear",extrap = T)
dum_sgref<-savgol(dum_intplref,31,2)
plot(dum_intplref[20:2270],dum_sg[20:2270] - dum_sgref[20:2270],type = "l")
lines(dum_intplref,dtemp_sam - dtemp_ref,type = "l",col="red")


dum_sgdtdt<-savgol(dum_intplsig,31,2,1)
plot(dum_intplref[20:2270],(dum_sg[20:2270] - dum_sgref[20:2270])/dum_sgdtdt[20:2270]*300,type = "l")
plot(dum_intplref[20:2270],dum_sgdtdt[20:2270]*300,type = "l")
lines(dum_intplref,dtemp_sam - dtemp_ref,type = "l",col="red")



dum_spldata<-splitData(dum_dtadata,T)
dum_rsdata<-rsmplData(dum_spldata,resist = T)
dum_sgdata<-sgsmoothdata(dum_rsdata,31,resist = T)

plot(dum_sgdata$Time[30:2260],dum_sgdata$Temp_sa[30:2260] - dum_rsdata$Temp_sa[30:2260],type = "l")

which(c("a","b","b")=="b")
seq_along(c("a","b","b"))
as.vector(t(outer(c("Tp_","Tb_","Tm_"), c("start","end"), paste0)))
