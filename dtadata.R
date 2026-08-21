wd<-commandArgs(trailingOnly=T)
setwd(wd[1])

sigmoid2d <- function(x) {
  y <- c()
  for (i in x) {
    if (i == 1 / 2) {
      x1 <- 1 / 2
    } else if (i < 0) {
      x1 <- 0
    } else if (i > 1) {
      x1 <- 1
    } else{
      x2 <- i - 1 / 2
      x1 <- (1 / 2 - 2 * (abs(x2) - 1 / 2)^2) * x2 / abs(x2) + 1 / 2
    }
    y <- c(y, x1)
  }
  return(y)
}

abconnect<-function(a1,b1,a2,b2,x0,x1,ximput,nzrange=NA){
  line1<-a1+b1*ximput
  line2<-a2+b2*ximput
  normx<-(ximput-x0)/(x1-x0)
  sigweight<-sigmoid2d(normx)
  y<-line1*(1-sigweight)+line2*sigweight
  if(length(nzrange>=2)){
    y<-as.numeric(ximput>=nzrange[1])*as.numeric(ximput<=nzrange[2])*y
  }
  return(y)
}

procData<-function(hoge,resist=F){
  #length adjustment
  rowex<-which(hoge$Channel==1)
  rowsa<-which(hoge$Channel==3)
  rowre<-which(hoge$Channel==4)
  minlen<-min(c(length(rowex),length(rowsa),length(rowre)))
  maxlen<-max(c(length(rowex),length(rowsa),length(rowre)))
  rowex<-rowex[1:minlen]
  rowsa<-rowsa[1:minlen]
  rowre<-rowre[1:minlen]
  
  if(resist){
    rowma<-which(hoge$Channel==5)
    minlen<-min(c(length(rowex),length(rowsa),length(rowre),length(rowma)))
    maxlen<-max(c(length(rowex),length(rowsa),length(rowre),length(rowma)))
    rowex<-rowex[1:minlen]
    rowsa<-rowsa[1:minlen]
    rowre<-rowre[1:minlen]
    rowma<-rowma[1:minlen]
  }
  
  #data splitting
  
  data_ex<-hoge[rowex,]
  data_sa<-hoge[rowsa,]
  data_re<-hoge[rowre,]
  tmp_ex<-data_ex$Reading
  tmp_sa<-data_sa$Reading
  tmp_re<-data_re$Reading
  tim_ex<-data_ex$Relative.Time
  tim_sa<-data_sa$Relative.Time
  tim_re<-data_re$Relative.Time
  
  if(resist){
    data_ma<-hoge[rowma,]
    res_ma<-data_ma$Reading
    tim_ma<-data_ma$Relative.Time
  }
  
  #bspline interpolation
  
  sim_sa<-LSM_bspline(tim_sa,tmp_sa,xi=0.001,dt_in = 6)
  sim_re<-LSM_bspline(tim_re,tmp_re,xi=0.001,dt_in = 6)
  
  if(resist){
    sim_ma<-LSM_bspline(tim_ma,res_ma,xi=0.001,dt_in = 6)
  }
  
  #differential
  
  dtim_re<-c()
  dtmpdtim_re<-c()
  
  for(i0 in 1:(minlen-1)){
    dtim_re<-c(dtim_re,(tim_re[i0]+tim_re[i0+1])/2)
    dtmpdtim_re<-c(dtmpdtim_re,60*(tmp_re[i0+1]-tmp_re[i0])/(tim_re[i0+1]-tim_re[i0]))
  }
  
  sim_dps<-seq(2,nrow(sim_re)-1,by=2)
  sim_dtim<-sim_re[sim_dps,1]
  sim_dtmpdtim<-c()
  
  #dtout=0.05
  #"*6" means dt=0.05*2,K/min->*60
  for(i1 in sim_dps){
    sim_dtmpdtim<-c(sim_dtmpdtim,(sim_re[i1+1,2]-sim_re[i1-1,2])*600)
  }

  #summarise
  rawdata<-data.frame(tmp_ex=tmp_ex,tmp_sa=tmp_sa,tmp_re=tmp_re,
                      tim_ex=tim_ex,tim_sa=tim_sa,tim_re=tim_re,
                      dtmp=tmp_sa-tmp_re)
  
  intpl<-data.frame(tmp_sa=sim_sa[,2],tmp_re=sim_re[,2],
                    tim_sa=sim_sa[,1],tim_re=sim_re[,1],
                    dtmp=sim_sa[,2]-sim_re[,2])
  dif_exp<-data.frame(dtim=dtim_re,dtdt=dtmpdtim_re)
  dif_sim<-data.frame(dtim=sim_dtim,dtdt=sim_dtmpdtim)
  if(resist){
    rawdata$res<-res_ma
    rawdata$tim_ma<-tim_ma
    intpl$res<-sim_ma[,2]
    intpl$tim_ma<-sim_ma[,1]
  }
  
  return(list(rawdata,intpl,dif_exp,dif_sim))
                 
}

bspline<-function(z,h){
  bvals<-c()
  for(i in z){
    zx<-2-abs(i/h-2)
    
    if(zx>0&zx<=1){
      bval<-zx^3/6
    }else if(zx>1&zx<=2){
      bval<-(-zx)^3/2+2*zx^2-2*zx+2/3
    }else{
      bval<-0
    }
    bvals<-c(bvals,bval)
  }
  return(bvals)
}

idMat<-function(N){
  x<-matrix(0,nrow = N,ncol = N)
  diag(x)<-1
  return(x)
}

LSM_bspline<-function(Xtime,Ytemp,dt_in=1,dt_out=0.05,xi){
  #Ymean<-mean(Ytemp)
  #Ytemp<-Ytemp-Ymean
  tmax<-max(Xtime)
  tmin<-min(Xtime)
  nnmin<-floor(tmin/dt_in)
  nnmax<-floor(tmax/dt_in)+1
  dim_node<-floor(tmax/dt_in)-nnmin+4
  dim_obs<-length(Xtime)
  
  Hmat<-matrix(0,nrow = dim_obs,ncol = dim_node)
  #return(c(dim_obs,dim_node,nnmin))
  for(i0 in 1:dim_obs){
    nnbase<-floor(Xtime[i0]/dt_in)-nnmin
    trel<-Xtime[i0]-floor(Xtime[i0]/dt_in)*dt_in
    #print(nnbase)
    for(i1 in 1:4){
      Hmat[i0,nnbase+i1]<-bspline(trel+(4-i1)*dt_in,dt_in)
    }
  }
  ansbeta<-solve(t(Hmat)%*%Hmat+xi^2*idMat(dim_node))%*%(t(Hmat)%*%Ytemp)
  ansbeta<-as.vector(ansbeta)
  
  ranl<-nnmin*dt_in
  ranh<-nnmax*dt_in-1
  
  tsim<-seq(ranl,ranh,by=dt_out)
  tempsim<-c()
  for(i2 in tsim){
    tempsim0<-0
    nnbsim<-floor(i2/dt_in)-nnmin
    trelsim<-i2-floor(i2/dt_in)*dt_in
    for(i3 in 1:4){
      tempsim0<-tempsim0+bspline(trelsim+(4-i3)*dt_in,dt_in)*ansbeta[nnbsim+i3]
    }
    tempsim<-c(tempsim,tempsim0)
  }
  # plot(Xtime,Ytemp)
  # lines(tsim,tempsim)
  return(cbind(tsim,tempsim))
  
  # tempsim<-c()
  # for(i4 in Xtime){
  #   tempsim0<-0
  #   nnbsim<-floor(i4/dt_in)-nnmin
  #   trelsim<-i4-floor(i4/dt_in)*dt_in
  #   for(i5 in 1:4){
  #     tempsim0<-tempsim0+bspline(trelsim+(4-i5)*dt_in,dt_in)*ansbeta[nnbsim+i5]
  #   }
  #   tempsim<-c(tempsim,tempsim0)
  # }
  # plot(Xtime,Ytemp-tempsim,ylim=c(-0.01,0.01))
  
  
}

#ここから実行

setting<-read.table("setting.txt",comment.char = "#")
datlist<-read.csv("datlist.csv")
csvname<-datlist$CSVNAME
exppress<-datlist$PRESSURE
horc<-datlist$COMMENT
mintemp<-datlist$MINTEMP
maxtemp<-datlist$MAXTEMP

lptimes<-length(csvname)

alldata<-list()
svdata<-list()
for (i6 in 1:lptimes) {
  dtadata<-read.csv(csvname[i6],skip = 8)
  alldata<-procData(dtadata,resist = as.logical(setting[3,1]))
  if(horc[i6]=="h"){
    dirhead<-"heating"
  }else{
    dirhead<-"cooling"
  }
  mintemp[i6]<-min(alldata[[1]]$tmp_re)
  maxtemp[i6]<-max(alldata[[1]]$tmp_re)
  write.csv(alldata[[1]],paste0(dirhead,"/r_",csvname[i6]),row.names = F)
  write.csv(alldata[[2]],paste0(dirhead,"/s_",csvname[i6]),row.names = F)
  write.csv(alldata[[3]],paste0(dirhead,"/dr_",csvname[i6]),row.names = F)
  write.csv(alldata[[4]],paste0(dirhead,"/ds_",csvname[i6]),row.names = F)
  
  svdata[[i6]]<-alldata
  
}
datlist$MINTEMP<-mintemp
datlist$MAXTEMP<-maxtemp

# x MPa to a+bx kbar
MPa2kba<-0
MPa2kbb<-0.2


write.csv(datlist,"datlist.csv",row.names = F)

save(svdata,setting,datlist,lptimes,file = "continue.RData")
rm(list = ls())

#dtagraph
load("continue.RData")
#svdata:
#setting:species,no,resistance
#datlist:csvname,exppress,horc,mintemp,maxtemp




#(ax+b)(1-s((x-c)/d))+

dtacurve<-function(x,pms){
  bl1<-(pms[1]*x+pms[2])*(1-sigmoid2d((x-pms[5])/pms[6]))
  bl2<-(pms[3]*x+pms[4])*(sigmoid2d((x-pms[5])/pms[6]))
  pk<-pms[7]*sigmoid2d((x-pms[8])/pms[9])*(1-sigmoid2d((x-pms[10])/pms[11]))
  return(bl1+bl2+pk)
}



#which data is heating or cooling
hdlist<-c()
cdlist<-c()
horc<-datlist$COMMENT
#inpress<-datlist$INPRESS

for(i7 in 1:lptimes){
  if(horc[i7]=="h"){
    hdlist<-c(hdlist,i7)
    fcolor<-"#ff00ff"
  }else{
    cdlist<-c(cdlist,i7)
    fcolor<-"#00ffff"
  }
  windows(8,5)
  
  
  
}
par(mygp)
mygp$ps<-15
mygp$lwd<-2

