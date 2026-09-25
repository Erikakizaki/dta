mirror_x <- scale_x_continuous(sec.axis = dup_axis(labels = NULL, name = NULL))
mirror_y <- scale_y_continuous(sec.axis = dup_axis(labels = NULL, name = NULL))

mytheme <- theme_classic(base_family = "serif", base_size = 16) +
  theme(
    axis.ticks.length = unit(-0.13, "inches"),
    axis.text.x.top = element_blank(),
    axis.text.y.right = element_blank()
  )
dtabl<-function(dtadataK,method,t1,t2,t3,t4){
  if(t1<t3 || t2>t4){
    warning("ピーク範囲がベースライン範囲に含まれることを推奨します")
  }
  TempK<-dtadataK$Temp_re
  HeatFlow<-dtadataK$hf
  
  basel<-HeatFlow
  params<-rep(0,5)
  # T1, T2に最も近いデータ点を抽出
  pt1 <- which.min(abs(TempK - t1))
  pt2 <- which.min(abs(TempK - t2))
  pt3 <- which.min(abs(TempK - t3))
  pt4 <- which.min(abs(TempK - t4))
  #poly用section
  sec1<-which(TempK>=TempK[pt3]&TempK<=TempK[pt4]&(TempK<=TempK[pt1]|TempK>=TempK[pt2]))
  
  if(method=="2 points"){
    # 2点間の直線式パラメータ
    params[2] <- (HeatFlow[pt2] - HeatFlow[pt1]) / (TempK[pt2] - TempK[pt1])
    params[1] <- HeatFlow[pt1] - params[2] * TempK[pt1]
    
    #T<t1: HF at t1
    basel[seq_len(basel)]<-HeatFlow[pt1]
    #t1<=T<=t2: line between pt1 and pt2
    sec1<-which(TempK>=TempK[pt1])
    basel[sec1]<-params[1]+params[2]*TempK[sec1]
    #T>t2: HF at pt2
    sec2<-which(TempK>TempK[pt2])
    basel[sec2]<-HeatFlow[pt2]
    
    method_id<-0
  }else if(method=="line"){
    fit_sml<-lm(HeatFlow[sec1]~TempK[sec1])
    
    params[c(1,2)]<-fit_sml$coefficients
    basel<-params[1]+params[2]*TempK
    method_id<-1
  }else if(method=="poly(2nd)"){
    fit_sml<-lm(HeatFlow[sec1]~poly(TempK[sec1],2,raw = T))
    
    params[c(1,2,3)]<-fit_sml$coefficients
    basel<-params[1]+params[2]*TempK+params[3]*TempK^2
    method_id<-2
  }else if(method=="poly(3rd)"){
    fit_sml<-lm(HeatFlow[sec1]~poly(TempK[sec1],3,raw = T))
    
    params[c(1,2,3,4)]<-fit_sml$coefficients
    basel<-params[1]+params[2]*TempK+params[3]*TempK^2+params[4]*TempK^3
    method_id<-3
  }else if(method=="poly(4th)"){
    fit_sml<-lm(HeatFlow[sec1]~poly(TempK[sec1],4,raw = T))
    
    params[c(1,2,3,4,5)]<-fit_sml$coefficients
    basel<-params[1]+params[2]*TempK+params[3]*TempK^2+params[4]*TempK^3+params[5]*TempK^4
    method_id<-4
  }else if(method=="2 line 2nd connection"){
    cat("準備中")
    basel<-basel-basel
    method_id<-5
  }else{
    errorCondition("不正なmethodです")
  }
  
  
  
  
  HeatFlow_0<-HeatFlow-basel
  res_df<-cbind(dtadataK,basel,HeatFlow_0)
  temp_sec<-c(t1,t2,t3,t4)
  
  return(list(method=method,method_id=method_id,dta_df=res_df,temp_sec=temp_sec,coefs=params))
}

#return(list(method=method,method_id=method_id,dta_df=res_df,temp_sec=temp_sec,coefs=params))
int_dta<-function(res_list,Scoef=1){
  int_list<-res_list
  
  dta_df<-int_list$dta_df
  t1<-int_list$temp_sec[1]
  t2<-int_list$temp_sec[2]
  pt1 <- which.min(abs(dta_df$Temp_re - t1))
  pt2 <- which.min(abs(dta_df$Temp_re - t2))
  
  if(pt1<pt2){
    #heating
    intsegs<-(pt1:pt2)
    is.cooling<-F
  }else{
    #coolingは積分値を低温側0にする処理を加える
    intsegs<-(pt2:pt1)
    is.cooling<-T
  }
  
  nom_HF<-dta_df$HeatFlow_0*Scoef
  dta_df<-cbind(dta_df,nom_HF)
  
  dta_intseg<-dta_df[intsegs,]
  H_int<-cumtrapz(dta_intseg$Temp_re,-dta_intseg$nom_HF)
  S_int<-cumtrapz(dta_intseg$Temp_re,-1000*dta_intseg$nom_HF/dta_intseg$Temp_re)
  DH<-H_int[length(H_int)]
  DS<-S_int[length(S_int)]
  
  
  dta_df$H<-0
  dta_df$S<-0
  dta_df$H[intsegs]<-H_int
  dta_df$S[intsegs]<-S_int
  dta_df$H[max(pt1,pt2):nrow(dta_df)]<-DH
  dta_df$S[max(pt1,pt2):nrow(dta_df)]<-DS
  
  if(is.cooling){
    dta_df$H<-dta_df$H-DH
    dta_df$S<-dta_df$S-DS
    DH<- -DH
    DS<- -DS
  }
  
  int_list$dta_df<-dta_df
  int_list$DH<-DH
  int_list$DS<-DS
  return(int_list)
}

calcoef<-function(res_list, stdentropy){
  int_list<-res_list
  
  dta_df<-int_list$dta_df
  t1<-int_list$temp_sec[1]
  t2<-int_list$temp_sec[2]
  pt1 <- which.min(abs(dta_df$Temp_re - t1))
  pt2 <- which.min(abs(dta_df$Temp_re - t2))
  
  if(pt1<pt2){
    #heating
    intsegs<-(pt1:pt2)
    is.cooling<-F
  }else{
    #coolingは積分値を低温側0にする処理を加える
    intsegs<-(pt2:pt1)
    is.cooling<-T
  }
  
  
  #dT/dt
  #dtdt<-gradient(dta_df$TempK,dta_df$Time)
  #nom_HF<-dta_df$HeatFlow_0/dtdt
  
  dta_intseg<-dta_df[intsegs,]
  #integral under unit [s], not [min]
  #H_int<-cumtrapz(dta_intseg$Time*60,-dta_intseg$HeatFlow_0)
  S_int<-cumtrapz(dta_intseg$Temp_re,-1000*dta_intseg$HeatFlow_0/dta_intseg$Temp_re)
  #DH<-H_int[length(H_int)]
  DS<-S_int[length(S_int)]
  
  # dta_df<-cbind(dta_df,dtdt,nom_HF)
  # dta_df$H<-0
  # dta_df$S<-0
  # dta_df$H[intsegs]<-H_int
  # dta_df$S[intsegs]<-S_int
  # dta_df$H[max(pt1,pt2):nrow(dta_df)]<-DH
  # dta_df$S[max(pt1,pt2):nrow(dta_df)]<-DS
  
  if(is.cooling){
    # dta_df$H<-dta_df$H-DH
    # dta_df$S<-dta_df$S-DS
    # DH<- -DH
    DS<- -DS
  }
  
  # int_list$dta_df<-dta_df
  # int_list$DH<-DH
  # int_list$DS<-DS
  return(stdentropy/DS)
}



splitData<-function(dtadata, resist=F, C2K=T){
  #length adjustment
  rowex<-which(dtadata$Channel==1)
  rowsa<-which(dtadata$Channel==3)
  rowre<-which(dtadata$Channel==4)
  minlen<-min(c(length(rowex),length(rowsa),length(rowre)))
  maxlen<-max(c(length(rowex),length(rowsa),length(rowre)))
  rowex<-rowex[1:minlen]
  rowsa<-rowsa[1:minlen]
  rowre<-rowre[1:minlen]
  
  if(resist){
    rowma<-which(dtadata$Channel==5)
    minlen<-min(c(length(rowex),length(rowsa),length(rowre),length(rowma)))
    maxlen<-max(c(length(rowex),length(rowsa),length(rowre),length(rowma)))
    rowex<-rowex[1:minlen]
    rowsa<-rowsa[1:minlen]
    rowre<-rowre[1:minlen]
    rowma<-rowma[1:minlen]
  }
  
  #data splitting
  
  data_ex<-dtadata[rowex,]
  data_sa<-dtadata[rowsa,]
  data_re<-dtadata[rowre,]
  tmp_ex<-data_ex$Reading
  tmp_sa<-data_sa$Reading
  tmp_re<-data_re$Reading
  tim_ex<-data_ex$Relative.Time
  tim_sa<-data_sa$Relative.Time
  tim_re<-data_re$Relative.Time
  
  if(C2K){
    tmp_ex<-tmp_ex+273.15
    tmp_sa<-tmp_sa+273.15
    tmp_re<-tmp_re+273.15
  }
  
  spldata<-cbind(tim_ex,tmp_ex,tim_sa,tmp_sa,tim_re,tmp_re)
  colnames(spldata)<-c("Time_ex","Temp_ex","Time_sa","Temp_sa","Time_re","Temp_re")
  if(resist){
    data_ma<-dtadata[rowma,]
    res_ma<-data_ma$Reading
    tim_ma<-data_ma$Relative.Time
    
    splres<-cbind(tim_ma,res_ma)
    colnames(splres)<-c("Time_ma","Resist")
    
    spldata<-cbind(spldata,splres)
  }
  
  return(as.data.frame(spldata))
}

rsmplData<-function(spldata, dt=0.2, resist=F){
  #resampling
  tim_rs<-seq(0,max(spldata$Time_ex),by=dt)
  
  tmp_exrs<-signal::interp1(spldata$Time_ex,spldata$Temp_ex,xi = tim_rs,method = "linear",extrap = T)
  tmp_sars<-signal::interp1(spldata$Time_sa,spldata$Temp_sa,xi = tim_rs,method = "linear",extrap = T)
  tmp_rers<-signal::interp1(spldata$Time_re,spldata$Temp_re,xi = tim_rs,method = "linear",extrap = T)
  
  rsdata<-cbind(tim_rs,tmp_exrs,tmp_sars,tmp_rers)
  colnames(rsdata)<-c("Time","Temp_ex","Temp_sa","Temp_re")
  
  if(resist){
    Resist<-signal::interp1(spldata$Time_ma,spldata$Resist,xi = tim_rs,method = "linear",extrap = T)
    rsdata<-cbind(rsdata,Resist)
  }
  
  return(as.data.frame(rsdata) )
  
  
}

sgsmoothdata<-function(rsdata,sgwindow,dt=0.2,resist=F){
  #rsdata consists of 4 (if resist=T, 5) columns; Time, Temp_ex, Temp_sa, Temp_re (, Resist).
  #SG filter
  
  sgdata<-rsdata
  sgdata$Temp_ex<-savgol(rsdata$Temp_ex,sgwindow,2)
  sgdata$Temp_sa<-savgol(rsdata$Temp_sa,sgwindow,2)
  sgdata$Temp_re<-savgol(rsdata$Temp_re,sgwindow,2)
  
  if(resist){
    sgdata$Resist<-savgol(rsdata$Resist,sgwindow,2)
  }
  
  
  #dtdt
  #*5 means sg filter span (1) / actual data span
  dtdt<-savgol(rsdata$Temp_re,sgwindow,2,1)/dt
  
  #heat flow
  hf<-(rsdata$Temp_sa-rsdata$Temp_re)/dtdt
  
  sgdata<-cbind(sgdata,dtdt,hf)
  
  dum_nrow<-nrow(sgdata)
  
  sgdata<-sgdata[sgwindow:(dum_nrow - sgwindow),]
  
  
  return(as.data.frame(sgdata))
}


