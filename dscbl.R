dscbl<-function(dscdataK,method,t1,t2,t3,t4){
  if(t1<t3 || t2>t4){
    warning("ピーク範囲がベースライン範囲に含まれることを推奨します")
  }
  TempK<-dscdataK$TempK
  HeatFlow<-dscdataK$HeatFlow
  
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
    basel<-HeatFlow[pt1]
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
    basel<-0
    method_id<-5
  }else{
    errorCondition("不正なmethodです")
  }
  
  
  
  
  HeatFlow_0<-HeatFlow-basel
  res_df<-cbind(dscbl,basel,HeatFlow_0)
  temp_sec<-c(t1,t2,t3,t4)
  
  return(list(method=method,method_id=method_id,dsc_df=res_df,temp_sec=temp_sec,coefs=params))
}

x<-1:100
which(x<=80&x>=60)
