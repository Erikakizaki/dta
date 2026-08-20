library(pracma)
#return(list(method=method,method_id=method_id,dsc_df=res_df,temp_sec=temp_sec,coefs=params))
int_dsc<-function(res_list){
  int_list<-res_list
  
  dsc_df<-int_list$dsc_df
  
  t1<-int_list$temp_sec[1]
  t2<-int_list$temp_sec[2]
  
  pt1 <- which.min(abs(dsc_df$TempK - t1))
  pt2 <- which.min(abs(dsc_df$TempK - t2))
  
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
  dtdt<-gradient(dsc_df$TempK,dsc_df$Time)
  nom_HF<-dsc_df$HeatFlow_0/dtdt
  
  dsc_intseg<-dsc_df[intsegs,]
  H_int<-cumtrapz(dsc_intseg$Time*60,-dsc_intseg$HeatFlow_0)
  S_int<-cumtrapz(dsc_intseg$Time*60,-1000*dsc_intseg$HeatFlow_0/dsc_intseg$TempK)
  DH<-H_int[length(H_int)]
  DS<-S_int[length(S_int)]
  
  dsc_df<-cbind(dsc_df,dtdt,nom_HF)
  dsc_df$H<-0
  dsc_df$S<-0
  dsc_df$H[intsegs]<-H_int
  dsc_df$S[intsegs]<-S_int
  dsc_df$H[pt2:nrow(dsc_df)]<-DH
  dsc_df$S[pt2:nrow(dsc_df)]<-DS
  
  if(is.cooling){
    dsc_df$H<-dsc_df$H-DH
    dsc_df$S<-dsc_df$S-DS
    DH<- -DH
    DS<- -DS
  }
  
  int_list$dsc_df<-dsc_df
  int_list$DH<-DH
  int_list$DS<-DS
  
  return(int_list)
}