library(shiny)
library(ggplot2)
library(tidyr)
library(patchwork)



# ggplot2用の共通テーマ・軸設定 (functions_2.R より)
mirror_x <- scale_x_continuous(sec.axis = dup_axis(labels = NULL, name = NULL))
mirror_y <- scale_y_continuous(sec.axis = dup_axis(labels = NULL, name = NULL))

mytheme <- theme_classic(base_family = "serif", base_size = 16) +
  theme(
    axis.ticks.length = unit(-0.13, "inches"),
    axis.text.x.top = element_blank(),
    axis.text.y.right = element_blank()
  )
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
  res_df<-cbind(dscdataK,basel,HeatFlow_0)
  temp_sec<-c(t1,t2,t3,t4)
  
  return(list(method=method,method_id=method_id,dsc_df=res_df,temp_sec=temp_sec,coefs=params))
}
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
  #integral under unit [s], not [min]
  H_int<-cumtrapz(dsc_intseg$Time*60,-dsc_intseg$HeatFlow_0)
  S_int<-cumtrapz(dsc_intseg$Time*60,-1000*dsc_intseg$HeatFlow_0/dsc_intseg$TempK)
  DH<-H_int[length(H_int)]
  DS<-S_int[length(S_int)]
  
  dsc_df<-cbind(dsc_df,dtdt,nom_HF)
  dsc_df$H<-0
  dsc_df$S<-0
  dsc_df$H[intsegs]<-H_int
  dsc_df$S[intsegs]<-S_int
  dsc_df$H[max(pt1,pt2):nrow(dsc_df)]<-DH
  dsc_df$S[max(pt1,pt2):nrow(dsc_df)]<-DS
  
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


function(input, output, session) {
  
  # 1. アップロードされたファイルのパース処理 (tbを含む)
  parsed_data <- reactive({
    req(input$file)
    file_path <- input$file$datapath
    
    alldata <- readLines(file(file_path, encoding = "UTF-16LE"))
    
    # サンプル・実験情報
    expinfoid <- charmatch(c("#IDENTITY", "#DATE", "#SAMPLE:", "#SAMPLE MASS"), alldata)
    expinfo_vals <- alldata[expinfoid] |>
      strsplit(split = ",") |>
      unlist() |>
      _[seq(2, 8, by = 2)] |>
      gsub(" ", "", x = _)
    
    expinfo_df <- data.frame(
      Terms = c("Sample No.", "Meas. Date", "Chem. Species", "Sample Mass (mg)"),
      Exp.info = expinfo_vals,
      stringsAsFactors = FALSE
    )
    
    # セグメント情報
    segs <- which(startsWith(alldata, "#SEG."))
    lsegs <- length(segs)
    
    seglst <- alldata[segs] |>
      gsub("/", ",", x = _) |>
      gsub("°C", "", x = _) |>
      gsub("K,", "K/", x = _) |>
      strsplit(split = ",")
    
    segdata <- data.frame(matrix("", nrow = lsegs, ncol = 4))
    colnames(segdata) <- c("from", "to", "expl", "mode")
    
    for (i0 in seq_len(lsegs)) {
      segdata[i0, 1] <- as.numeric(seglst[[i0]][2])
      segdata[i0, 2] <- as.numeric(seglst[[i0]][4])
      segdata[i0, 3] <- seglst[[i0]][3]
      
      from_val <- as.numeric(segdata[i0, 1])
      to_val   <- as.numeric(segdata[i0, 2])
      
      if (!is.na(from_val) && !is.na(to_val)) {
        if (from_val > to_val) {
          segdata[i0, 4] <- "cooling"
        } else if (from_val < to_val) {
          segdata[i0, 4] <- "heating"
        } else {
          segdata[i0, 4] <- "isothermal"
        }
      } else {
        segdata[i0, 4] <- "unknown"
      }
    }
    
    # --- functions_1.R 相当の tb (数値データ) 作成処理 ---
    starttb <- charmatch("##", alldata)
    endtb   <- length(alldata)
    
    # ヘッダー行とデータ行の分離
    raw_tb <- alldata[(starttb + 1):endtb]
    # 空行等の除去
    raw_tb <- raw_tb[raw_tb != ""]
    
    # カンマ区切りデータフレーム化
    tb_list <- strsplit(raw_tb, split = ",")
    tb <- as.data.frame(do.call(rbind, tb_list), stringsAsFactors = FALSE)
    
    # 列名の設定（NETZSCHの標準出力構成: Temp, Time, HeatFlow, Segment等に対応）
    # ※列数が異なる場合は数値型変換を行いつつ整形
    for (col_idx in seq_len(ncol(tb))) {
      tb[[col_idx]] <- as.numeric(trimws(tb[[col_idx]]))
    }
    
    # 通常のDSC構成: 1列目:Temp(℃), 2列目:Time(min), 3列目:HeatFlow(mW/mg), ..., 最終列:Segment番号
    # 既存コードで `seg` 列が存在することを想定
    if (ncol(tb) >= 5) {
      colnames(tb)[1:5] <- c("Temp", "Time", "HeatFlow", "DSC_sens", "seg")
    }
    
    #絶対温度の列を挿入、順番は時間、セ氏度、絶対温度、熱流、電圧、セグメント
    TempK<-tb$Temp+273.15
    tb<-cbind(tb[,c(2,1)],TempK,tb[,3:5])
    
    sep_tb<-list()
    for(i1 in seq_len(lsegs)){
      sep_tb[[i1]]<-tb[tb$seg==i1,]
    }
    
    list(
      expinfo = expinfo_df,
      segdata = segdata,
      tb = sep_tb
    )
  })
  
  # 実験情報表示
  output$exp_info_table <- renderTable({
    req(parsed_data())
    parsed_data()$expinfo
  })
  
  # segdata表示
  output$segdata_table <- renderTable({
    req(parsed_data())
    segdf <- parsed_data()$segdata
    cbind(Segment = seq_len(nrow(segdf)), segdf)
  })
  
  # セグメント選択用チェックボックス (heating / cooling のみ)
  output$segment_select_ui <- renderUI({
    req(parsed_data())
    segdf <- parsed_data()$segdata
    segdf$orig_id <- seq_len(nrow(segdf))
    
    filtered_segdf <- segdf[segdf$mode %in% c("heating", "cooling"), ]
    
    choices_vec <- setNames(
      filtered_segdf$orig_id,
      paste0("SEG ", filtered_segdf$orig_id, ": ", filtered_segdf$from, "℃ → ", filtered_segdf$to, "℃ (", filtered_segdf$mode, ")")
    )
    
    checkboxGroupInput(
      inputId  = "selected_segments",
      label    = "解析対象セグメントの選択 (heating/cooling):",
      choices  = choices_vec,
      selected = NULL
    )
  })
  
  # 状態保持用 (reactiveValues)
  state <- reactiveValues(
    current_tb = NULL,
    active_seg = NULL,
    out_list = list(),
    ndata = 0
  )
  
  # ファイル読み込み時に state$current_tb を初期化
  observe({
    req(parsed_data())
    state$current_tb <- parsed_data()$tb
  })
  
  # 「選択したセグメントを表示・解析」ボタン押下時
  observeEvent(input$run_btn, {
    req(input$selected_segments)
    # 選択された最初のセグメントをデフォルトアクティブにする
    state$active_seg <- as.numeric(input$selected_segments[1])
  })
  
  # 表示対象セグメント切り替えUI
  output$segment_tab_selector <- renderUI({
    req(input$selected_segments)
    selectInput(
      "active_seg_select",
      "プロット・補正対象のセグメント切替:",
      choices = input$selected_segments,
      selected = state$active_seg
    )
  })
  
  observeEvent(input$active_seg_select, {
    state$active_seg <- as.numeric(input$active_seg_select)
  })
  
  #ベースラインフィッティング関数のチョイス
  output$baseline_method<-renderUI({
    req(state$current_tb, state$active_seg)
    selectInput(
      "baseline_method",
      "ベースラインフィッティングの方法",
      choices = c("2 line 2nd connection",
                  "line",
                  "poly(2nd)",
                  "poly(3rd)",
                  "poly(4th)",
                  "2 points"),
      selected = "2 line 2nd connection"
    )
  })
  
  # 選択されたセグメントの温度範囲スライダーUI生成
  output$baseline_controls_ui <- renderUI({
    req(state$current_tb, state$active_seg)
    
    sub_tb <- state$current_tb[[state$active_seg]]
    req(nrow(sub_tb) > 0)
    
    min_t <- min(sub_tb$TempK, na.rm = TRUE)
    max_t <- max(sub_tb$TempK, na.rm = TRUE)
    space_t<-(max_t-min_t)/5
    
    tagList(
      h4("ベースライン範囲指定"),
      sliderInput(
        "bg_temp_range",
        "ピーク温度範囲",
        min = floor(min_t),
        max = ceiling(max_t),
        value = c(floor(min_t+2*space_t), ceiling(max_t-2*space_t)),
        step = 0.5
      ),
      sliderInput(
        "fit_temp_range",
        "ベースライン温度範囲",
        min = floor(min_t),
        max = ceiling(max_t),
        value = c(floor(min_t+space_t), ceiling(max_t-space_t)),
        step = 0.5
      )
    )
  })
  
  # 1. Heat Flow - T プロット & ベースラインプレビュー描画 (functions_2.R 準拠)
  output$dsc_plot <- renderPlot({
    req(state$current_tb, state$active_seg)
    
    targ_seg <- state$active_seg
    targ_tb  <- state$current_tb[[targ_seg]]
    
    req(nrow(targ_tb) > 0)
    
    p <- ggplot(targ_tb, aes(x = TempK, y = HeatFlow)) +
      geom_line(color = "black", linewidth = 1) +
      labs(x = "Temperature (K)", y = "Heat Flow (mW/mg)", title = paste("Segment", targ_seg)) +
      mirror_x + mirror_y + mytheme
    
    # スライダーが設定されている場合、ベースライン直線をプレビュー描画
    
    #ベースラインフィット方式別
    if (!is.null(input$bg_temp_range)) {
      t1 <- input$bg_temp_range[1]
      t2 <- input$bg_temp_range[2]
      t3 <- input$fit_temp_range[1]
      t4 <- input$fit_temp_range[2]

      blresult<-dscbl(targ_tb,input$baseline_method,t1,t2,t3,t4)
      p <- p +
        geom_vline(xintercept = blresult$temp_sec[c(1,2)], color = "red",linetype="dashed") +
        geom_vline(xintercept = blresult$temp_sec[c(3,4)], color = "blue",linetype="dashed") +
        geom_line(data=blresult$dsc_df, aes(x = TempK, y = basel),
                     color = "red", linewidth = 1)
    }
    
    p
  },
  res = 120)
  
  # 2 & 3. ベースライン引き算処理 & tbへの列追加 (cbind)
  observeEvent(input$apply_sub_btn, {
    req(state$current_tb, state$active_seg, input$bg_temp_range)
    
    targ_seg <- state$active_seg
    t_df     <- state$current_tb
    state$ndata <- state$ndata+1
    
    # 対象セグメントの抽出
    targ_tb <- t_df[[targ_seg]]
    
    t1 <- input$bg_temp_range[1]
    t2 <- input$bg_temp_range[2]
    t3 <- input$fit_temp_range[1]
    t4 <- input$fit_temp_range[2]
    
    res_list<-dscbl(targ_tb,input$baseline_method,t1,t2,t3,t4)
    res_list$seg<-targ_seg
    #ここにエントロピー計算を加える
    res_list<-int_dsc(res_list)
    
    state$out_list[[state$ndata]]<-res_list
   
    output$status_msg <- renderText({
      paste0("Segment ", targ_seg, " に対するベースライン引き算処理が完了し、`tb` に 'HF_sub' 列が追加/更新されました。\n",
             "データ数:",state$ndata)

    })
  })
  
  #プレビューする解析番号
  output$prev_choice<-renderUI({
    req(state$out_list)
    selectInput(
      "prev_choice",
      "プレビュー(解析番号)",
      choices = as.character(1:state$ndata),
      selected = "1"
    )
  })
  
  # プレビュー
  #return(list(method=method,method_id=method_id,dsc_df=res_df,temp_sec=temp_sec,coefs=params))
  output$prev_txt <- renderText({
    req(state$out_list,input$prev_choice)
    prevlst<-state$out_list[[as.numeric(input$prev_choice)]]
    paste0("Segment No.: ",prevlst$seg,"\n",
           "Baseline method: ",prevlst$method,"\n",
           "Peak section: ",prevlst$temp_sec[1]," K ~ ",prevlst$temp_sec[2]," K\n",
           "DH: ",round(prevlst$DH,2),"kJ kg-1, DS: ",round(prevlst$DS,2),"J kg-1 K-1")
  })
  
  output$prev_dtdt <- renderPlot({
    req(state$out_list,input$prev_choice)
    prevlst<-state$out_list[[as.numeric(input$prev_choice)]]
    p <- ggplot(prevlst$dsc_df, aes(x=Time,y=dtdt)) +
      geom_line(color = "black", linewidth = 1) +
      labs(x = "Meas. Time (min)", y = expression(Heating~rate~(K~min^{-1}))) +
      mirror_x + mirror_y + mytheme
    p
  },
  res = 120)
  
  output$prev_main <- renderPlot({
    req(state$out_list,input$prev_choice)
    prevdf<-state$out_list[[as.numeric(input$prev_choice)]]$dsc_df
    
    p_HF<-ggplot(prevdf, aes(x=TempK,y=HeatFlow_0)) +
      geom_line(color = "red", linewidth = 1) +
      labs(x = "", y = expression(Heat~Flow~(W~g^{-1}))) +
      mirror_x + mirror_y + mytheme+
      theme(
        axis.text.x.bottom = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin()
      )
    
    p_H<-ggplot(prevdf, aes(x=TempK,y=H)) +
      geom_line(color = "red", linewidth = 1) +
      labs(x = "", y = expression(italic(H)~(kJ~kg^{-1}))) +
      mirror_x + mirror_y + mytheme +
      theme(
        axis.text.x.bottom = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin()
      )
    
    p_S<-ggplot(prevdf, aes(x=TempK,y=S)) +
      geom_line(color = "red", linewidth = 1) +
      labs(x = "Temperature (K)", y = expression(italic(S)~(J~K^{-1}~kg^{-1}))) +
      mirror_x + mirror_y + mytheme+
      theme(
        plot.margin = margin()
      )
    
    p_HF/p_H/p_S
  },
  res = 120)
  
  
  output$download_csv <- downloadHandler(
    filename = function(){
      expinfo<-parsed_data()$expinfo[,2]
      segnum<-state$out_list[[as.numeric(input$prev_choice)]]$seg
      paste0(expinfo[3],"_",expinfo[1],"_SEG",segnum,"_",as.numeric(input$prev_choice),".csv")
    },
    content = function(file) {
      data<-state$out_list[[as.numeric(input$prev_choice)]]$dsc_df
      write.csv(data, file,row.names = F)
    }
  )
  
  output$download_report <- downloadHandler(
    filename = function(){
      expinfo<-parsed_data()$expinfo[,2]
      paste0(expinfo[3],"_",expinfo[1],"_report.html")
    },
    content = function(file){
      
    }
  )
}