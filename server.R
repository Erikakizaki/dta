library(shiny)
library(ggplot2)
library(tidyr)
library(patchwork)



# ggplot2用の共通テーマ・軸設定 (functions_2.R より)

library(pracma)
source("functions.R")


function(input, output, session) {
  
  fileinfos <- reactive({
    req(input$file,input$resist)
    resist<-if(input$resist=="Yes"){T}else{F}
    nf <- nrow(input$file)
    rowfiles<-list()
    spldatas<-list()
    rsdatas<-list()
    
    for(i01 in 1:nf){
      rowfiles[[i01]]<-read.csv(input$file$datapath[i01],skip = 8)
      spldatas[[i01]]<-splitData(rowfiles[[i01]],resist = resist)
      rsdatas[[i01]]<-rsmplData(spldatas[[i01]],resist = resist)
    }
    
    list(
      file = input$file,
      nf = nf,
      rowfiles = rowfiles,
      spldatas = spldatas,
      rsdatas = rsdatas
      )
    })
  
  
  output$exp_info_table <- renderTable({
    req(fileinfos)
    fileinfos()$file
  })
  output$nf<-renderText({
    req(fileinfos)
    paste0("The number of loaded files is ",
           fileinfos()$nf
           )
    
  })
  

#smooth
  
  
  sgdatas<- reactive({
    req(fileinfos,input$sgwindow)
    sgdts<-list()
    sgwindow<-input$sgwindow
    rsdatas<-fileinfos()$rsdatas
    nf<-fileinfos()$nf
    resist<-if(input$resist=="Yes"){T}else{F}
    
    for(i02 in 1:nf){
      sgdts[[i02]]<-sgsmoothdata(rsdatas[[i02]],sgwindow = sgwindow, resist = resist)
    }
    
    sgdts
  })
  
   output$testtable<-renderTable({
     req(sgdatas)
     sgdatas()[[1]]
   })
  
  output$smoothplot<-renderPlot({
    req(sgdatas,fileinfos,input$sgwindow)
    sgwindow<-input$sgwindow
    sgdts<-sgdatas()
    #rsdts<-as.data.frame(fileinfos()$rsdatas)
    nf<-fileinfos()$nf
    pdata<-sgdts[[1]]
    #print(pdata)
    testplot<-ggplot(pdata,aes(x=Temp_re,y=hf))+
      geom_line()+
      labs(x="Temperature (K)",y="Heat flow (arb.unit)")+
      mytheme + mirror_x +mirror_y
    
    testplot
  })
  
  
#Baseline
  #UI作成
  output$data_analysis_ui <- renderUI({
    req(fileinfos(), sgdatas())
    files <- fileinfos()$file$name
    
    if (length(files) == 0) return(NULL)
    
    # ファイルごとのタブを作成
    tabs <- lapply(seq_along(files), function(i) {
      fname <- files[i]
      df <- sgdatas()[[i]]
      
      # Temp_re の範囲を取得してスライダーの初期範囲を設定
      min_t <- min(df$Temp_re, na.rm = TRUE)
      max_t <- max(df$Temp_re, na.rm = TRUE)
      
      tabPanel(
        title = fname,
        br(),
        fluidRow(
          colWidth = 4,
          wellPanel(
            h5("ピーク範囲指定 (ベースライン)"),
            
            
            # 各ファイル固有のIDを付与 (例: peak_range_1, peak_range_2...)
            sliderInput(
              inputId = paste0("peak_range_", i),
              label = "ピーク温度範囲:",
              min = floor(min_t),
              max = ceiling(max_t),
              value = c(floor(min_t + (max_t - min_t)*0.4), 
                        ceiling(min_t + (max_t - min_t)*0.6)),
              step = 0.5
            ),
            sliderInput(
              inputId = paste0("fit_range_", i),
              "ベースライン温度範囲",
              min = floor(min_t),
              max = ceiling(max_t),
              value = c(floor(min_t + (max_t - min_t)*0.2), 
                        ceiling(min_t + (max_t - min_t)*0.8)),
              step = 0.5
            )
          )
        ),
        fluidRow(
          colWidth = 12,
          plotOutput(outputId = paste0("dta_plot_", i), height = "400px")
        )
      )
    })
    
    do.call(tabsetPanel, tabs)
  })
  
  # ファイルのデータ処理・プレビュー描画
  observe({
    req(fileinfos(), sgdatas() )
    files <- fileinfos()$file$name
    
    lapply(seq_along(files), function(i) {
      plot_id <- paste0("dta_plot_", i)
      range_id <- paste0("peak_range_", i)
      fitrng_id <- paste0("fit_range_", i)
      
      output[[plot_id]] <- renderPlot({
        df <- sgdatas()[[i]]
        rng <- input[[range_id]]
        fitrng <- input[[fitrng_id]]
        
        
        # 基本のプロット (hf vs. Temp_re)
        p <- ggplot(df, aes(x = Temp_re, y = hf)) +
          geom_line(color = "black") +
          labs(
            title = paste("データ:", files[i]),
            x = "Temperature (K)",
            y = "Heat Flow (arb. unit)"
          ) +
          mytheme + mirror_x + mirror_y
        
        # スライダー入力値が存在する場合、ベースラインのプレビューを描画
        if (!is.null(rng)) {
          t1 <- rng[1]
          t2 <- rng[2]
          t3 <- fitrng[1]
          t4 <- fitrng[2]
          
          blresult <- dtabl(df,input$baseline_method,t1,t2,t3,t4)



          # ベースライン直線のプロット (赤点と赤破線)
          p <- p +
            geom_vline(xintercept = blresult$temp_sec[c(1,2)], color = "red",linetype="dashed") +
            geom_vline(xintercept = blresult$temp_sec[c(3,4)], color = "blue",linetype="dashed") +
            geom_line(data=blresult$dta_df, aes(x = Temp_re, y = basel),
                      color = "red", linewidth = 1)
        }
        
        return(p)
      })
    })
  })
  
  #データ保存
  
  
  out_list<-reactiveValues()
  observeEvent(input$apply_sub_btn, {
    req(fileinfos(), sgdatas())
    
    files <- fileinfos()$file$name
    
    lapply(seq_along(files), function(i) {
      data_id <- paste0("data_", i)
      range_id <- paste0("peak_range_", i)
      fitrng_id <- paste0("fit_range_", i)
      
      df <- sgdatas()[[i]]
      rng <- input[[range_id]]
      fitrng <- input[[fitrng_id]]
      if (!is.null(rng)) {
        t1 <- rng[1]
        t2 <- rng[2]
        t3 <- fitrng[1]
        t4 <- fitrng[2]
        out_list[[data_id]] <- dtabl(df,input$baseline_method,t1,t2,t3,t4)
      }
    })
    output$status_msg <- renderText({
      paste0(" ベースライン引き算処理が完了しました。")
    })
  })
  
  #基準データ選択、圧力入力用UI
  output$stddata <- renderUI({
    req(fileinfos())
    files <- fileinfos()$file$name
    
    
    if (length(files) == 0) return(NULL)
    grid_content <- unlist(lapply(seq_along(files), function(i) {
      list(
        tags$label(paste0(files[i], ":")),
        textInput(paste0("press_", i), label = NULL, value = "0")
      )
    }), recursive = FALSE)
    
    list(
      selectInput(
        "stdfile",
        "エントロピー較正基準データ",
        choices = files
      ),
      textInput("stdS","基準エントロピー(J/K/kg)"),
      br(),
      h5("圧力計読み値(MPa)"),
      # CSS Grid の div で囲む
      do.call(div, c(list(class = "grid-container"), grid_content)),
      actionButton("apply_stddata_btn", "適用", class = "btn-success")
    )
  })
  
  

  
  
  #変換定数計算
  
  expsetting<-reactiveValues()
  observeEvent(input$apply_stddata_btn, {
    req(fileinfos(),out_list)
    
    
    expsetting$stdfile<-input$stdfile
    expsetting$stdS<-as.numeric(input$stdS)
    
    files <- fileinfos()$file$name
    
    if(length(names(out_list))==0){
      output$status_msg2 <- renderText({
        paste0("先にベースライン引きを実行してください")
      })
    }else{
      stdid<-which(files==expsetting$stdfile)
      Scoef<-calcoef(out_list[[paste0("data_",stdid)]],expsetting$stdS)
      
      expsetting$Scoef<-Scoef
      
      lapply(seq_along(files), function(i) {
        press_id<-paste0("press_", i)
        
        expsetting[[press_id]]<-input[[press_id]]
      })
      
      output$status_msg2 <- renderText({
        paste0("実験データを設定しました; ",expsetting$stdfile,", coefficient; ",Scoef)
      })
    }
    
    
  })
  
  #エントロピー計算・まとめ
  observeEvent(input$summary_btn, {
    req(fileinfos(), out_list, expsetting)
    files <- fileinfos()$file$name
    Scoef <- expsetting$Scoef
    # -------------------------------------------------------------
    # 【1】データの計算および out_list の更新処理（renderUI の外で実施）
    # -------------------------------------------------------------
    for (i in seq_along(files)) {
      local({
        my_i <- i
        data_id <- paste0("data_", my_i)
        plot_id <- paste0("plotf_", my_i)
        
        # 既存のデータを計算して上書き
        res <- int_dta(out_list[[data_id]], Scoef)
        out_list[[data_id]] <- res
        
        # プロットを作成
        p <- ggplot(res$dta_df, aes(x = Temp_re, y = nom_HF)) +
          geom_line(color = "black") +
          labs(
            x = "Temperature (K)",
            y = "Heat Flow (arb. unit)"
          ) +
          mytheme + mirror_x + mirror_y
        
        # local 内なので my_i / p の値が固定される
        output[[plot_id]] <- renderPlot({ p })
      })
      
    }
    
    # -------------------------------------------------------------
    # 【2】UIの描画処理（データの更新・書き込みは絶対に行わない）
    # -------------------------------------------------------------
    output$summary <- renderUI({
      n_reports <- length(files)
      
      # まだ解析が行われていない場合
      if (n_reports == 0) {
        return(p("まだ解析結果がありません。解析を実行してください。"))
      }
      
      # 1〜n回分の要素を lapply で動的に作成する
      report_elements <- lapply(seq_along(files), function(i) {
        data_id <- paste0("data_", i)
        
        # ★ここでは out_list の「参照（読み込み）」のみを行う
        res <- out_list[[data_id]]
        plot_id <- paste0("plotf_", i)
        
        
        
        div(
          style = "border: 1px solid #ccc; padding: 15px; margin-bottom: 20px; border-radius: 5px;",
          h4(paste0("解析回数: #", i)),
          
          # ピーク範囲およびエントロピー変化のテキスト表示
          p(strong("ピーク範囲: "), paste0(res$temp_sec[1], " K ~ ", res$temp_sec[2], " K")),
          p(strong("エントロピー変化: "), paste0(res$DS, " J K-1 kg-1")),
          br(),
          # グラフの描画用 output（必要に応じて有効化）
          plotOutput(plot_id)
        )
      })
      
      # リストをタグリストに変換して返す
      do.call(tagList, report_elements)
    })
  })
  

}