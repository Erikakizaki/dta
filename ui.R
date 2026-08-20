library(shiny)
library(ggplot2)
library(tidyr)
library(patchwork)

fluidPage(
  titlePanel("DSC 解析レポート"),

#Interface    
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "CSVファイルをアップロードしてください",
                accept = c(".csv", ".txt")),
      hr(),
      uiOutput("segment_select_ui"),
      br(),
      actionButton("run_btn", "解析するセグメントを選択", class = "btn-primary"),
      
      hr(),
      # ベースライン調整用UI（セグメントが選択されたら動的表示）
      uiOutput("baseline_method"),
      uiOutput("baseline_controls_ui")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("データ概要",
                 h4("実験・サンプル情報"),
                 tableOutput("exp_info_table"),
                 h4("セグメントデータ一覧 (`segdata`)"),
                 tableOutput("segdata_table")
        ),
        
        tabPanel("Heat Flow - T プロット & ベースライン補正",
                 br(),
                 # セグメント切替用のUI
                 uiOutput("segment_tab_selector"),
                 plotOutput("dsc_plot", height = "450px"),
                 hr(),
                 h4("ベースライン補正実行"),
                 actionButton("apply_sub_btn", "ベースライン引き算を適用してデータ(tb)に追加", class = "btn-success"),
                 br(), br(),
                 verbatimTextOutput("status_msg")
        ),
        
        tabPanel("レポート",
                 br(),
                 uiOutput("prev_choice"),
                 downloadButton("download_csv", "csvをダウンロード"),
                 verbatimTextOutput("prev_txt"),
                 plotOutput("prev_dtdt",width = "640px"),
                 br(),
                 plotOutput("prev_main",height = "900px",width = "640px")
        )
      )
    )
  )
)

