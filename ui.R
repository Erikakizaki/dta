library(shiny)
library(ggplot2)
library(tidyr)
library(patchwork)

fluidPage(
  titlePanel("DTA"),
  
  navlistPanel(
    tabPanel("データ・設定",
             h4("実験データ読込"),
             fileInput("file", "CSVファイルをアップロードしてください",
                       accept = c(".csv", ".txt"),multiple = T),
             radioButtons("resist","抵抗値測定",c("Yes","No"),"No"),
             tableOutput("exp_info_table"),textOutput("nf")
             ),
    tabPanel("補間・平滑化",
             h4("補間・平滑化"),
             sliderInput("sgwindow","SGフィルタの幅(1ならばフィルタなし)",1,51,31,2),
             #tableOutput("testtable"),
             plotOutput("smoothplot")
             ),
    tabPanel(
      "ベースライン",
      h4("ベースラインフィッティング"),
      br(),
      selectInput(
        "baseline_method",
        "ベースラインフィッティングの方法",
        choices = c(
          "line",
          "poly(2nd)",
          "poly(3rd)",
          "poly(4th)",
          "2 points",
          "2 line 2nd connection"
        ),
        selected = "line"
      ),
      actionButton("apply_sub_btn", "すべてのベースライン引き算を適用", class = "btn-success"),
      br(),
      verbatimTextOutput("status_msg"),
      br(),
      br(),
      uiOutput("data_analysis_ui")
    ),
    tabPanel("実験情報",
             h4("基準データ"),
             uiOutput("stddata")
             )


  )

#Interface  

  
   
)

