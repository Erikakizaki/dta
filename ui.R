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
#             sliderInput("sgwindow","SGフィルタの幅(1ならばフィルタなし)",1,51,31,2),
#             hr(style = "border-top: 1px solid #000000;"),
#             actionButton("cleardata","データをクリア"),

             
             tableOutput("exp_info_table"),
textOutput("nf")
             
    ),
  tabPanel("補間・平滑化",
         h4("補間・平滑化"),
         
         sliderInput("sgwindow","SGフィルタの幅(1ならばフィルタなし)",1,51,31,2),
         
         #tableOutput("testtable"),
         plotOutput("smoothplot")
         

         
),
tabPanel("ベースライン",
         h4("ベースラインフィッティング"),
         
         )

  )

#Interface  

  
   
)

