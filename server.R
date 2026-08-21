library(shiny)
library(ggplot2)
library(tidyr)
library(patchwork)



# ggplot2用の共通テーマ・軸設定 (functions_2.R より)

library(pracma)



function(input, output, session) {
  
  fileinfos <- reactive({
    req(input$file)
    input$file
  })
  
  output$exp_info_table <- renderTable({
    req(fileinfos())
    fileinfos()
  })
  

}