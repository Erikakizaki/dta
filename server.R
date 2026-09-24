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
    testplot<-ggplot(pdata[sgwindow:(nrow(pdata)-sgwindow),],aes(x=Temp_re,y=hf))+
      geom_line()+
      labs(x="Temperature (K)",y="Heat flow (arb.unit)")+
      mytheme + mirror_x +mirror_y
    
    testplot
  })
  
  
#Baseline
  
  

}