###################################################################################
# Shiny app for the diurnal transcriptome response to temperature stress in Sorghum
###################################################################################


# Packages ----
library(shiny)
library(DT)
library(shinymanager)
library(shinydashboard)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(scales)
library(readr)
library(leaflet)
library(gridExtra)
library(svglite)
library(png)
library(tidyr)
library(doBy)
library(reshape2)
library(stringr)
library(tidyverse)
library(shinythemes)

# Loading data ----
#setwd("C:/Users/tbonnot/Documents/sorghum/Shiny app/V1")

#========= Data Sorghum ========
#===============================
rlog <- read.table("rlog_expression_values.txt", header = T, stringsAsFactors = T)

rlog1 <- rlog[,c(1,50:length(rlog))]
rlog2 <- rlog[,c(1:49)]

#annot <- read.csv2("Sbicolor_454_v3.1.1.annotation_info_shiny.csv", header = T, sep = ";")
annot <- read.csv2(text = readLines("Sbicolor_454_v3.1.1.annotation_info_shiny.csv", warn = FALSE),header=T)
#annot <- read.table("Sbicolor_454_v3.1.1.annotation_info_shiny.txt")

Lai <- read.table("Lai_data.txt", header = T)

Genes <- read.table("Genes.txt", header = T)
Genes2 <- read.table("Genes2.txt", header = T)

#-----------------------------------------------------------------------------------------------#
# ui.R ---- contains information about the layout of the app as it appears in web browser ------
#-----------------------------------------------------------------------------------------------#

ui <- fluidPage(

  theme = shinytheme("flatly"),
  
  tags$head(tags$style(HTML('
  
        sidebar {background-color: #FF66CC;font-size: 14px}
        /* logo */
        .skin-blue .main-header .logo {
        background-color: #7f00caff;
        }
 
        /* navbar (rest of the header) */
        .skin-blue .main-header .navbar {
        background-color: #7f00caff;
        }
        
        .navbar-custom-menu {
       float: left!important;
        }
        
         .tabbable > .nav > li > a {background-color: #7300b6ff; color:white;font-size: 15px}

        # .shiny-notification {position: fixed; top: 20% ;left: 10%; right: 80%}
        #  .form-group {
        #      margin-bottom: 0 !important;
            }
        
        body, label, input, button, select {
          font-family: "Arial";
        }')
  )),
  dashboardPage(
    dashboardHeader(title= "Diurnal transcriptomic responses to cold and heat stress in Sorghum", titleWidth = 1000,
                    tags$li(class = "dropdown")),
    

    dashboardSidebar(disable=T),
    dashboardBody(
      tabsetPanel(type="tabs",
                  tabPanel(title = "Sorghum genes", 
                           sidebarLayout(
                             sidebarPanel(width = 2,
                                          style = "margin-bottom:15px"),
                             mainPanel(
                               fluidRow(box(width=20,title = span("Gene descriptions and best Arabidopsis and Rice hits", style = "color: #7f00caff; font-weight: bold; font-size: 20px"),
                                            dataTableOutput("Three"),
                                            htmlOutput("Legend_table")))
                             )
                           )
                  ),
                  
                  tabPanel(title = "Cold stress response", 
                           sidebarLayout(
                             sidebarPanel(selectizeInput(inputId = "AGI", choices = NULL, label = "Select a locus ID:"),
                                          width = 3, style = "margin-bottom:15px"),
                             mainPanel(
                               fluidRow(box(width=20,title = span("Diurnal transcriptomic response to cold stress", style = "color: #7f00caff; font-weight: bold; font-size: 20px"),
                                                    plotOutput("One", height="330px"), status="primary", solidHeader = F, collapsible = F,
                                                    div(style = "margin-top:20px"),        
                                                    htmlOutput("Legend_cold"),
                                                    div(style = "margin-top:20px"),
                                                    downloadButton("Diurnal_cold_stressPNG","Plot.png"),
                                                    downloadButton("Diurnal_cold_stressSVG","Plot.svg"),
                                                    downloadButton("Diurnal_cold_stressPDF","Plot.pdf"),
                                                    downloadButton("Diurnal_cold_stress", "Data.txt")))
                             )
                           )
                  ),
                  
                  tabPanel(title = "Heat stress response", 
                           sidebarLayout(
                             sidebarPanel(selectizeInput(inputId = "AGI2", choices = NULL, label = "Select a locus ID:"),
                                          width = 3, style = "margin-bottom:15px"),
                             mainPanel(
                               fluidRow(box(width=20,title = span("Diurnal transcriptomic response to heat stress", style = "color: #7f00caff; font-weight: bold; font-size: 20px"),
                                            plotOutput("Two", height="330px"), status="primary", solidHeader = F, collapsible = F,
                                            div(style = "margin-top:20px"),
                                            htmlOutput("Legend_heat"),
                                            div(style = "margin-top:20px"),   
                                            downloadButton("Diurnal_heat_stressPNG","Plot.png"),
                                            downloadButton("Diurnal_heat_stressSVG","Plot.svg"),
                                            downloadButton("Diurnal_heat_stressPDF","Plot.pdf"),
                                            downloadButton("Diurnal_heat_stress", "Data.txt")))
                             )
                           )
                  ),
                  
                  tabPanel(title = "Diurnal rhythmic transcriptome", 
                           sidebarLayout(
                             sidebarPanel(selectizeInput(inputId = "AGI3", choices = NULL, label = "Select a locus ID:"),
                                          width = 3, style = "margin-bottom:15px",
                                          div(HTML("<i>Only rhythmic genes can be selected</i>"),style = "margin-bottom:15px")),
                             mainPanel(
                               fluidRow(box(width=20,title = span("Transcript profile of diurnal rhythmic genes from Lai et al. (2020)", style = "color: #7f00caff; font-weight: bold; font-size: 20px"),
                                            plotOutput("Four"), status="primary", solidHeader = F, collapsible = F,
                                            div(style = "margin-top:20px"),
                                            htmlOutput("Legend_diurnal"),
                                            div(style = "margin-top:20px"),
                                            downloadButton("Diurnal_LaiPNG","Plot.png"),
                                            downloadButton("Diurnal_LaiSVG","Plot.svg"),
                                            downloadButton("Diurnal_LaiPDF","Plot.pdf"),
                                            downloadButton("Diurnal_Lai", "Data.txt")))
                             )
                           )
                  )
                  
      )
    )
  )
)

  


# server.R ---- contains information about the computation of the app, creating plots, tables, maps etc. using information provided by the user
server <- function(input, output, session) {
  
  updateSelectizeInput(session, inputId = "AGI", choices = Genes$AGI, server = TRUE)
  updateSelectizeInput(session, inputId = "AGI2", choices = Genes$AGI, server = TRUE)
  updateSelectizeInput(session, inputId = "AGI3", choices = Genes2$AGI, server = TRUE)
  
  # Select data
  #------------
  data1 <- eventReactive(input$AGI, {
    data1 <- rlog2[rlog2$AGI == input$AGI,]
    data1 <- melt(data1, id.vars = "AGI")
    data1 <- data1 %>% separate(variable, into = c("Exp","Genotype", "Time", "Response"))
    data1$Genotype <- str_replace_all(data1$Genotype, 'RTX', 'RTX430')
    data1$Time <- str_replace_all(data1$Time, 'T', '')
    data1$Response <- str_replace_all(data1$Response, 'Cont1', '30.1')
    data1$Response <- str_replace_all(data1$Response, 'Cont2', '30.2')
    data1$Response <- str_replace_all(data1$Response, 'Cont3', '30.3')
    data1$Response <- str_replace_all(data1$Response, 'Cold1', '10.1')
    data1$Response <- str_replace_all(data1$Response, 'Cold2', '10.2')
    data1$Response <- str_replace_all(data1$Response, 'Cold3', '10.3')
    data1$Response <- str_replace_all(data1$Response, 'Cold4', '10.3')
    data1 <- data1 %>% separate(Response, into = c("Temperature", "Replicate"))
  })
  
  data_mean1 <- reactive({
    data_mean1 <- data1()
    data_mean1 <- summaryBy(.~AGI+Genotype+Time+Temperature, keep.names = T,FUN = function(x){mean(x,na.rm=T)} , data= data_mean1)
  })
  
  data_sd1 <- reactive({
    data_sd1 <- data1()
    data_sd1 <- aggregate(value~ AGI+Genotype+Time+Temperature, data = data_sd1, FUN = sd)
  })
  
  data_plot1 <- reactive({
    data_mean1 <- data_mean1()
    data_sd1 <- data_sd1()
    data_plot1 <- merge.data.frame(data_mean1,data_sd1, by = c("AGI","Genotype","Time","Temperature"))
  })

  
  data2 <- eventReactive(input$AGI2, {
    data2 <- rlog1[rlog1$AGI == input$AGI2,]
    data2 <- melt(data2, id.vars = "AGI")
    data2 <- data2 %>% separate(variable, into = c("Exp","Genotype", "Time", "Response"))
    data2$Genotype <- str_replace_all(data2$Genotype, 'RTX', 'RTX430')
    data2$Time <- str_replace_all(data2$Time, 'T', '')
    data2$Response <- str_replace_all(data2$Response, 'C1', '30.1')
    data2$Response <- str_replace_all(data2$Response, 'C2', '30.2')
    data2$Response <- str_replace_all(data2$Response, 'C3', '30.3')
    data2$Response <- str_replace_all(data2$Response, 'H1', '42.1')
    data2$Response <- str_replace_all(data2$Response, 'H2', '42.2')
    data2$Response <- str_replace_all(data2$Response, 'H3', '42.3')
    data2 <- data2 %>% separate(Response, into = c("Temperature", "Replicate"))
  })
  
  data_mean2 <- reactive({
    data_mean2 <- data2()
    data_mean2 <- summaryBy(.~AGI+Genotype+Time+Temperature, keep.names = T,FUN = function(x){mean(x,na.rm=T)} , data= data_mean2)
  })
  
  data_sd2 <- reactive({
    data_sd2 <- data2()
    data_sd2 <- aggregate(value~ AGI+Genotype+Time+Temperature, data = data_sd2, FUN = sd)
  })
  
  data_plot2 <- reactive({
    data_mean2 <- data_mean2()
    data_sd2 <- data_sd2()
    data_plot2 <- merge.data.frame(data_mean2,data_sd2, by = c("AGI","Genotype","Time","Temperature"))
  })
  
  
  data3 <- eventReactive(input$AGI3, {
    data3 <- Lai[Lai$AGI == input$AGI3,]
  })
  

  # Draw plots
  #-------------
  Plot1 <- reactive({
    data_plot1 <- data_plot1()
    names(data_plot1)[5:6] <- c("Mean","SD")
    data_plot1$Genotype <- as.factor(data_plot1$Genotype)
    data_plot1$Time <- as.numeric(as.character(data_plot1$Time))
    data_plot1$Temperature <- as.factor(data_plot1$Temperature)
    
    g <- ggplot(data_plot1, aes(x = Time, y = Mean, group= Temperature)) +
      geom_rect(aes(xmin=12, xmax=16, ymin=-Inf, ymax=Inf), 
                fill="grey", alpha=0.16, color = NA)+
      geom_point(data=data_plot1,aes(shape = Temperature, colour = Temperature),
                 position=position_dodge(width=0.1), size = 3)+
      geom_line(data=data_plot1,position=position_dodge(width=0.1),
                aes(linetype = Temperature, colour = Temperature))+
      scale_shape_manual(values=c(18,16))+
      scale_linetype_manual(values=c("blank", "solid"))+
      geom_errorbar(data=data_plot1,aes(ymin=Mean-SD, ymax=Mean+SD, colour = Temperature),width=1.5,
                    position=position_dodge(width=0.1))+
      scale_color_manual(values=c("#0066FF","#666666"))+
      scale_x_continuous(breaks = seq(0,15,3))+
      theme_bw()+
      facet_wrap(~Genotype)+
      xlab("Time of day (ZT)")+
      ylab("Transcript abundance\n(rlog normalized counts)")+
      theme(panel.grid=element_blank(), axis.ticks.length=unit(-0.15, "cm"),
            axis.text=element_text(size=12, colour = "black"),axis.text.x=element_text(margin = unit(c(3, 0, 0, 0), "mm")),
            axis.text.y=element_text(margin=unit(c(0,3,0,0),"mm")),
            strip.text = element_text(size = 14),
            axis.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12))
  })
  
  Plot2 <- reactive({
    data_plot2 <- data_plot2()
    names(data_plot2)[5:6] <- c("Mean","SD")
    data_plot2$Genotype <- factor(data_plot2$Genotype, levels = c("RTX430","Macia"))
    data_plot2$Time <- as.numeric(as.character(data_plot2$Time))
    data_plot2$Temperature <- as.factor(data_plot2$Temperature)
    
    g <- ggplot(data_plot2, aes(x = Time, y = Mean, group= Temperature)) +
      geom_rect(aes(xmin=12, xmax=16, ymin=-Inf, ymax=Inf), 
                fill="grey", alpha=0.16, color = NA)+
      geom_point(data=data_plot2,aes(shape = Temperature, colour = Temperature),
                 position=position_dodge(width=0.1), size = 3)+
      geom_line(data=data_plot2,position=position_dodge(width=0.1),
                aes(linetype = Temperature, colour = Temperature))+
      scale_shape_manual(values=c(16,17))+
      scale_linetype_manual(values=c("solid", "blank"))+
      geom_errorbar(data=data_plot2,aes(ymin=Mean-SD, ymax=Mean+SD, colour = Temperature),width=1.5,
                    position=position_dodge(width=0.1))+
      scale_color_manual(values=c("#666666","#FF6600"))+
      scale_x_continuous(breaks = seq(0,15,3))+
      theme_bw()+
      facet_wrap(~Genotype)+
      xlab("Time of day (ZT)")+
      ylab("Transcript abundance\n(rlog normalized counts)")+
      theme(panel.grid=element_blank(), axis.ticks.length=unit(-0.15, "cm"),
            axis.text=element_text(size=12, colour = "black"),axis.text.x=element_text(margin = unit(c(3, 0, 0, 0), "mm")),
            axis.text.y=element_text(margin=unit(c(0,3,0,0),"mm")),
            strip.text = element_text(size = 14),
            axis.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12))
  })
  
  Plot3 <- reactive({
    data3 <- data3()
    data3$Time <- as.numeric(as.character(data3$Time))

    g <- ggplot(data3, aes(x = Time, y = value)) +
    geom_rect(aes(xmin=12, xmax=24, ymin=-Inf, ymax=Inf),
              fill="grey", alpha=0.08, color = NA)+
    geom_rect(aes(xmin=36, xmax=48, ymin=-Inf, ymax=Inf),
              fill="grey", alpha=0.08, color = NA)+
    geom_rect(aes(xmin=60, xmax=72, ymin=-Inf, ymax=Inf),
              fill="grey", alpha=0.08, color = NA)+
    geom_point(size = 3)+
    geom_line()+
    scale_x_continuous(breaks = seq(0,72,3))+
    theme_bw()+
    xlab("Time of day (ZT)")+
    ylab("Transcript abundance (FPKM)")+
    theme(panel.grid=element_blank(), axis.ticks.length=unit(-0.15, "cm"),
          axis.text=element_text(size=12, colour = "black"),axis.text.x=element_text(margin = unit(c(3, 0, 0, 0), "mm")),
          axis.text.y=element_text(margin=unit(c(0,3,0,0),"mm")),
          strip.text = element_text(size = 14),
          axis.title = element_text(size = 14))
  })
  
  output$One = renderPlot({print(Plot1())
  })
  
  output$Two = renderPlot({print(Plot2())
  })
  
  output$Four = renderPlot({print(Plot3())
  })
  
  output$Three = renderDataTable(
    DT::datatable({
      annot
    },
    extensions = c('Buttons', 'Scroller'),
    options = list(
      scroller = TRUE,
      deferRender = TRUE,
      scrollY = 400
      # dom = 'tB'
    ),
    fillContainer = T,
    class = "display"

  )
  )
  
  # Table and figure legends
  output$Legend_table <-renderText({
    paste("<p style='text-align:justify'>",
          "Data correspond to the <a href='https://data.jgi.doe.gov/refine-download/phytozome?organism=Sbicolor&expanded=454'>Sbicolor_454_v3.1.1.annotation_info</a>.",
          "</p>"
    )
  })
  
  output$Legend_cold <-renderText({
    paste("<p style='text-align:justify'>",
          "Sorghum plants of two different varieties - RTX430 (cold susceptible) and SC224 (cold tolerant) - have been cultivated in photocycles and thermocycles (12 h light at 30\u00B0C/12 h dark at 22\u00B0C) conditions. 
          Temperature stress for 1 h at 10\u00B0C has been applied at different times of day (ZT0-ZT1, ZT5-ZT6, ZT8-ZT9, ZT14-ZT15) on different sets of plants. 
          Control plants were maintained in control temperature conditions. 
          The first two leaves were collected at ZT1, ZT6, ZT9 and ZT15 to isolate mRNAs. Data correspond to mRNA-Seq, with normalized raw counts calculated using the rlog transformation from the <a href='https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8'>'DeSeq2'</a> R package. 
          Data are means +/- SD for n = 3 biological replicates. The grey areas represent the night period.",
          "</p>"
    )
  })
  
  output$Legend_heat <-renderText({
    paste("<p style='text-align:justify'>",
          "Sorghum plants of two different varieties - RTX430 (heat susceptible) and Macia (heat tolerant) - have been cultivated in photocycles and thermocycles (12 h light at 30\u00B0C/12 h dark at 22\u00B0C) conditions. 
          Temperature stress for 1 h at 42\u00B0C has been applied at different times of day (ZT0-ZT1, ZT5-ZT6, ZT8-ZT9, ZT14-ZT15) on different sets of plants. 
          Control plants were maintained in control temperature conditions. 
          The first two leaves were collected at ZT1, ZT6, ZT9 and ZT15 to isolate mRNAs. Data correspond to mRNA-Seq, with normalized raw counts calculated using the rlog transformation from the <a href='https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8'>'DeSeq2'</a> R package. 
          Data are means +/- SD for n = 3 biological replicates. The grey areas represent the night period.",
          "</p>"
    )
  })
  
  output$Legend_diurnal <-renderText({
    paste("<p style='text-align:justify'>",
          "Growth conditions were 12 h light and 12 h dark and daytime temperatures were between 26 and 28\u00B0C and nighttime temperature was 22\u00B0C.
          Fully expanded third leaves at the 3-leaf stage were sampled every 3 h over the course of 3 days (a total of 72 h), and mRNA were isolated. 
          Data correspond to mRNA-Seq, and FPKM values are used to represent transcript abundance. 
          Rhythmic oscillations were determined using the <a href='https://doi.org/10.1177/0748730410379711'>'JTK_Cycle'</a> algorithm. 
          The grey areas represent the night period. For more details, please see <a href='https://doi.org/10.1186/s12864-020-06824-3'>Lai et al. (2020)</a>.",
          "</p>"
    )
  })
  
  # Export cold stress plots
  output$Diurnal_cold_stressPNG <- downloadHandler(
    filename = "Diurnal_cold_stress.png.png",
    content = function(file){
      ggsave(file, plot = Plot1(), width = 9, height = 5, device = "png")
    })
  
  output$Diurnal_cold_stressSVG <- downloadHandler(
    filename = "Diurnal_cold_stress.svg.svg",
    content = function(file){
      ggsave(file, plot = Plot1(), width = 9, height = 5, device = "svg")
    })
  
  output$Diurnal_cold_stressPDF <- downloadHandler(
    filename = "Diurnal_cold_stress.pdf.pdf",
    content = function(file){
      ggsave(file, plot = Plot1(), width = 9, height = 5, device = "pdf")
    })
  
  # Export heat stress plots
  output$Diurnal_heat_stressPNG <- downloadHandler(
    filename = "Diurnal_heat_stress.png.png",
    content = function(file){
      ggsave(file, plot = Plot2(), width = 9, height = 5, device = "png")
    })
  
  output$Diurnal_heat_stressSVG <- downloadHandler(
    filename = "Diurnal_heat_stress.svg.svg",
    content = function(file){
      ggsave(file, plot = Plot2(), width = 9, height = 5, device = "svg")
    })
  
  output$Diurnal_heat_stressPDF <- downloadHandler(
    filename = "Diurnal_heat_stress.pdf.pdf",
    content = function(file){
      ggsave(file, plot = Plot2(), width = 9, height = 5, device = "pdf")
    })
  
  # Export diurnal data
  output$Diurnal_LaiPNG <- downloadHandler(
    filename = "Diurnal_Lai.png.png",
    content = function(file){
      ggsave(file, plot = Plot3(), width = 11, height = 5, device = "png")
    })
  
  output$Diurnal_LaiSVG <- downloadHandler(
    filename = "Diurnal_Lai.svg.svg",
    content = function(file){
      ggsave(file, plot = Plot3(), width = 11, height = 5, device = "svg")
    })
  
  output$Diurnal_LaiPDF <- downloadHandler(
    filename = "Diurnal_Lai.pdf.pdf",
    content = function(file){
      ggsave(file, plot = Plot3(), width = 11, height = 5, device = "pdf")
    })
  
  
  # Reactive value for selected dataset
  datasetInputOne <- reactive({
    rlog2[rlog2$AGI == input$AGI,]
  })
  
  datasetInputTwo <- reactive({
    rlog1[rlog1$AGI == input$AGI2,]
  })
  
  datasetInputThree <- reactive({
    Lai[Lai$AGI == input$AGI3,]
  })
  
  # Downloadable txt of selected dataset
  output$Diurnal_cold_stress <- downloadHandler(
    filename = "Diurnal_cold_stress.txt.txt",
    content = function(file) {
      write.table(datasetInputOne(), file, row.names = FALSE)
    }
  )
  
  output$Diurnal_heat_stress <- downloadHandler(
    filename = "Diurnal_heat_stress.txt.txt",
    content = function(file) {
      write.table(datasetInputTwo(), file, row.names = FALSE)
    }
  )
  
  output$Diurnal_Lai <- downloadHandler(
    filename = "Diurnal_Lai.txt.txt",
    content = function(file) {
      write.table(datasetInputThree(), file, row.names = FALSE)
    }
  )


}

# Run the app ----
shinyApp(ui = ui, server = server)



