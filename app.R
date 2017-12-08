#library(rsconnect)
#setwd("C:\Users\mjohn\Documents\R\shinyapp")
#rsconnect::deployApp('path/to/your/app')
#2.0.0
library(shiny)
library(shinydashboard)
library(gdata)
library('R.matlab')
library('ggplot2')
library('reshape2')
library("R.matlab")
library("rgl")
library("plot3D")
library("plot3Drgl")
library("scatterplot3d")
library("tibble")
library(stringr)
library("plyr")
require(MASS)
library(gridExtra)

ui <- dashboardPage(
  skin = "black",
  dashboardHeader(title = "NasTracks", titleWidth = 350),
  dashboardSidebar(
    width = 350,
    sidebarMenu(
      menuItem("Movement Analysis", icon = icon("area-chart"), tabName = "MA"),
      menuItem("Mating Prediction", icon = icon("area-chart"), tabName = "Mating_Prediction"),
      fileInput("file1", "Upload trajectory File (.txt)",accept = c("text/csv",".txt")),
      actionButton("goButton", "Run Analysis!")
    )
  ),
    
  dashboardBody(
    tags$head(tags$style(HTML("
                               @import url('https://fonts.googleapis.com/css?family=Averia+Serif+Libre');
                              .main-header .logo {
                              font-family: 'Averia Serif Libre', cursive;
                              font-weight: bold;
                              font-size: 32px;
                              font-color: red;
                              }
                              "))),
    #attributes
    fluidRow(
    valueBoxOutput(width = 6, "count_wasps"),
    valueBoxOutput(width = 6, "frames")
    ),
    
    tabItems(
      tabItem("MA",
        fluidRow(
          box(
            width = 6, status = "warning",
            title = "2D Plot",
            plotOutput("plot"),
            numericInput("wasp_id_2dplot", label = h4("Enter wasp number"), value = 1),
            selectInput("color", "Color:",
                        c("Red" = "red",
                          "Green" = "green",
                          "Blue" = "blue",
                          "Black" = "black"))
          ),
          box(
            width = 6, status = "warning",
            title = "Heatmap",
            plotOutput("heatmap"),
            h5("Two-dimensional kernel density estimation"),
            numericInput("wasp_id_heatmap", label = h4("Enter wasp number"), value = 1),
            sliderInput("bandwidths", "Bandwidths:",
                        min = 1, max = 200, value = 75
            ),
            sliderInput("gridpoints", "Grid points:",
                        min = 1, max = 200, value = 50
            ),
            numericInput("x_pixels", label = h5("Resolution (x axis)"), value = 1800),
            numericInput("y_pixels", label = h5("Resolution (y axis)"), value = 1200)
          )
        )
      ),

      tabItem("Mating_Prediction",
        fluidRow(
          tabBox(
            title = "Analysis",
            tabPanel("Neighbor attraction time-series", plotOutput("mate_prediction")),
            tabPanel("Histogram", plotOutput("hist_plots"),
                     h4("Peak data"),
                     tableOutput("peak_attr"))
          ),
          
          box(
            width = 6, status = "warning",
            title = "Parameters",
            "Select two wasps", br(),
            numericInput("wasp1_id_social", label = h5("Wasp ID"), value = 1),
            numericInput("wasp2_id_social", label = h5("Wasp ID"), value = 2),
            hr(),
            "Smoothing algorithm: Nadaraya-Watson kernel regression", br(),
            sliderInput("neck", label = h5("peak threshold"), min = 1, 
                        max = 1000, value = 500),
            checkboxInput("smooth", label = "Enable smoothing", value = TRUE),
            sliderInput("bw", label = h5("bandwidth"), min = 1, 
                        max = 1000, value = 20),
            checkboxInput("filter", label = "floor threshold (beta)", value = FALSE),
            sliderInput("thresh", label = h5("Threshold"), min = 1, 
                        max = 1000, value = 200),
            checkboxInput("flip", label = "flip axis", value = TRUE),
            h4("Peak finding"),
            h5("Higher m values are stricter"),
            sliderInput("m", label = h5("m value"), min = 1, 
                        max = 5000, value = 1000)
          )
        )

       )
      )
    )
  )
        
#########################################################################################

server <- function(input, output) {
  
  output$debug <- renderTable({
    return(NULL)
  })
  
  idTracker <- eventReactive(input$goButton, {
    inFile <- input$file1
    if (is.null(inFile))
      return(NULL)
    read.delim(inFile$datapath, header = TRUE, strip.white = TRUE, sep = "")})
  
  output$count_wasps <- renderValueBox({
    count<-dim(idTracker())[2]/3
    
    valueBox(
      value = count,
      subtitle = "Total number of wasps",
      icon = icon("list-ol"),
      color = "olive" 
    )
  })
  
  
  count_wasps <- reactive({
    count<-dim(idTracker())[2]/3
    return(count)})
  
  output$frames <- renderValueBox({
    frames <- dim(idTracker())[1]
    
    valueBox(
      value = frames,
      subtitle = "Total number of frames",
      icon = icon("sort-numeric-asc"),
      color = "olive"
    )
  })
  
  frames <- renderText({
    frames <- dim(idTracker())[1]
    return(frames)})
  
  observeEvent(input$goButton, {
    frames <- as.numeric(frames())
    if (input$goButton == 1){
      insertUI("#goButton", "afterEnd",
               numericInput("FOI_e", "End frame", value = frames)
      )
      insertUI("#goButton", "afterEnd",
               numericInput("FOI_s", "Start frame", value = 1)
      )
      insertUI("#goButton", "afterEnd",
               ui = "Choose frames of interest")
    }})
  
  time_series <- reactive({time_series <- (1:frames()/input$fps)/60})
  
  ff <- reactive({ff <- input$FOI_s})
  lf <- reactive({lf <- input$FOI_e})
  
  output$plot <- renderPlot({
    
    a <- as.double(input$wasp_id_2dplot)
    x <- a+(2*(a-1)) #thanks kim!
    y <- x+1
    ff <- as.double(ff())
    lf <- as.double(lf())
    x <- idTracker()[ff:lf,x]
    y <- idTracker()[ff:lf,y]
    p<-data.frame(x,y)
    # plot <- ggplot(idTracker()[ff:lf,], aes(x,y)) +
    #          geom_point(alpha = input$alpha, color = input$color) + coord_fixed(ratio = 1)
    plot <- plot(x,y, type="l", col=input$color, asp=1, ylim = rev(range(y)))
    
    return(plot)})
  
  output$heatmap <- renderPlot({
    a <- as.double(input$wasp_id_heatmap)
    xp <- as.double(input$x_pixels)
    yp <- as.double(input$y_pixels)
    x <- a+(2*(a-1))
    y <- x+1
    ff <- as.double(ff())
    lf <- as.double(lf())
    size <- lf-ff
    x <- idTracker()[1:size,x]
    y <- idTracker()[1:size,y]
    x <- x[!is.na(x)]
    y <- y[!is.na(y)]
    map <- list(x,y)
    d <- structure(map, .Names = c("X", "Y"), class = "data.frame", row.names = c(NA, -size-1))
    dens <- kde2d(d$X, d$Y, h=input$bandwidths, n=input$gridpoints, lims = c(1, xp, 1, yp))
    filled.contour(dens, color=terrain.colors, asp=1, ylim = rev(range(y)) )})
  
  output$distance_by_time <- renderPlot({
    ff <- as.double(ff())
    lf <- as.double(lf())
    low <- as.numeric(input$prox_range[1])
    high <- as.numeric(input$prox_range[2])
    size <- dim(idTracker())[1]
    w1 <- as.double(input$wasp1_id_social)
    x1 <- w1+(2*(w1-1))
    y1 <- x1+1
    w2 <- as.double(input$wasp2_id_social)
    x2 <- w2+(2*(w2-1))
    y2 <- x2+1
    x1 <- idTracker()[2:size,x1]
    y1 <- idTracker()[2:size,y1]     
    x2 <- idTracker()[2:size,x2]
    y2 <- idTracker()[2:size,y2]
    compX <- abs(x1-x2)
    compY <- abs(y1-y2)
    x <- data.frame(compX)
    y <- data.frame(compY)
    hypot <- sqrt(x^2+y^2)
    # plot <- ggplot(data=idTracker()[2:size,]) +
    #          geom_line(mapping=aes(x=2:size,y=hypot), color = "black") + 
    #          xlim(ff,lf)+
    #          geom_segment(mapping=aes(x=ff,y=low,xend=lf,yend=low),color="red", linetype = 2)+
    #          geom_segment(mapping=aes(x=ff,y=high,xend=lf,yend=high),color="red", linetype = 2)
    timeline <- c(ff,lf)
    plot(x=1:length(hypot), y=hypot, type="l", xlim=timeline)
    abline(h=low, col="red", lty=2)
    abline(h=high, col="red", lty=2)
    xlab("Frames")
    ylab("proximity (px)")
    return(plot)})
  
  
  output$mate_prediction <- renderPlot({
    #peak_attr <- peak_attr()
    #peak_attr <- as_tibble(peak_attr)
    #hist(peak_attr$time.of.peak)
    idTracker <- idTracker()
    idTracker<-idTracker[complete.cases(idTracker), ]
    ff <- as.double(ff())
    lf <- as.double(lf())
    timeline <- c(ff,lf)
    low <- as.numeric(input$prox_range[1])
    high <- as.numeric(input$prox_range[2])
    size <- dim(idTracker)[1]
    w1 <- as.double(input$wasp1_id_social)
    x1 <- w1+(2*(w1-1))
    y1 <- x1+1
    w2 <- as.double(input$wasp2_id_social)
    x2 <- w2+(2*(w2-1))
    y2 <- x2+1
    x1 <- idTracker[2:size,x1]
    y1 <- idTracker[2:size,y1]     
    x2 <- idTracker[2:size,x2]
    y2 <- idTracker[2:size,y2]
    compX <- abs(x1-x2)
    compY <- abs(y1-y2)
    x <- data.frame(compX)
    y <- data.frame(compY)
    hypot <- sqrt(x^2+y^2)
    y <- hypot
    x <- c(1:length(y))
    y_out <- smooth.me()(x,y, smooth = input$smooth, input$bw , filter = input$filter, input$thresh, flip = input$flip)
    peaks <- find.peaks()(y_out[,2],input$m)
    plot(y_out, type="l", xlim=timeline)
    for (i in 1:length(peaks)) {points(peaks[i], y_out[peaks[i],2], col="blue", pch=16)}
    
    #indexing
    forward <- list()
    reverse <- list()
    peak_location <- list()
    peak_height <- list()
    neck <- input$neck
    abline(h=-input$neck, col = "red", lty=2)
    frames <- dim(y_out)[1]
    for (p in peaks){
      if (neck + y_out[p,2] < 0) next #skips peak if its below neck threshold
      f <- 0
      r <- 0
      b <- 0
      a <- 0
      name <- paste0("peak_", p, sep="")
      peak_location[[name]] <- p
      peak_height[[name]] <- y_out[p,2]
      while(f < neck) 
      {
        a <- a + 1
        if (p+a==frames) {f <- NA; break}
        f <- 0-y_out[p+a,2]
      }
      points(p+a,y_out[p+a,2], col="red", pch=16)
      name <- paste0("peak_", p, sep="")
      forward[[name]] <- a
      
      while(r < neck)
      { name <- paste0("peak_", p, sep="")
      b <- b + 1
      if (p-b == 1) {r <- NA; break}
      r <- 0-y_out[p-b,2]
      }
      points(p-b,y_out[p-b,2], col="red", pch=16)
      reverse[[name]] <- b
    }
    
    
    
  })
  
  peak_attr <- reactive({
    idTracker <- idTracker()
    idTracker<-idTracker[complete.cases(idTracker), ]
    ff <- as.double(ff())
    lf <- as.double(lf())
    low <- as.numeric(input$prox_range[1])
    high <- as.numeric(input$prox_range[2])
    size <- dim(idTracker)[1]
    w1 <- as.double(input$wasp1_id_social)
    x1 <- w1+(2*(w1-1))
    y1 <- x1+1
    w2 <- as.double(input$wasp2_id_social)
    x2 <- w2+(2*(w2-1))
    y2 <- x2+1
    x1 <- idTracker[2:size,x1]
    y1 <- idTracker[2:size,y1]     
    x2 <- idTracker[2:size,x2]
    y2 <- idTracker[2:size,y2]
    compX <- abs(x1-x2)
    compY <- abs(y1-y2)
    x <- data.frame(compX)
    y <- data.frame(compY)
    hypot <- sqrt(x^2+y^2)
    y <- hypot
    x <- c(1:length(y))
    y_out <- smooth.me()(x,y, smooth = input$smooth, input$bw , filter = input$filter, input$thresh, flip = input$flip)
    peaks <- find.peaks()(y_out[,2],input$m)
    #plot(y_out, type="l")
    #for (i in 1:length(peaks)) {points(peaks[i], y_out[peaks[i],2], col="blue", pch=16)}
    
    #indexing
    forward <- list()
    reverse <- list()
    peak_location <- list()
    peak_height <- list()
    neck <- input$neck
    #abline(h=-input$neck, col = "red", lty=2)
    frames <- dim(y_out)[1]
    for (p in peaks){
      if (neck + y_out[p,2] < 0) next #skips peak if its below neck threshold
      f <- 0
      r <- 0
      b <- 0
      a <- 0
      name <- paste0("peak_", p, sep="")
      peak_location[[name]] <- p
      peak_height[[name]] <- y_out[p,2]
      while(f < neck) 
      {
        a <- a + 1
        if (p+a==frames) {f <- NA; break}
        f <- 0-y_out[p+a,2]
      }
      #points(p+a,y_out[p+a,2], col="red", pch=16)
      name <- paste0("peak_", p, sep="")
      forward[[name]] <- a
      
      while(r < neck)
      { name <- paste0("peak_", p, sep="")
      b <- b + 1
      if (p-b == 1) {r <- NA; break}
      r <- 0-y_out[p-b,2]
      }
      #points(p-b,y_out[p-b,2], col="red", pch=16)
      reverse[[name]] <- b
    }
    #Data frame building
    neck_duration <- lapply(seq_along(forward), function(i)
      unlist(forward[i])+unlist(reverse[i]))
    
    peak_location <- unlist(peak_location)
    neck_duration <- unlist(neck_duration)
    start_time <- unlist(forward)
    end_time <- unlist(reverse)
    peak_height <- unlist(peak_height)
    abs_start_time <- peak_location - start_time
    abs_end_time <- peak_location + end_time
    
    peak_attr <- data.frame(
      "time of peak" = peak_location,
      "duration of peak" = neck_duration, 
      "start time" = start_time,
      "end time" = end_time,
      "absolute start time" = abs_start_time,
      "absolute end time" = abs_end_time,
      "peak height" = peak_height)
    return(peak_attr)
    
  })
  
  output$peak_attr <- renderTable({
    peak_attr()
  })
  
  output$hist_plots <- renderPlot({
    par(mfrow=c(2,2))
    peak_attr <- peak_attr()
    hist(peak_attr$duration.of.peak, main ="Histogram of mate attempt duration")
    hist(peak_attr$time.of.peak, main="Histogram of mate attempt time")
    hist(peak_attr$peak.height, main="Histogram of max closeness")
    hist(peak_attr$absolute.start.time, main="Histogram of mate attempt initiation time")
  })
  
  
  output$prox_calc <- renderText({
    ff <- as.double(ff())
    lf <- as.double(lf())
    size <- dim(idTracker())[1]
    w1 <- as.double(input$wasp1_id_social)
    x1 <- w1+(2*(w1-1))
    y1 <- x1+1
    w2 <- as.double(input$wasp2_id_social)
    x2 <- w2+(2*(w2-1))
    y2 <- x2+1
    x1 <- idTracker()[ff:lf,x1]
    y1 <- idTracker()[ff:lf,y1]     
    x2 <- idTracker()[ff:lf,x2]
    y2 <- idTracker()[ff:lf,y2]
    both <- data_frame(x1,y1)
    compX <- abs(x1-x2)
    compY <- abs(y1-y2)
    x <- data.frame(compX)
    y <- data.frame(compY)
    hypot <- sqrt(x^2+y^2)
    low <- as.numeric(input$prox_range[1])
    high <- as.numeric(input$prox_range[2])
    values <- sum(hypot>low & hypot<high, na.rm = TRUE)
    return(values)})
  
  output$probs <- renderPlot({
    num <- seq(3, ncol(idTracker()), by=3)
    prob <- idTracker()[num]
    return(NULL)})
  
  output$version_log <- renderText({
    log<-c("
           0.9.4
           8/14/2017
           * heatmap bug fixed
           * ability to input video resolution
           * change the frames of interest with own input values
           ")
    return(log)})
  
  smooth.me <- reactive({
    smooth.me <- function(x, y, smooth, bw, filter, thresh, flip){
      if (!smooth){
        #potential_hits <- ifelse(y<thresh,y,0)
        return(y)
      }
      else{
        if (!flip){
          if (!filter){
            y <- ksmooth(x,y,"normal", bandwidth = bw)
            y <- data.frame(y)
            return(y)
          }
          else{
            y <- lapply(y, function(x) if(x>thresh){x<-thresh}else{x<-x})
            y <- ksmooth(x,y,"normal", bandwidth = bw)
            y <- data.frame(y)
            return(y)
          }
        }
        else{#flipping
          max_value <- max(y)
          min_value <- min(y)
          if (!filter){
            y_flip <- lapply(y, function(x) min_value-x)
            y <- ksmooth(x,y_flip,"normal", bandwidth = bw)
            y <- data.frame(y)
            return(y)
          }
          else{
            y <- lapply(y, function(x) if(x>thresh){x<-thresh}else{x<-x})
            y_flip <- lapply(y, function(x) min_value-x)
            y <- ksmooth(x,y_flip,"normal", bandwidth = bw)
            y <- data.frame(y)
            return(y)
          }
        }
      }
    }
  })
  
  find.peaks <- reactive({
    find.peaks <- function (x, m = 3){
      shape <- diff(sign(diff(x, na.pad = FALSE)))
      pks <- sapply(which(shape < 0), FUN = function(i){
        z <- i - m + 1
        z <- ifelse(z > 0, z, 1)
        w <- i + m + 1
        w <- ifelse(w < length(x), w, length(x))
        if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
      })
      pks <- unlist(pks)
      pks
    }
  })
}

shinyApp(ui=ui, server=server)
