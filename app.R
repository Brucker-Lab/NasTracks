#library(rsconnect)
#setwd("C:/Users/mjohn/Documents/R/shinyapp")
#rsconnect::deployApp()
#2.0.0
library(shiny)
library(shinydashboard)
library(ggplot2)
library(reshape2)
library(stringr)
library(MASS)
library(tibble)
library(RCurl)
library(dplyr)
library(dbscan)

ui <- dashboardPage(
  skin = "black",
  dashboardHeader(title = "NasTracks", titleWidth = 350),
  dashboardSidebar(
    width = 350,
    sidebarMenu(
      menuItem("Movement Analysis", icon = icon("area-chart"), tabName = "MA"),
      menuItem("Mating Prediction", icon = icon("area-chart"), tabName = "Mating_Prediction"),
      menuItem("Chasing Prediction (beta)", icon = icon("area-chart"), tabName = "Chase_Prediction"),
      fileInput("file1", "Upload trajectory File (.txt)",accept = c("text/csv",".txt")),
      #actionButton("demo", "Use demo data"),
      checkboxInput("demo", "Use demo data", FALSE),
      actionButton("goButton", "Run analysis!")
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
                  tableOutput("debug"),
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
                              max = 1000, value = 100),
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
      ),
      
      tabItem("Chase_Prediction",
              fluidRow(
                tabBox(
                  title = "Analysis",
                  tabPanel("Chasing events", plotOutput("chase_prediction"))
                ),
                
                box(
                  width = 6, status = "warning",
                  title = "Parameters",
                  "Select two wasps", br(),
                  numericInput("SPD", label = h5("Speed"), value = 2),
                  numericInput("BOX", label = h5("Proximity"), value = 75),
                  numericInput("CONE", label = h5("Theta"), value = 0.35)
                )
              )
      )
    )
    )
    )


#########################################################################################

server <- function(input, output) {
  
  idTracker <- eventReactive(input$goButton, {
    if (input$demo){
      x <- getURL('https://gist.githubusercontent.com/mjoh223/7d66d0108a03e6fe65de32f213ba86a0/raw/f382d93e5a1b42e8c9f5e2d6c6872ca59502b7b7/gistfile1.txt')
      demo <- read.delim(text = x, header = TRUE, strip.white = TRUE, sep = "")
      return(demo)
    } else{
      inFile <- input$file1
      if (is.null(inFile)) {
        return(NULL)} else{
          read.delim(inFile$datapath, header = TRUE, strip.white = TRUE, sep = "")}
    }})
  
  output$debug <- renderTable({
    return(NULL)
  })
  
  idTrackerFOI <- eventReactive({c(input$ID,ff(),lf())}, {
    a <- as.double(input$ID)
    x <- a+(2*(a-1))
    ff <- as.double(ff())
    lf <- as.double(lf())
    wx <- idTracker()[ff():lf(),x]
    wy <- idTracker()[ff():lf(),x+1]
    p<-data.frame(wx,wy)
    return(p)
  })
  
  output$plot <- renderPlot({
    plot <- plot(idTrackerFOI(), type="l", col=input$color, asp=1, xlab = "x (px)", ylab = "y (px)", main = "Animal trajectory")
    return(plot)
  })
  
  output$heatmap <- renderPlot({
    xp <- as.double(input$x_pixels)
    yp <- as.double(input$y_pixels)
    map <- na.omit(idTrackerFOI())
    d <- structure(map, .Names = c("X", "Y"), class = "data.frame", row.names = c(NA, -nrow(map)-1))
    dens <- kde2d(d$X, d$Y, h=input$bandwidths, n=input$gridpoints, lims = c(1, xp, 1, yp))
    filled.contour(dens, color=terrain.colors, asp=1)
  })
  
  output$chase_prediction <- renderPlot({
    NUM_WASPS <- dim(idTracker())[2]/3
    SPD <- as.numeric(input$SPD)
    BOX <- as.numeric(input$BOX)
    CONE <- as.numeric(input$CONE)
    dat <- na.omit(idTracker())
    dat_rad <- list()
    dat_spd <- list()
    cat <- list()
    n <- 1:NUM_WASPS
    x<-n+(2*(n-1))
    y<-x+1
    dy <- as.matrix(dat[-1,y]-dat[-nrow(dat),y])
    dx <- as.matrix(dat[-1,x]-dat[-nrow(dat),x])
    dat_rad <- as_tibble(atan2(dy, dx))
    dat_spd <- as_tibble(sqrt(abs(dx^2)+abs(dy^2)))
    # #Ensures parameters are met
    chase <- vector('list', length(dat))
    combn_index <- combn(1:NUM_WASPS,2)
    chase <- apply(combn_index, 2, function(x) {
      x1<-x[1]+(2*(x[1]-1))
      y1<-x1+1
      x2<-x[2]+(2*(x[2]-1))
      y2<-x2+1
      p1<-x[1]*3
      p2<-x[2]*3
      mutate(dat[-1,], idx = 2:nrow(dat)) %>%
        dplyr::filter(
          abs(dat_rad[x[1]] - dat_rad[x[2]]) < CONE &
            dat_spd[x[1]] > SPD &
            dat_spd[x[2]] > SPD &
            abs(as.matrix(dat[-1,x1]) - as.matrix(dat[-1,x2])) < BOX &
            abs(as.matrix(dat[-1,y1]) - as.matrix(dat[-1,y2])) < BOX) %>%
        dplyr::select(idx,x1,y1,p1,x2,y2,p2)
    })
    #check chase dimensions and replace no chases with NA
    for (c in seq_along(chase)){
      #print(paste0(c, " ", cc, " ", nrow(chase[[c]]), sep=""))
      if (nrow(chase[[c]])==0){
        chase[[c]] <- NA
      }
    }
    #clusters
    cluster <- list()
    for (k in seq_along(chase)){
      idx <- chase[[k]][1]
      chase.y <- dat_spd[as.matrix(idx), k]
      if (length(chase[[k]]) > 1 ){
        t <- data.frame(idx, chase.y)
        cluster[[k]] <-dbscan(as.matrix(t), eps = 25, minPts = 15, borderPoints = T)["cluster"]
      } else {
        cluster[[k]] <- NA
      }
    }
    #master data frame with location, cluster, and prob
    master <- list()
    for (k in seq_along(chase)){
      master[[k]] <- cbind(cluster[[k]], chase[[k]])
    }
    #remove 0 grp
    master_v2 <- list()
    for (k in seq_along(master)){
      if (length(chase[[k]])!=1){
        master_v2[[k]] <- filter(master[[k]], cluster!=0)
      } else {
        master_v2[[k]] <- list(NA)
      }
    }
    #check master v2 dimensions
    for (c in seq_along(master_v2)){
      if (length(unlist(master_v2[[c]])) < 1){
        master_v2[[c]] <- list(NA)
      }
      #print(paste0(c, " ", cc, " ", length(unlist(master_v2[[c]][[cc]])), sep=""))
    }
    #graph for verification
    i <- as.numeric(input$combn)
    combn_index <- combn(1:NUM_WASPS,2)
    idv1 <- combn_index[1,i]
    idv2 <- combn_index[2,i]
    p <- ggplot() +
      geom_point(data = master_v2[[i]], mapping = aes(x=master_v2[[i]][,3], y=master_v2[[i]][,4], color = paste0("wasp ", idv1)), color="orange") +
      geom_point(data = master_v2[[i]], mapping = aes(x=master_v2[[i]][,6], y=master_v2[[i]][,7], color = paste0("wasp ", idv2)), color="darkolivegreen") +
      facet_wrap(~cluster, scales = "free") +
      xlab("x (px)") +
      ylab("y (px)") +
      theme(legend.title=element_blank())
    return(p)
  })
  
  output$count_wasps <- renderValueBox({
    count<-dim(idTracker())[2]/3
    valueBox(
      value = count,
      subtitle = "Total number of wasps",
      icon = icon("list-ol"),
      color = "olive" 
    )
  })
  
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
      insertUI("#goButton", "afterEnd",
               numericInput("ID", "animal ID", value = 1)
      )
      insertUI("#goButton", "afterEnd",
               ui = "Choose ID of animal of interest")
      insertUI("#CONE", "afterEnd",
               #numericInput("combn", "pairwise index", value = 1),
               selectInput("combn", "Choose comparision:",
                           c("Wasp1:Wasp2" = 1,
                             "Wasp1:Wasp3" = 2,
                             "Wasp2:Wasp3" = 3))
      )
    }})
  
  time_series <- reactive({time_series <- (1:frames()/input$fps)/60})
  
  ff <- reactive({ff <- input$FOI_s})
  lf <- reactive({lf <- input$FOI_e})
  
  output$mate_prediction <- renderPlot({
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
    plot(y_out, type="l", xlim=timeline, xlab = "time (frames)", ylab = "neighbor proximity (px)")
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
  
  output$version_log <- renderText({
    log<-c("
           1.0.0
           1/27/2018
           * plotting bug fixed
           * quicker analysis
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

