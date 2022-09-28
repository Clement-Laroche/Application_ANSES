############### packages #############
library(shiny)
library(udunits2)
library(dplyr)
library(htmltools)
library(stringr)
library(sp)
library(sf)
library(ggmap)
library(zoo)
library(Rcpp)
library(lubridate)
library(leaflet)
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(leafpop)
library(leaflet.minicharts)
library(ClustGeo)
library(randomcoloR)
library(rlist)
library(transport)
library(rPref)
library(viridis)
############# pretraitement necessaire ##########
load("Appli_prosulfocarb.Rdata")
input_max$LOQ <- NA
for(i in 1:nrow(input_max))
{
  input_max$LOQ[i] <- max(tab_res$LOQ[which((tab_res$d == input_max$d[i])&(tab_res$C == input_max$C[i]))])
}
input_max$LOQ[which(is.na(input_max$LOQ))] <- 0.005
# load("Prosulfo_appli.Rdata")
sta_crs_o <- clust_results
sta_crs_o$comp <- as.factor(sta_crs_o$comp)
data <- tab_res
map_eau <- select( map_eau_new, -c(SEG_DEB,SEG_FIN))
xmin <- min(st_coordinates(map_reg)[,1])
xmax <- max(st_coordinates(map_reg)[,1])
ymin <- min(st_coordinates(map_reg)[,2])
ymax <- max(st_coordinates(map_reg)[,2])
sta_crs_o <- st_transform(sta_crs_o,4326)
g <- ggmap(ggmap = get_map(c(xmin,ymin,xmax,ymax)),extent = "device")
g <- g+geom_sf(data = map_reg_her,inherit.aes = FALSE,alpha = 0)
# g <- g+geom_line(data=essai,
#                  aes(x=lon, y=lat, group=V5),
#                  col = "blue",alpha = 0.5,lwd = 0.5)+
#   xlab("Longitude")+ylab("Latitude")

##### Obtention des clusterings #####

# pos <- which(unlist(lapply(lapply(KLUST_MAT,"labels"),"length"))>1)
# SCORE_FUN <- rep(NA,length(KK))
# for(i in 1:length(KK))
# {
#   SCORE_FUN[i] <- sum(mapply(FUN = "withindiss",KLUST_MAT[pos],KK[[i]]))
# }
# elb <- cbind(SCORE_FUN,1:length(SCORE_FUN)+8)
# elb <- as.data.frame(elb)
# L_1 <- c()
# L_2 <- c()
# big_elb <- c()
# for(l1 in 2:(length(SCORE_FUN)-2))
# {
#   for(l2 in (l1+1):(length(SCORE_FUN)-1))
#   {
#     L_1 <- c(L_1,l1)
#     L_2 <- c(L_2,l2)
#     ml1 <- lm(formula = log(elb$SCORE_FUN)[1:l1]~elb$V2[1:l1])
#     ml2 <- lm(formula = log(elb$SCORE_FUN)[l1:l2]~elb$V2[l1:l2])
#     ml3 <- lm(formula = log(elb$SCORE_FUN)[l2:length(SCORE_FUN)]~elb$V2[l2:length(SCORE_FUN)])
#     big_elb <- c(big_elb,sum(ml1$residuals^2)+sum(ml2$residuals^2)+sum(ml3$residuals^2))
#   }
# }
# l1 <- L_1[which.min(big_elb)]
# l2 <- L_2[which.min(big_elb)]
# ml1 <- lm(formula = log(elb$SCORE_FUN)[1:l1]~elb$V2[1:l1])
# ml2 <- lm(formula = log(elb$SCORE_FUN)[l1:l2]~elb$V2[l1:l2])
# ml3 <- lm(formula = log(elb$SCORE_FUN)[l2:length(elb$SCORE_FUN)]~elb$V2[l2:length(elb$SCORE_FUN)])
pos <- which(substr(colnames(sta_crs_o),start = 1,stop = 5) ==  "clust")
l1 <- substr(colnames(sta_crs_o)[min(pos)],start = nchar(colnames(sta_crs_o)[min(pos)])-1,stop = nchar(colnames(sta_crs_o)[min(pos)]))
l1 <- as.numeric(l1)
l2 <- substr(colnames(sta_crs_o)[max(pos)],start = nchar(colnames(sta_crs_o)[max(pos)])-1,stop = nchar(colnames(sta_crs_o)[max(pos)]))
l2 <- as.numeric(l2)
##### direction des rivières #####
# récupération d'un segment par cours d'eau

flow_dir <- map_eau[order(map_eau$ID_C_EAU),]
pos <- which(!is.na(flow_dir$ID_C_EAU))
flow_dir <- flow_dir[pos,]
pos <- summary(as.factor(flow_dir$ID_C_EAU),maxsum = length(unique(flow_dir$ID_C_EAU)))
pos <- as.numeric(pos)
pos <- cumsum(pos)
# pos1 : si l'on veut plus de segments par cours d'eau
# pos1 <- sort(c(pos-round(diff(c(1,pos))/2),pos))
flow_dir <- flow_dir[pos,]
flow_dir <- as.data.table(st_coordinates(x = flow_dir$geometry))
direction <- flow_dir %>% group_by(flow_dir$L1) %>% summarize(N = n())
direction <- cumsum(direction$N)
# on trace la flèche de direction entre le premier
# point du segment et le dernier point du segment
flow_dir <- flow_dir[sort(c(c(1,direction[-length(direction)]+1),direction)),]
#flow_dir <- flow_dir[sort(c(direction-1,direction)),]
flow_dir <- st_as_sf(x = flow_dir,coords = c(1,2))
st_crs(flow_dir) <- 4326
flow_dir <- st_transform(x = flow_dir,crs = st_crs(map_reg_her))


##### initialisation de la carte #####

# pos <- which(coord_sta$id_station %in% data$sta[data$d %in% seq(input_max$d[1],input_max$d[crops_weibull$segments[[1]][1]-1],by = "day")])
# inf_map <- data.frame(matrix(NA,ncol = 4, nrow = length(pos)))
# inf_map[,1] <- coord_sta$id_station[pos]
# tb <- data[which(data$d %in% seq(input_max$d[1],input_max$d[crops_weibull$segments[[1]][1]-1],by = "day")),]
# n <- tb %>% group_by(sta) %>% summarise(n = n(),n_quant = sum(1-col))
# inf_map[,2] <- n$n[match(inf_map$X1,n$sta)]
# inf_map[,4] <- n$n_quant[match(inf_map$X1,n$sta)]

sourceCpp("Weibull2.cpp")
# for(i in 1:nrow(inf_map))
# {
#     inf_map[i,3] <- costfcpp(data = tb[tb$sta == inf_map[i,1],],tstart = 0,tstop = inf_map[i,2]-1,init = "mean",q_init = 0.5,k = 1/2,prec = 10^-8,Nmax = 100,lmax = 10)[1]
# }
# l_init <- costfcpp(data = input_max[,2:3],tstart = 0,tstop = crops_weibull$segments[[1]][1]-1,init = "mean",q_init = 0.5,k = 0.5,prec = 10^-8,Nmax = 100,lmax = 10)[2]
# Quantification <- as.factor(1-input_max$col)
# crops_exp$beta<-round(crops_exp$beta,2)
# crops_weibull$beta<-round(crops_weibull$beta,2)
# crops_frechet$beta<-round(crops_frechet$beta,2)

K_l50 <- unlist(lapply(crops_weibull$segments,FUN = "length"))
# slope <- rep(NA,length(K_l50)-2)
# for(i in 2:(length(K_l50)-1))
# {
#   ml1 <- lm(crops_weibull$cost[1:i]~K_l50[1:i])
#   ml2 <- lm(crops_weibull$cost[(i):length(K_l50)]~K_l50[(i):length(K_l50)])
#   slope[i-1] <- sum((ml1$residuals)^2)+sum((ml2$residuals)^2)
# }
# K_star <- (2:length(K_l50))[which.min(slope)]
ml1 <- lm(crops_weibull$cost[1:K_star]~K_l50[1:K_star])
ml2 <- lm(crops_weibull$cost[K_star:length(K_l50)]~K_l50[K_star:length(K_l50)])

 m <- leaflet() %>%
   setView(lng = (max(st_coordinates(sta_crs_o)[,1])+min(st_coordinates(sta_crs_o)[,1]))/2,
           lat = (max(st_coordinates(sta_crs_o)[,2])+min(st_coordinates(sta_crs_o)[,2]))/2,
           zoom = 5) %>%
   addTiles() %>%
   addPolygons(data = map_reg_her,
               color = "black",
               fillOpacity = 0,
               opacity = 1,
               stroke = 0.5,
               weight = 3) %>%
   addPolylines(data = st_zm(x = map_eau$geometry,drop = T,what = "ZM"),
                weight = 2,color = "darkblue"
   ) %>%
   addFlows(lng0 = as.numeric(st_coordinates(flow_dir)[1:nrow(flow_dir)%%2 == 0,1]),
            lat0 = as.numeric(st_coordinates(flow_dir)[1:nrow(flow_dir)%%2 == 0,2]),
            lng1 = as.numeric(st_coordinates(flow_dir)[1:nrow(flow_dir)%%2 == 1,1]),
            lat1 = as.numeric(st_coordinates(flow_dir)[1:nrow(flow_dir)%%2 == 1,2]),
            dir = -1,minThickness = 1.2,maxThickness = 1.2,color = "darkblue",opacity = 1)


#################### UI #######################

# Define UI for application that draws a histogram
ui <- navbarPage("Surveillance de concentrations",
    tabPanel(title = "Accueil",
             column(width = 12,
                    fluidRow(style = "border: 1px solid gray;",
                             h5(strong("Informations générales :")),
                             verbatimTextOutput("Desc_mod1")),
                    fluidRow(style = "border: 1px solid gray;",
                             h5(strong("Tracé des maximum journaliers :")),
                             plotOutput(outputId = "Plot_Max",
                                        height = 250)),
                    fluidRow(style = "border: 1px solid gray;",
                             h5(strong("Carte des stations :")),
                             plotOutput(outputId = "Plot_comp"))
                    )
            ),
    tabPanel(title = "Détection",
             column(width = 12,
               fluidRow(style = "border: 1px solid gray;",
                        h4(strong("1. Détection temporelle : ")),
               column(width = 3,
                      # fluidRow(uiOutput("sliders")),
                      sliderInput(inputId = "pen",label =  "Valeur de pénalité : ", min = round(min(crops_weibull$beta)+10^-2,2),
                                  max = round(max(crops_weibull$beta)+10^-2,2), value = round(crops_weibull$beta[K_star]+10^-2,2),round = 2,step = 0.05
                      ),
                      fluidRow(plotOutput(outputId = "Plot_part",height = 180)),
                      fluidRow(checkboxInput(inputId = "button",label =  "Echelle Log",value = FALSE))
                      ),
               column(width = 9,
                      fluidRow(plotOutput("Plot_Seg", click = "plot_click", height=300))
             )),
             fluidRow(style = "border: 1px solid gray;",
                      column(width = 12,
                             h4(strong("2. Informations sur le segment selectionné")),
                             column(width = 4,verbatimTextOutput("Desc1_mod1")),
                             column(width = 4,plotOutput("Plot_Dens", height=250)),
                             column(width = 4,plotOutput("Sais_mod1", height=250))
                     )),
             fluidRow(style = "border: 1px solid gray;",
                      column(width = 12,
                             h4(strong("3. Détection spatiale")),
                             column(width = 4,
                                    fluidRow(sliderInput("Kluster", "Nombre de cluster : ", min = l1,max = l2, value = l1,step = 1),
                                    selectInput(inputId = "inf",label =  "Information carte",
                                                choices = c("Stations actives","Composantes hydrographiques","Clusters spatiaux","Valeurs du front de Pareto"),
                                                selected = "Stations actives"))
                      ),
                      column(width = 8,fluidRow(leafletOutput("mymap", height=300)))
                      )),
             fluidRow(style = "border: 1px solid gray;",
                      column(width = 12,
                             h4(strong("4. Informations complémentaires")),
                             column(width = 4,
                                    fluidRow(plotOutput(outputId = "Plot_sta",click = "mymap_marker_click",height = 250))
                                    ),
                             column(width = 4,
                                    fluidRow(plotOutput(outputId = "Plot_pareto",click = "mymap_marker_click",height = 250))
                             ),
                             column(width = 4,
                                    fluidRow(verbatimTextOutput("Klust_info"))
                             )
                             )
                      # ,
                      #                 column(width = 6,
                      #                        plotOutput(outputId = "Plot_sta",click = "mymap_marker_click",height = 250)),
                      #                 column(width = 6,
                      #                        plotOutput(outputId = "Hist_sta",click = "mymap_marker_click",height = 250))
             )
             )
    #     sidebarLayout(
    #        sidebarPanel(
    #           #  fluidRow(selectInput(inputId = "Model",label = "Choix du modèle",
    #           #        choices = c("Weibull","Fréchet","Exponentiel","MultRank")
    #           # )),
    #           fluidRow(uiOutput("sliders")),
    #           fluidRow(h5("Coûts des segmentations :"),
    #                    plotOutput(outputId = "Plot_part",height = 250)),
    #        fluidRow(h5("Informations cartographie :"),
    #                  plotOutput(outputId = "Plot_sta",click = "mymap_marker_click",height = 250))
    #     ),
    # mainPanel(
    #        fluidRow(leafletOutput("mymap", height=300),align = "center"),
    #        fluidRow(plotOutput("Plot_Seg", click = "plot_click", height=250),
    #                 h5("Informations sur le segment selectionné"),
    #                 verbatimTextOutput("Desc1_mod1"),
    #                 h5("Distribution du segment selectionné"),
    #                 plotOutput("Plot_Dens", height=250),
    #                 h5("Saisonnalité du segment selectionné"),
    #                 plotOutput("Sais_mod1", height=250))
    #     )
      # )
    ),
  tabPanel(title = "Readme",
           htmlOutput(outputId = "notice")
           # ,fluidRow(
           #   column(width = 4,
           #          selectInput(inputId = "Model_comp1",label = "Modélisation 1 :",
           #                      choices = c("Weibull","Fréchet","Exponentiel","MultRank")
           #          ),
           #          uiOutput("sliders_comp1")),
           #   column(width = 8,
           #          plotOutput("Plot_Seg_comp1", click = "plot_click1", height=150)
           #          )
           #   ),
           # fluidRow(
           #   column(width = 4,
           #          selectInput(inputId = "Model_comp2",label = "Modélisation 2 :",
           #                      choices = c("Weibull","Fréchet","Exponentiel","MultRank")
           #          ),
           #          uiOutput("sliders_comp2")),
           #   column(width = 8,
           #          plotOutput("Plot_Seg_comp2", click = "plot_click2", height=150)
           #          )
           #   ),
           # fluidRow(
           #   column(width = 12,
           #          align = "center",
           #          tableOutput(outputId = "table_comp"))
           # ),fluidRow(
           #   column(width = 6,
           #   leafletOutput("mymap_comp1", height=200)
           #     ),
           #   column(width = 6,
           #          leafletOutput("mymap_comp2", height=200)
           #   )
           #   )
           )
)

################## Server #################

 
server <- function(input, output) {
  
########################
#### PAGE D'ACCUEIL ####
########################
  
  output$Desc_mod1 <- renderText(
    paste(
      paste("Substance :", "Prosulfocarbe"),
      paste("Période d'observation :", as.character(min(input_max$d)),"-",as.character(max(input_max$d))),
      paste("Région d'observation :", "Centre-Val de Loire"),
      paste("Nombre de prélèvements :", as.character(nrow(data))),
      paste("% de quantification des prélèvements  :",as.character(round(sum(1-(data$col))/nrow(data)*100,3))),
      paste("Nombre de jours de prélèvements :", as.character(nrow(input_max))),
      paste("% de quantification des maximum journaliers :", as.character(round(sum(1-(input_max$col))/nrow(input_max)*100,3))),
      paste("Nombre de stations actives :", as.character(nrow(sta_crs_o))),
      sep = "\n"
    )
  )
  output$Plot_Max <- renderPlot(
    ggplot()+geom_point(data = input_max,
                        aes(x = d,
                            y = C,
                            col = as.factor(col==0))
                        )+
      theme(legend.position = "bottom")+
      scale_color_discrete("Quantification")+
      scale_x_continuous(name = "Date",
                         breaks = seq(min(input_max$d),max(input_max$d),by = "year"),
                         labels = year(seq(min(input_max$d),max(input_max$d),by = "year")))+
      scale_y_continuous(name = "Concentrations (µg/L)",
                         breaks = seq(0,max(input_max$C)+0.1),
                         labels = seq(0,max(input_max$C)+0.1))
  )
  output$Plot_comp <- renderPlot(
    g+geom_point(data = sta_crs_o
      ,aes(x = st_coordinates(geometry)[,1],
                     y = st_coordinates(geometry)[,2],
                     col = comp),
      size = 3)+
      scale_color_discrete("Composante hydrographique")+
      theme(legend.position = "bottom")
  )
  ########################
  #### PAGE DETECTION ####
  ########################  

   ts <- reactive({
          # if(length(which(crops_weibull$beta<=input$pen+10^-5)) == 0)
          # {
          #   ts_mod <- crops_weibull$segments[[1]]
          # }else
          # {
            ts_mod <- crops_weibull$segments[[max(which(crops_weibull$beta<=input$pen))]]
          # }
          ts_mod <- c(0,ts_mod,nrow(input_max))
          return(ts_mod)
       })
   tstp <- reactive({
           tstp_mod <- input_max$d[ts()[-1]]
           return(tstp_mod)
         })
   tstrt <- reactive({
             tstrt_mod <- input_max$d[(ts()+1)[-length(ts())]]
             return(tstrt_mod)
         })
   Ksel <- reactive({
       Ksel1 <- max(which(crops_weibull$beta<=input$pen+10^-5))
       return(Ksel1)
       })
   
   # slidinfo <- reactive({
   #       list("pen","Valeur de pénalité",
   #            crops_weibull$beta)
   #   })
            
    # output$sliders <- renderUI({
    #    inputID <- slidinfo()[[1]]
    #    inputName <- slidinfo()[[2]]
    #    sliderInput(inputID, inputName, min = round(min(slidinfo()[[3]])+10^-2,2),
    #                max = round(max(slidinfo()[[3]])+10^-2,2), value = round(min(slidinfo()[[3]])+10^-2,2),round = 2,step = 0.05
    #                )
    #    print(input$pen)
    #  })
    
    
    x_value <- reactiveValues()
    x_value$sel <- 1
    observeEvent(input$plot_click,{
      val <- as.Date(input$plot_click$x,origin = "1970-01-01")
      x_value$sel <- max(which(tstrt()<val))
    })
    observeEvent(input$pen,{
      x_value$sel <- 1
    })
    output$Plot_Seg <- renderPlot({
           ggplot()+geom_point(aes(x = input_max$d,y = seg()$log_C,col = as.factor(input_max$col == 0)))+
               geom_rect(aes(xmin = tstrt(),xmax = tstp(),ymin = min(seg()$log_C)-0.1,ymax = max(seg()$log_C)+0.1),col = "white",alpha = 0.3,show.legend = FALSE)+ 
               geom_rect(aes(xmin = tstrt()[x_value$sel],xmax = tstp()[x_value$sel],ymin = min(seg()$log_C)-0.1,ymax = max(seg()$log_C)+0.1),alpha = 0.0,col = "black",size = 1.5)+
               scale_color_discrete("Quantification")+
               ylim(c(min(seg()$log_C)-0.1,max(seg()$log_C)+0.1))+ylab("Concentrations (µg/L)")+xlab("Temps")+
               ggtitle("Résultat de segmentation")+
               theme(legend.position = "bottom")
    })
    x_value$part_abs <- unlist(lapply(crops_weibull$segments,FUN = "length"))
    x_value$part_ord <- crops_weibull$cost
    output$Plot_part <- renderPlot({
      ggplot()+geom_point(aes(x = x_value$part_abs,
                              y = x_value$part_ord))+
        geom_point(aes(x_value$part_abs[Ksel()],
                       x_value$part_ord[Ksel()]),
                   col = "red")+
        geom_abline(slope = ml1$coefficients[2],intercept = ml1$coefficients[1],col = "red")+
        geom_abline(slope = ml2$coefficients[2],intercept = ml2$coefficients[1],col = "red")+
        geom_vline(xintercept = x_value$part_abs[K_star])+
        xlab("Nombre de ruptures")+
        ylab("Coût de la segmentation")
    })
    
    seg <- reactive({
      pos <- which(input_max$d %in% seq(tstrt()[x_value$sel],tstp()[x_value$sel],by = "day"))
      print(c(min(input_max$d[pos]),max(input_max$d[pos])))
      tab <- input_max[pos,]
      lambda <- costfcpp(data = tab[,2:3],tstart = 0,tstop = nrow(tab)-1,init = "mle",q_init = 0.5,k = sigma,prec = 10^-8,Nmax = 100,lmax = 100000)[1]
      x_plot <- seq(0.001,max(tab$C),by = 0.001)
      a_plot <- unique(tab$LOQ)
      # print(a_plot)
      tab$weight <- NA
      p <- c()
      for(i in 1:length(a_plot))
      {
        p <- c(p,sum((tab$LOQ == a_plot[i])))
        tab$weight[tab$LOQ == a_plot[i]] <- p[i]
      }
      # print(p)
      # print(nrow(tab))
      p <- p/nrow(tab)
      tab$weight <- tab$weight/nrow(tab)
      y_plot <- rep(0,length(x_plot))
      for(i in 1:length(p))
      {
        p1 <- rep(0,length(x_plot))
        p1[x_plot >= a_plot[i]] <- (1-exp(-(lambda*x_plot[x_plot >= a_plot[i]])^(sigma)))
        # p1[x_plot > a_plot[i]] <- ((lambda*sigma)*(lambda*x_plot[x_plot > a_plot[i]])^(sigma-1)*exp(-(lambda*x_plot[x_plot > a_plot[i]])^sigma))*exp(-(lambda*a_plot[i])^(sigma))
        y_plot <- y_plot + p1 * p[i]
      }
      # y_plot <- y_plot#/sum(y_plot)
      # ay_plot <- y_plot[which(as.character(x_plot) %in% as.character(a_plot))]
      # y_plot[x_plot %in% a_plot] <- NA
      sais_C <- tab$C
      year_C <- rep(paste(as.character(unique(year(tab$d))),collapse = "-"),nrow(tab))
      span <- as.numeric(tab$d[nrow(tab)]-tab$d[1])
      if(span >= 365)
      {
        sais_C <- tab$C
      }else
      {
        span2 <- seq(tab$d[1],tab$d[nrow(tab)],by = "day")
        for(i in unique(year(input_max$d))) #[-length(unique(year(input_max$d)))]
        {
          add_sais <- span2[1]
          year(span2[1]) <- i
          add_sais <- seq(add_sais,by = "day",length.out = length(span2))
          sais_C <- c(sais_C,input_max$C[input_max$d %in% add_sais])
          year_C <- c(year_C,rep(paste(as.character(unique(year(add_sais))),collapse = "-"),length(input_max$C[input_max$d %in% add_sais])))
        }
      }
      sta_seg_tab <- sta_crs_o[,c("CdStationMesureEauxSurface","comp","geometry")]
      colnames(sta_seg_tab) <- c("sta","comp","geometry")
      sta_seg_tab <- st_transform(sta_seg_tab,4326) 
      pos <- which(data$d %in% seq(tstrt()[x_value$sel],tstp()[x_value$sel],by = "day"))
      sta_tab <- data[pos,]
      sta_tab <- sta_tab %>% group_by(sta) %>%
        summarize(N = n()) 
      sta_seg_tab$N <- 0
      sta_seg_tab$N[match(sta_tab$sta,sta_seg_tab$sta)] <- sta_tab$N
      # part_arg <- KK[[which(elb$V2 == input$Kluster)+1]]
      # part_arg <- mapply(FUN = "cbind",KK[[which(elb$V2 == input$Kluster)]],1:length(KK[[which(elb$V2 == input$Kluster)]]))
      # part_arg <- list.rbind(part_arg)
      # new_clust<- as.numeric(as.factor(paste(part_arg[,1],part_arg[,2],sep = "")))
      # names(new_clust) <- row.names(part_arg)
      # sta_seg_tab$clust <- new_clust[match(sta_seg_tab$sta,names(new_clust))]
      # sta_seg_tab$clust[is.na(sta_seg_tab$clust)] <- max(sta_seg_tab$clust,na.rm = TRUE)+1:length(which(is.na(sta_seg_tab$clust)))
      sta_seg_tab[["clust"]] <- sta_crs_o[[which(substr(x = colnames(sta_crs_o),start = nchar(colnames(sta_crs_o))-1,stop = nchar(colnames(sta_crs_o))) == as.character(input$Kluster))]][match(sta_seg_tab$sta,sta_crs_o$CdStationMesureEauxSurface)]
      # print(sta_seg_tab[["clust"]])
      C_MAT <- as.list(rep(NA,max(sta_seg_tab$clust)))
      for(k in 1:max(sta_seg_tab$clust))
      {
        pos <- which(sta_seg_tab$clust == k)
        sta <- sta_seg_tab[pos,]
        pos <- which(sta$N != 0)
        if(length(pos) > 0)
        {
          sta <- sta[pos,]
          mat_c <- matrix(NA,nc = nrow(sta),nr = nrow(sta))
          for(i in 1:nrow(mat_c))
          {
            for(j in 1:nrow(mat_c))
            {
              if(nrow(mat_c)>1)
              {
                # dij <- as.matrix(KLUST_MAT[[sta$comp[i]]])[labels(KLUST_MAT[[sta$comp[i]]])==sta$sta[i],labels(KLUST_MAT[[sta$comp[i]]])==sta$sta[j]]
                # 1/dij* 
                mat_c[i,j] <- 
                  wasserstein1d(a = data$C[(data$sta == sta$sta[i])&(data$d <= tstp()[x_value$sel])&(data$d >= tstrt()[x_value$sel])],
                                b = data$C[(data$sta == sta$sta[j])&(data$d <= tstp()[x_value$sel])&(data$d >= tstrt()[x_value$sel])],
                                p = 1)
              }else
              {
                mat_c[i,j] <- 0
              }
              diag(mat_c) <- 0
            }
          }
          row.names(mat_c) <- sta$sta
          colnames(mat_c) <- sta$sta
          C_MAT[[k]] <- mat_c
        }
      }
      sta_seg_tab$abs <- NA
      sta_seg_tab$ord <- NA
      for(i in 1:length(C_MAT))
      {
        if(is.matrix(C_MAT[[i]]))
        {
          pos <- which(sta_seg_tab$sta%in%row.names(C_MAT[[i]]))
          abs <- mean(apply(X = C_MAT[[i]],MARGIN = 1,FUN = "sum"))
          ord <- costfcpp(data = data[(data$sta %in% sta_seg_tab$sta[pos])&(data$d <= tstp()[x_value$sel])&(data$d >= tstrt()[x_value$sel]),2:3],
                          tstart = 0,
                          tstop = nrow(data[(data$sta %in% sta_seg_tab$sta[pos])&(data$d <= tstp()[x_value$sel])&(data$d >= tstrt()[x_value$sel]),2:3])-1,
                          init = "mle",q_init = 0.5,k = sigma,prec = 10^-8,
                          Nmax = 100,
                          lmax = 100000)[1]
          sta_seg_tab$abs[pos] <- abs
          sta_seg_tab$ord[pos] <- ord
        }
      }
      sta_seg_tab$ord <- 1/sta_seg_tab$ord
      # print(unique(sta_seg_tab$abs[which(sta_seg_tab$clust == 2)]))
      # D1 <- sta_seg_tab[which(!duplicated(sta_seg_tab$clust)),c("abs","ord","clust")]
      D1 <- data.table(matrix(NA,nrow = input$Kluster,ncol = 3))
      colnames(D1) <- c("abs","ord","clust")
      D1$clust <- 1:input$Kluster
      for(i in 1:nrow(D1))
      {
        v <- unique(sta_seg_tab$abs[which(sta_seg_tab$clust == i)])
        w <- unique(sta_seg_tab$ord[which(sta_seg_tab$clust == i)])
        pos <- which(!is.na(v))
        pos2 <- which(!is.na(w))
        if(length(pos) > 0)
        {
          D1$abs[i] <- v[pos]
          D1$ord[i] <- w[pos2]
        }
      }
      # print(D1)
      # D1 <- data.table(D1)
      # D1 <- D1[,geometry:=NULL]
      D1 <- D1[which(!is.na(D1$abs)),]
      D1$clust <- as.factor(D1$clust)
      p <- high(D1$abs) * high(D1$ord)
      D1 <- psel(D1, p, top = nrow(D1))
      D1$.level <- as.factor(D1$.level)
      sta_seg_tab$front <- NA
      sta_seg_tab$front <- D1$.level[match(sta_seg_tab$clust,D1$clust)]
      sta_seg_tab$front <- as.numeric(as.character(sta_seg_tab$front))
      if(input$inf == "Stations actives")
      {
        sta_seg_tab$to_plot <- sta_seg_tab$N != 0
        pal <- distinctColorPalette(k = 2)
        title <- "Activité"
      }else if(input$inf == "Composantes hydrographiques")
      {
        sta_seg_tab$to_plot <- sta_seg_tab$comp
        pal <- distinctColorPalette(k = length(unique(sta_seg_tab$to_plot)))
        title <- "Composante"
      }else if(input$inf == "Clusters spatiaux")
      {
        sta_seg_tab$to_plot <- sta_seg_tab$clust
        pal <- distinctColorPalette(k = length(unique(sta_seg_tab$to_plot)))
        title <- "Clusters"
      }else if(input$inf == "Valeurs du front de Pareto")
      {
        sta_seg_tab$to_plot <- sta_seg_tab$front
        # pal <- distinctColorPalette(k = length(unique(sta_seg_tab$to_plot[which(!is.na(sta_seg_tab$to_plot))])))
        # pal <- brewer.pal(n = length(unique(sta_seg_tab$to_plot[which(!is.na(sta_seg_tab$to_plot))])),name =  "YlOrRd")
        # pal <- rev(pal)
        pal <- viridis_pal(option = "C")(length(levels(as.factor(sta_seg_tab$to_plot))))
        pal <- rev(pal)
        title <- "Niveau de front de Pareto"
      }
      log_C = input_max$C 
      sais_C2 = sais_C
      if(input$button == TRUE)
      {
        log_C = log(log_C)
        sais_C = log(sais_C2)
      }else
      {
        log_C = input_max$C 
        sais_C = sais_C2
      }
      return(list(tab = tab,
                  log_C = log_C,
                  lambda = lambda,
                  x_plot = x_plot,
                  y_plot = y_plot,
                  a_plot = a_plot,
                  sais_C = sais_C,
                  year_C = year_C,
                  sta_seg_tab = sta_seg_tab,
                  pal = pal,
                  title = title
                  ))
    })
    
    output$Desc1_mod1 <- renderText(
      paste(
          paste("Dates des ruptures : \n",
                as.character(tstrt()[x_value$sel]),";",as.character(tstp()[x_value$sel])),
          paste("Nombre de données : ",as.character(nrow(seg()$tab))),
          paste("% de données quantifiées : ",as.character(round(sum(1-seg()$tab[3])/nrow(seg()$tab)*100,3))),
          paste("Nombre de stations actives : ",as.character(sum(seg()$sta_seg_tab$N!=0))),
          paste("Min : ",as.character(min(seg()$tab$C))),
          paste("Moyenne : ",as.character(mean(seg()$tab$C))),
          paste("Médiane : ",as.character(median(seg()$tab$C))),
          paste("Max : ",as.character(max(seg()$tab$C))),
          sep = "\n"
       )
    )
      
    output$Plot_Dens <- renderPlot({
        f <- ecdf(seg()$tab$C)
        ggplot()+
        geom_line(aes(x = seg()$x_plot,
                      y = seg()$y_plot,
                      col = "red"))+
          geom_step(aes(x = seg()$tab$C,
                        y = f(seg()$tab$C),
                        col = "blue"))+
          geom_vline(aes(xintercept = seg()$a_plot,
                         col = "black"),show.legend = T)+
          ggtitle("Fonctions de répartition")+
          scale_color_identity(name = "",
                               breaks = c("red", "blue","black"),
                               labels = c("Théorique", "Empirique","LOQ"),
                               guide = "legend")+
          xlab("")+
          ylab("")+
          theme(legend.position = "bottom",
                legend.direction = "horizontal")
       })

       output$Sais_mod1 <- renderPlot(
          ggplot()+ylim(c(min(seg()$sais_C)-0.1,max(seg()$sais_C)+0.1))+
             geom_violin(aes(y = seg()$sais_C,
                                    x = seg()$year_C),
                                trim = FALSE)+
               stat_ydensity(aes(y = seg()$sais_C,
                                 x = seg()$year_C))+
               geom_violin(aes(y = seg()$sais_C[1:length((ts()[x_value$sel]+1):ts()[x_value$sel+1])],
                               x = seg()$year_C[1:length((ts()[x_value$sel]+1):ts()[x_value$sel+1])]),
                           col = "red",trim = FALSE)+
               stat_ydensity(aes(y = seg()$sais_C[1:length((ts()[x_value$sel]+1):ts()[x_value$sel+1])],
                                 x = seg()$year_C[1:length((ts()[x_value$sel]+1):ts()[x_value$sel+1])]),col = "red")+
               xlab("Années")+
               ylab("Concentrations (µg/L)")+
               theme(axis.text.x = element_text(angle = 90))
       )
    
       output$mymap <- renderLeaflet({
                   m
                })
       
       observeEvent({input$pen
         input$inf
         input$plot_click
         input$Kluster
         input$mymap_marker_click
         1
       },
       {
         if(input$inf == "Stations actives")
         {
           leafletProxy("mymap") %>%
             clearMarkers() %>%
             clearControls() %>%
             addCircleMarkers(layerId = seg()$sta_seg_tab$sta,
                              lng = st_coordinates(seg()$sta_seg_tab)[,1],
                              lat = st_coordinates(seg()$sta_seg_tab)[,2],
                              color = "black",
                              stroke = TRUE,
                              weight = 1,
                              opacity = 1,
                              fillOpacity = 1,
                              fillColor = seg()$pal[as.numeric(as.factor(seg()$sta_seg_tab$to_plot))],
                              radius = 4,
                              popup = seg()$sta_seg_tab$sta)%>%
             addCircleMarkers(layerId = seg()$sta_seg_tab$sta[which(seg()$sta_seg_tab$clust == map_clicker$clust_sel)],
                              lng = st_coordinates(seg()$sta_seg_tab)[which(seg()$sta_seg_tab$clust == map_clicker$clust_sel),1],
                              lat = st_coordinates(seg()$sta_seg_tab)[which(seg()$sta_seg_tab$clust == map_clicker$clust_sel),2],
                              color = "red",
                              stroke = TRUE,
                              weight = 1,
                              opacity = 1,
                              fillOpacity = 1,
                              fillColor = seg()$pal[as.numeric(as.factor(seg()$sta_seg_tab$to_plot[which(seg()$sta_seg_tab$clust == map_clicker$clust_sel)]))],
                              radius = 4,
                              popup = seg()$sta_seg_tab$sta[which(seg()$sta_seg_tab$clust == map_clicker$clust_sel)]) %>%
             addLegend(position = 'topright',
                       colors = seg()$pal,
                       labels = levels(as.factor(seg()$sta_seg_tab$to_plot)),
                       values = seg()$pal,
                       opacity = 1,
                       title = "Stations actives")
         }else if(input$inf == "Composantes hydrographiques")
         {
           leafletProxy("mymap") %>%
                  clearMarkers() %>%
                  clearControls() %>%
                  addCircleMarkers(layerId = seg()$sta_seg_tab$sta[which(seg()$sta_seg_tab$N == 0)],
                                   lng = st_coordinates(seg()$sta_seg_tab)[which(seg()$sta_seg_tab$N == 0),1],
                                   lat = st_coordinates(seg()$sta_seg_tab)[which(seg()$sta_seg_tab$N == 0),2],
                                             color = "black",
                                             stroke = TRUE,
                                             weight = 1,
                                             opacity = 1,
                                             fillOpacity = 1,
                                             fillColor = "black",
                                             radius = 2,
                                             popup = seg()$sta_seg_tab$sta[which(seg()$sta_seg_tab$N == 0)]) %>%
                            addCircleMarkers(layerId = seg()$sta_seg_tab$sta[which(seg()$sta_seg_tab$N != 0)],
                                             lng = st_coordinates(seg()$sta_seg_tab)[which(seg()$sta_seg_tab$N != 0),1],
                                                     lat = st_coordinates(seg()$sta_seg_tab)[which(seg()$sta_seg_tab$N != 0),2],
                                                     color = "black",
                                                     stroke = TRUE,
                                                     weight = 1,
                                                     opacity = 1,
                                                     fillOpacity = 1,
                                                     fillColor = seg()$pal[as.numeric(as.factor(seg()$sta_seg_tab$to_plot[which(seg()$sta_seg_tab$N != 0)]))],
                                                     radius = 4,
                                                     popup = seg()$sta_seg_tab$sta[which(seg()$sta_seg_tab$N != 0)])%>%
             addCircleMarkers(layerId = seg()$sta_seg_tab$sta[which(seg()$sta_seg_tab$clust == map_clicker$clust_sel)],
                              lng = st_coordinates(seg()$sta_seg_tab)[which(seg()$sta_seg_tab$clust == map_clicker$clust_sel),1],
                              lat = st_coordinates(seg()$sta_seg_tab)[which(seg()$sta_seg_tab$clust == map_clicker$clust_sel),2],
                              color = "red",
                              stroke = TRUE,
                              weight = 1,
                              opacity = 1,
                              fillOpacity = 1,
                              fillColor = seg()$pal[seg()$sta_seg_tab$to_plot[which(seg()$sta_seg_tab$clust == map_clicker$clust_sel)]],
                              radius = 4,
                              popup = seg()$sta_seg_tab$sta[which(seg()$sta_seg_tab$clust == map_clicker$clust_sel)]) %>%
                                    addLegend(position = 'topright',
                                              colors = seg()$pal,
                                              labels = levels(as.factor(seg()$sta_seg_tab$to_plot)),
                                              values = seg()$pal,
                                              opacity = 1,
                                              title = "Composante")
         }else if(input$inf == "Clusters spatiaux")
         {
           leafletProxy("mymap") %>%
             clearMarkers() %>%
             clearControls() %>%
             addCircleMarkers(layerId = seg()$sta_seg_tab$sta[which(seg()$sta_seg_tab$N == 0)],
                              lng = st_coordinates(seg()$sta_seg_tab)[which(seg()$sta_seg_tab$N == 0),1],
                              lat = st_coordinates(seg()$sta_seg_tab)[which(seg()$sta_seg_tab$N == 0),2],
                              color = "black",
                              stroke = TRUE,
                              weight = 1,
                              opacity = 1,
                              fillOpacity = 1,
                              fillColor = "black",
                              radius = 2,
                              popup = seg()$sta_seg_tab$sta[which(seg()$sta_seg_tab$N == 0)]) %>%
             addCircleMarkers(layerId = seg()$sta_seg_tab$sta[which(seg()$sta_seg_tab$N != 0)],
                              lng = st_coordinates(seg()$sta_seg_tab)[which(seg()$sta_seg_tab$N != 0),1],
                              lat = st_coordinates(seg()$sta_seg_tab)[which(seg()$sta_seg_tab$N != 0),2],
                              color = "black",
                              stroke = TRUE,
                              weight = 1,
                              opacity = 1,
                              fillOpacity = 1,
                              fillColor = seg()$pal[as.numeric(as.factor(seg()$sta_seg_tab$to_plot[which(seg()$sta_seg_tab$N != 0)]))],
                              radius = 4,
                              popup = seg()$sta_seg_tab$sta[which(seg()$sta_seg_tab$N != 0)])%>%
             addCircleMarkers(layerId = seg()$sta_seg_tab$sta[which(seg()$sta_seg_tab$clust == map_clicker$clust_sel)],
                              lng = st_coordinates(seg()$sta_seg_tab)[which(seg()$sta_seg_tab$clust == map_clicker$clust_sel),1],
                              lat = st_coordinates(seg()$sta_seg_tab)[which(seg()$sta_seg_tab$clust == map_clicker$clust_sel),2],
                              color = "red",
                              stroke = TRUE,
                              weight = 1,
                              opacity = 1,
                              fillOpacity = 1,
                              fillColor = seg()$pal[seg()$sta_seg_tab$to_plot[which(seg()$sta_seg_tab$clust == map_clicker$clust_sel)]],
                              radius = 4,
                              popup = seg()$sta_seg_tab$sta[which(seg()$sta_seg_tab$clust == map_clicker$clust_sel)]) %>%
             addLegend(position = 'topright',
                       colors = seg()$pal,
                       labels = levels(as.factor(seg()$sta_seg_tab$to_plot)),
                       values = seg()$pal,
                       opacity = 1,
                       title = "Clusters")
         }else if(input$inf == "Valeurs du front de Pareto")
         {
           leafletProxy("mymap") %>%
             clearMarkers() %>%
             clearControls() %>%
             addCircleMarkers(layerId = seg()$sta_seg_tab$sta[which(seg()$sta_seg_tab$N == 0)],
                              lng = st_coordinates(seg()$sta_seg_tab)[which(seg()$sta_seg_tab$N == 0),1],
                              lat = st_coordinates(seg()$sta_seg_tab)[which(seg()$sta_seg_tab$N == 0),2],
                              color = "black",
                              stroke = TRUE,
                              weight = 1,
                              opacity = 1,
                              fillOpacity = 1,
                              fillColor = "black",
                              radius = 2,
                              popup = seg()$sta_seg_tab$sta[which(is.na(seg()$sta_seg_tab$to_plot))]) %>%
             addCircleMarkers(layerId = seg()$sta_seg_tab$sta[which(seg()$sta_seg_tab$N != 0)],
                              lng = st_coordinates(seg()$sta_seg_tab)[which(seg()$sta_seg_tab$N != 0),1],
                              lat = st_coordinates(seg()$sta_seg_tab)[which(seg()$sta_seg_tab$N != 0),2],
                              color = "black",
                              stroke = TRUE,
                              weight = 1,
                              opacity = 1,
                              fillOpacity = 1,
                              fillColor = seg()$pal[as.factor(seg()$sta_seg_tab$to_plot)[which(seg()$sta_seg_tab$N != 0)]],
                              radius = 4,
                              popup = seg()$sta_seg_tab$sta[which(seg()$sta_seg_tab$N != 0)])%>%
             addCircleMarkers(layerId = seg()$sta_seg_tab$sta[which(seg()$sta_seg_tab$clust == map_clicker$clust_sel)],
                              lng = st_coordinates(seg()$sta_seg_tab)[which(seg()$sta_seg_tab$clust == map_clicker$clust_sel),1],
                              lat = st_coordinates(seg()$sta_seg_tab)[which(seg()$sta_seg_tab$clust == map_clicker$clust_sel),2],
                              color = "red",
                              stroke = TRUE,
                              weight = 1,
                              opacity = 1,
                              fillOpacity = 1,
                              fillColor = seg()$pal[as.factor(seg()$sta_seg_tab$to_plot)[which(seg()$sta_seg_tab$clust == map_clicker$clust_sel)]],
                              radius = 4,
                              popup = seg()$sta_seg_tab$sta[which(seg()$sta_seg_tab$clust == map_clicker$clust_sel)]) %>%
             addLegend(position = 'topright',
                       colors = seg()$pal,
                       labels = levels(as.factor(seg()$sta_seg_tab$to_plot)),
                       values = seg()$pal,
                       opacity = 1,
                       title = "Niveau du front de Pareto")
         }
       })
       
       map_clicker <- reactiveValues()
       map_clicker$d <- as.Date("1970-01-01")
       map_clicker$C <- 0
       map_clicker$col <- as.factor(1)
       
       
         observeEvent({
           input$mymap_marker_click
           }, 
           {
             # print(input$mymap_marker_click)
              clickId <- input$mymap_marker_click$id
              n <- data[which(data$d %in% seq(tstrt()[x_value$sel],tstp()[x_value$sel],by = "day")),] 
              map_clicker$clust_sel <- unique(seg()$sta_seg_tab$clust[seg()$sta_seg_tab$sta == clickId]) 
              kl <- seg()$sta_seg_tab$sta[which(seg()$sta_seg_tab$clust == map_clicker$clust_sel)]
              map_clicker$klC <- n[which(n$sta %in% kl),c("sta","C","col","LOQ")]
              map_clicker$mostq <- c(NA,NA,NA)
              temp <- map_clicker$klC%>%group_by(sta)%>%summarize(N=n(),q = round(sum(1-col)/n()*100,2))
              pos <- which.max(temp$q)
              map_clicker$mostq <- c(temp$sta[pos],temp$q[pos],temp$N[pos]) 
              pos <- which(n$sta == clickId)
              n <- n[pos,]
              X_sta <- n$d
              Y_sta <- as.numeric(as.character(n$C))
              col_sta <- n$col
              map_clicker$Id <- clickId
              map_clicker$d <- X_sta
              map_clicker$C <- Y_sta
              map_clicker$col <- as.factor(1-col_sta)
          })
         
         output$Klust_info <- renderText({
           mod <- paste(c("Valeurs de LOQ : ",
                          paste(as.character(sort(unique(map_clicker$klC$LOQ))),collapse = ", ")),
                        collapse = "")
           mod2 <- names(table(map_clicker$klC$LOQ))[which.max(table(map_clicker$klC$LOQ))]
           pos <- str_locate(string = mod,pattern = mod2)[2]
           mod <- paste(substring(mod, c(1,pos+1), c(pos,nchar(mod))), collapse="*")
           paste("Informations sur les concentrations dans le cluster selectionné",
                 paste("Nombre de données : ",
                       as.character(length(map_clicker$klC$C))),
                 paste("% de quantification : ",
                       as.character(round(sum(1-map_clicker$klC$col)/length(map_clicker$klC$col)*100,2))),
                 paste("Nombre de stations : ",
                       as.character(length(unique(map_clicker$klC$sta)))),
                 paste("Min : ",
                       as.character(min(map_clicker$klC$C))),
                 paste("Mean : ",
                       as.character(mean(map_clicker$klC$C))),
                 paste("Median : ",
                       as.character(median(map_clicker$klC$C))),
                 paste("Max : ",
                       as.character(max(map_clicker$klC$C))),
                 mod,
                 paste("Station la plus quantifiée du cluster : ",
                       as.character(max(map_clicker$mostq[1])),"avec",as.character(max(map_clicker$mostq[2])),"% de quantification sur",as.character(max(map_clicker$mostq[3])),"données."),
                 sep = "\n"
           )
         })
         
         output$Plot_pareto <- renderPlot({
           cols <- rep(0,length(seg()$sta_seg_tab$abs))
           cols[seg()$sta_seg_tab$clust == map_clicker$clust_sel] <- 1
           ggplot()+geom_point(aes(x = seg()$sta_seg_tab$abs,
                                   y = seg()$sta_seg_tab$ord,
                                   fill = as.factor(seg()$sta_seg_tab$front),
                                   col = cols),
                               shape=21, 
                               stroke = 2,
                               size = 4)+
             scale_fill_viridis(discrete = TRUE,option = "C",direction = -1)+labs(fill = "Valeurs du front de Pareto")+
             theme(legend.position = "bottom")+ylab("Paramètre estimé")+xlab("Homogénéité des concentrations des stations")+
             ggtitle("Analyse des clusters")+
             scale_color_gradient(low = "black",high = "red")+
             guides(colour = "none")
             # guides(guide_legend(list(color = "#000000")))
           # +
           #   geom_point(aes(x = seg()$sta_seg_tab$abs[seg()$sta_seg_tab$clust == map_clicker$clust_sel],
           #                  y = seg()$sta_seg_tab$ord[seg()$sta_seg_tab$clust == map_clicker$clust_sel],
           #                  fill = as.factor(seg()$sta_seg_tab$front[seg()$sta_seg_tab$clust == map_clicker$clust_sel])),
           #              color = "red",
           #              stroke = 2,
           #              shape=21, 
           #              size = 4)
         })
         
           output$Plot_sta <- renderPlot({
                ggplot()+geom_point(aes(x = as.Date(map_clicker$d),
                                        y = map_clicker$C,
                                        col = map_clicker$col))+
                    xlim(c(tstrt()[x_value$sel],tstp()[x_value$sel]))+
                    ylim(c(min(map_clicker$C),max(map_clicker$C)))+
                    scale_color_discrete('Quantification')+
                    xlab("Temps")+ylab("Concentration")+
                    ggtitle(paste("Concentration de la station",map_clicker$Id))+
               theme(legend.position = "bottom")
            })

           
#########################
##### ONGLET NOTICE #####  
#########################
           

           output$notice <- renderUI(renderDocument(x = htmlTemplate(filename = "notice.html")))
           
           
           
       # observeEvent({input$plot_click
       #   input$Kluster
       #   input$inf},
       #   {
       #     leafletProxy("mymap") %>%
       #         clearMarkers() %>%
       #         clearControls() %>%
       #         addCircleMarkers(lng = st_coordinates(seg()$sta_seg_tab)[,1],
       #                          lat = st_coordinates(seg()$sta_seg_tab)[,2],
       #                          color = "black",
       #                          stroke = TRUE,
       #                          weight = 1,
       #                          opacity = 1,
       #                          fillOpacity = 1,
       #                          fillColor = seg()$pal[match(seg()$sta_seg_tab$to_plot,unique(seg()$sta_seg_tab$to_plot))],
       #                          radius = 4) %>%
       #         addLegend(position = 'topright',
       #                   colors = seg()$pal,
       #                   labels = unique(seg()$sta_seg_tab$to_plot),
       #                   values = unique(match(seg()$sta_seg_tab$to_plot,unique(seg()$sta_seg_tab$to_plot))),
       #                   opacity = 1,
       #                   title = "Activité")
       #     # if(input$inf == "Composantes hydrographiques")
       #     # {
       #     #   pal <- distinctColorPalette(k = length(KLUST_MAT))
       #     #   leafletProxy("mymap") %>%
       #     #     clearMarkers() %>%
       #     #     clearControls() %>%
       #     #     addCircleMarkers(lng = st_coordinates(seg()$sta_seg_tab)[,1],
       #     #                      lat = st_coordinates(seg()$sta_seg_tab)[,2],
       #     #                      color = "black",
       #     #                      stroke = TRUE,
       #     #                      weight = 1,
       #     #                      opacity = 1,
       #     #                      fillOpacity = 1,
       #     #                      fillColor = pal[seg()$sta_seg_tab$comp],
       #     #                      radius = 4) %>%
       #     #     addLegend(position = 'topright',
       #     #               pal = pal,
       #     #               values = unique(seg()$sta_seg_tab$comp),
       #     #               title = "Composante hydrographique",
       #     #               opacity = 1)
       #     # }
       #      
       #     # %>%
       #     #            addCircleMarkers(layerId = x_value$info1,
       #     #                             lng = st_coordinates(coord_sta)[x_value$station,1],
       #     #                             lat = st_coordinates(coord_sta)[x_value$station,2],
       #     #                             color = "black",
       #     #                             stroke = TRUE,
       #     #                             weight = 1,
       #     #                             opacity = 1,
       #     #                             fillColor = x_value$lambdaCol(x_value$info3),
       #     #                             fillOpacity = 1,
       #     #                             radius = sqrt(x_value$info2*7.5)*(!is.na(x_value$info3))+sqrt(x_value$info2*2)*(is.na(x_value$info3)),
       #     #                             popup = paste("Code station :",as.character(x_value$info1),"<br>",
       #     #                                           "n =",as.character(x_value$info2),"<br>",
       #     #                                           "n_quant :",as.character(x_value$info6),"<br>",
       #     #                                           "lambda =",as.character(x_value$popupinfo3),"<br>",
       #     #                                           "Seg_Hydro :",as.character(x_value$seg_hydro))) %>%
       #     #            addLegend('topright', pal = x_value$lambdaCol, values = x_value$info3,
       #     #                        title = x_value$title,
       #     #                        opacity = 1)
       # })
       
 #  m <- leaflet() %>%
 #    setView(lng = (max(st_coordinates(coord_sta)[,1])+min(st_coordinates(coord_sta)[,1]))/2,
 #            lat = (max(st_coordinates(coord_sta)[,2])+min(st_coordinates(coord_sta)[,2]))/2,
 #            zoom = 5) %>%
 #    addTiles() %>%
 #    addPolygons(data = map_reg_her,
 #                color = "black",
 #                fillOpacity = 0,
 #                opacity = 1,
 #                stroke = 0.5,
 #                weight = 3) %>%
 #    addPolylines(data = st_zm(x = map_eau,drop = T,what = "ZM"),
 #                 weight = 2,color = "darkblue"
 #    ) %>%
 #    addFlows(lng0 = as.numeric(st_coordinates(flow_dir)[1:nrow(flow_dir)%%2 == 0,1]),
 #             lat0 = as.numeric(st_coordinates(flow_dir)[1:nrow(flow_dir)%%2 == 0,2]),
 #             lng1 = as.numeric(st_coordinates(flow_dir)[1:nrow(flow_dir)%%2 == 1,1]),
 #             lat1 = as.numeric(st_coordinates(flow_dir)[1:nrow(flow_dir)%%2 == 1,2]),
 #             dir = -1,minThickness = 1.2,maxThickness = 1.2,color = "darkblue",opacity = 1)
 #    
 #    ##### Onglet 1 #####
 #  
 #    tab <- input_max[1:(crops_weibull$segments[[1]][1]-1),]
 #    pos <- which(data$d %in% seq(input_max$d[1],input_max$d[crops_weibull$segments[[1]][1]-1],by = "day"))
 #    tab_2 <- data[pos,]
 #    a <- unique(tab$C[tab$col == 1])
 #    p <- c()
 #    for(i in 1:length(a))
 #    {
 #      p <- c(p,sum((tab$C == a[i])&(tab$col==1)))
 #    }
 #    p <- (p+sum(tab$col==0)/length(a))/nrow(tab)
 #    absi <- seq(0.001,max(tab$C),length.out = max(tab$C)*1000)
 #    dens_mel <- rep(0,length(absi))
 #    p1 <- rep(0,length(absi))
 #    p1[absi == a[i]] <- (1-exp(-(l_init*a[i])^(0.5)))
 #    p1[absi > a[i]] <- (l_init*0.5)*(l_init*absi[absi > a[i]])^(-0.5)*exp(-(l_init*absi[absi > a[i]])^0.5)
 #    dens_mel <- dens_mel + p1 * p[i]
 #    x_value <- reactiveValues()
 #    x_value$sel <- 1
 #    x_value$lambda  <- costfcpp(data = tab,tstart = 0,tstop = nrow(tab)-1,init = "mean",q_init = 0.5,k = 0.5,prec = 10^-8,Nmax = 100,lmax = 10)[1]
 #    x_value$station <- which(coord_sta$id_station %in% tab_2$sta)
 #    x_value$info1 <- inf_map[,1]
 #    x_value$info2 <- inf_map[,2]
 #    x_value$info3 <- inf_map[,3]
 #    x_value$info4 <- "50"
 #    x_value$info6 <-  inf_map[,4]
 #    x_value$lambdaCol <- colorNumeric(palette = 'Oranges', inf_map[,3])
 #    x_value$y_plot <- dens_mel
 #    x_value$x_plot <- absi
 #    x_value$sais_C <- tab$C
 #    x_value$year_C <- rep("2007",nrow(tab_2))
 #    x_value$title <- "Distance entre le lambda station<br>et le lambda segment"
 #    x_value$popupinfo3 <- ""
 #    x_value$a_plot <- a
 #      # unique(data$C[ts3[1]:(ts3[2]-1)][data$col[ts3[1]:(ts3[2]-1)] == 1])
 #    x_value$seg_hydro <- ""
 #    x_value$part_abs <- unlist(lapply(crops_weibull$segments,FUN = "length"))
 #    x_value$part_ord <- crops_weibull$cost
 #    x_value$Kopt <- 8
 #    map_clicker <- reactiveValues()
 #    map_clicker$d <- tab_2$d[tab_2$sta == tab_2$sta[1]]
 #    map_clicker$C <- tab_2$C[tab_2$sta == tab_2$sta[1]]
 #    map_clicker$col <- as.factor(1-tab_2$col[tab_2$sta == tab_2$sta[1]])
 # 
 #   ts <- reactive({
 #          if(input$Model == "Fréchet")
 #          {
 #            if(length(which(crops_frechet$beta<=input$pen+10^-5)) == 0)
 #            {
 #              ts_mod <- crops_frechet$segments[[1]]
 #            }else
 #            {
 #              ts_mod <- crops_frechet$segments[[max(which(crops_frechet$beta<=input$pen+10^-5))]]
 #            }
 #          }else if(input$Model == "Weibull")
 #          {
 #            if(length(which(crops_weibull$beta<=input$pen+10^-5)) == 0)
 #            {
 #              ts_mod <- crops_weibull$segments[[1]]
 #            }else
 #            {
 #              ts_mod <- crops_weibull$segments[[max(which(crops_weibull$beta<=input$pen+10^-5))]]
 #            }
 #          }else if(input$Model == "Exponentiel")
 #          {
 #            if(length(which(crops_exp$beta<=input$pen+10^-5)) == 0)
 #            {
 #              ts_mod <- crops_exp$segments[[1]]
 #            }else
 #            {
 #              ts_mod <- crops_exp$segments[[max(which(crops_exp$beta<=input$pen+10^-5))]]
 #            }
 #          }else if(input$Model == "MultRank")
 #          {
 #            if(length(which(1:length(ts_Multrank)<=input$pen+0.1)) == 0)
 #            {
 #              ts_mod <- ts_Multrank[[1]]
 #            }else
 #            {
 #              ts_mod <- ts_Multrank[[floor(input$pen+0.1)]] 
 #            }
 #          }
 #          ts_mod <- c(0,ts_mod,nrow(input_max))
 #          return(ts_mod)
 #      })
 #   
 #   tstp <- reactive({
 #          tstp_mod <- input_max$d[ts()[-1]]
 #          return(tstp_mod)
 #        })
 #   
 #   tstrt <- reactive({
 #            tstrt_mod <- input_max$d[(ts()+1)[-length(ts())]]
 #            return(tstrt_mod)
 #        })
 #    
 #   slidinfo <- reactive({
 #      if(input$Model == "Fréchet")
 #      {
 #        list("pen","Valeur de pénalité",
 #             crops_frechet$beta)
 #      }else if(input$Model == "Weibull")
 #      {
 #        list("pen","Valeur de pénalité",
 #             crops_weibull$beta)
 #      }else if(input$Model == "Exponentiel")
 #      {
 #        list("pen","Valeur de pénalité",
 #             crops_exp$beta)
 #      }else if(input$Model == "MultRank")
 #      {
 #        list("pen","Nombre de ruptures",
 #             1:length(ts_Multrank))
 #      }
 #    })
 #    
 #   output$sliders <- renderUI({
 #      inputID <- slidinfo()[[1]]
 #      inputName <- slidinfo()[[2]]
 #      sliderInput(inputID, inputName, min = min(slidinfo()[[3]]), 
 #                  max = max(slidinfo()[[3]])+10^-5, value = min(slidinfo()[[3]])+10^-5,round = -2,step = 0.05
 #                  )
 #    })
 #   
 #   observeEvent(input$Model,{
 #                     if(input$Model == "Weibull")
 #                     {
 #                       p_a <- unlist(lapply(crops_weibull$segments,FUN = "length"))
 #                       p_o <- crops_weibull$cost
 #                       val5 <- "Paramètre de la station"
 #                       val6 <- "50"
 #                     }else if(input$Model == "Fréchet")
 #                     {
 #                       p_a <- unlist(lapply(crops_frechet$segments,FUN = "length"))
 #                       p_o <- crops_frechet$cost
 #                       val5 <- "Paramètre de la station"
 #                       val6 <- "50"
 #                     }else if(input$Model == "Exponentiel")
 #                     {
 #                       p_a <- unlist(lapply(crops_exp$segments,FUN = "length"))
 #                       p_o <- crops_exp$cost
 #                       val5 <- "Paramètre de la station"
 #                       val6 <- "50"
 #                     }else
 #                     {
 #                       p_a <- 1:length(ts_Multrank)
 #                       p_o <- rep(NA,length(ts_Multrank))
 #                       for(i in 1:length(p_o))
 #                       {
 #                         p_o[i] <- sum(mapply(FUN = "costnp",c(1,ts_Multrank[[i]]+1),c(ts_Multrank[[i]],nrow(input_max)),MoreArgs = list(sumstat = as.data.frame(input_max[2:3]))))
 #                       }
 #                       val5 <- "Log P valeur du test de Wilcoxon"
 #                       val6 <- "Non applicable"
 #                     }
 #                     K_opt <- rep(NA,length(p_a)-2)
 #                     for(i in 2:(length(p_a)-1))
 #                     {
 #                       ml1 <- lm(p_o[1:i]~p_a[1:i])
 #                       ml2 <- lm(p_o[i:length(p_a)]~p_a[i:length(p_a)])
 #                       K_opt[i-1] <- sum(ml1$residuals^2)+sum(ml2$residuals^2)
 #                     }
 #                     K_opt <- p_a[(2:(length(p_a)-1))][which.min(K_opt)]
 #                     # print(K_opt)
 #                     x_value$Kopt <- K_opt 
 #                     x_value$part_abs <- p_a
 #                     x_value$part_ord <- p_o
 #                     x_value$title <- val5
 #                     x_value$info4 <- val6
 #                 })
 #    
 #   observeEvent(input$plot_click, {
 #        val <- as.Date(input$plot_click$x,origin = "1970-01-01")
 #        x_value$sel <- max(which(tstrt()<val))
 #        if(input$Model == "Exponentiel")
 #        {
 #          sourceCpp("RCPPpropre3.cpp")
 #          val2 <- costfcpp(data = input_max[,2:3],tstart = ts()[x_value$sel],tstop = ts()[x_value$sel+1]-1,init = "mean",q_init = 0.5,prec = 10^-8,Nmax = 100,lmax = 100)[1]
 #        }else if(input$Model == "Fréchet")
 #        {
 #          sourceCpp("Frechetbis.cpp")
 #          val2 <- costfcpp(data = input_max[,2:3],tstart = ts()[x_value$sel],tstop = ts()[x_value$sel+1]-1,q_init = 0.75,prec = 10^-8,Nmax = 100,lmax = 10)[1]
 #        }else if(input$Model == "Weibull")
 #        {
 #          sourceCpp("Weibull.cpp")
 #          val2 <- costfcpp(data = input_max[,2:3],tstart = ts()[x_value$sel],tstop = ts()[x_value$sel+1]-1,init = "mean",q_init = 0.5,prec = 10^-8,Nmax = 100,lmax = 10,k = 1/2)[1]
 #        }else if(input$Model == "MultRank")
 #        {
 #          val2 <- "Non applicable"
 #        }
 #        x_value$lambda <- val2
 #        val3 <- which(coord_sta$id_station %in% data$sta[data$d %in% seq(tstrt()[x_value$sel],tstp()[x_value$sel],by = "day")]) #(ts()[x_value$sel]+1):ts()[x_value$sel+1]
 #        x_value$station <- val3
 #        x_value$info1 <- coord_sta$id_station[val3]
 #        # n <- data[(ts()[x_value$sel]+1):ts()[x_value$sel+1],] %>% group_by(sta) %>% summarise(n = n(),n_quant = sum(1-col),cdloc = unique(cdloc))
 #        n <- data[which(data$d %in% seq(tstrt()[x_value$sel],tstp()[x_value$sel],by = "day")),] %>% group_by(sta) %>% summarise(n = n(),n_quant = sum(1-col),cdloc = unique(cdloc))
 #        x_value$info2 <- n$n[match(x_value$info1,n$sta)]
 #        x_value$info6 <- n$n_quant[match(x_value$info1,n$sta)]
 #        val_hydro <- n$cdloc[match(x_value$info1,n$sta)]
 #        # n <- data[(ts()[x_value$sel]+1):ts()[x_value$sel+1],]
 #        n <- data[which(data$d %in% seq(tstrt()[x_value$sel],tstp()[x_value$sel],by = "day")),]
 #        val4 <- rep(0,length(x_value$info1))
 #        if(input$Model == "Fréchet")
 #        {
 #          for(i in 1:length(val4))
 #          {
 #            if(sum(1-n$col[n$sta == x_value$info1[i]]) != 0)
 #            {
 #              val4[i] <- costfcpp(data = n[n$sta == x_value$info1[i],],tstart = 0,tstop = x_value$info2[i]-1,q_init = 0.75,prec = 10^-8,Nmax = 100,lmax = 10)[1]
 #            }else
 #            {
 #              val4[i] <- NA
 #            }
 #          }
 #        }else if(input$Model == "Exponentiel")
 #        {
 #          for(i in 1:length(val4))
 #          {
 #            if(sum(1-n$col[n$sta == x_value$info1[i]]) != 0)
 #            {
 #              val4[i] <- costfcpp(data = n[n$sta == x_value$info1[i],],tstart = 0,tstop = x_value$info2[i]-1,init = "mean",q_init = 1/2,prec = 10^-8,Nmax = 100,lmax = 100)[1]
 #            }else
 #            {
 #              val4[i] <- NA
 #            }
 #          }
 #        }else if(input$Model == "Weibull")
 #        {
 #          for(i in 1:length(val4))
 #          {
 #            if(sum(1-n$col[n$sta == x_value$info1[i]]) != 0)
 #            {
 #              val4[i] <- costfcpp(data = n[n$sta == x_value$info1[i],],tstart = 0,tstop = x_value$info2[i]-1,init = "mean",q_init = 0.5,prec = 10^-8,Nmax = 100,lmax = 10,k = 1/2)[1]
 #            }else
 #            {
 #              val4[i] <- NA
 #            }
 #          }
 #        }else
 #        {
 #          for(i in 1:length(val4))
 #          {
 #            val4[i] <- UStatcens(data = c(n$C[n$sta == x_value$info1[i]],n$C[n$sta != x_value$info1[i]]),
 #                                 cens = c(1-n$col[n$sta == x_value$info1[i]],1-n$col[n$sta != x_value$info1[i]]),
 #                                 t = sum(n$sta == x_value$info1[i]))
 #          }
 #        }
 #        # if(input$Model != "MultRank")
 #        # {
 #        #     val4 <- abs(x_value$lambda - val4)
 #        # }
 #        x_value$seg_hydro <- val_hydro
 #        pale <- brewer.pal(n = 9,name =  "Oranges")
 #        pale_rev <- rev(brewer.pal(9, "Oranges"))
 #        if(input$Model == "MultRank")
 #        {
 #            x_value$info3 <- log(1-pchisq(val4,df=1))
 #            x_value$popupinfo3 <- log(1-pchisq(val4,df=1))
 #            x_value$lambdaCol <- colorNumeric(palette = pale_rev,x_value$info3)
 #        }
 #        else if(input$Model == "Fréchet")
 #        {
 #          x_value$info3 <- val4
 #          x_value$popupinfo3 <- val4
 #          x_value$lambdaCol <- colorNumeric(palette = pale, x_value$info3)
 #        }
 #        else{
 #          x_value$info3 <- val4
 #          x_value$popupinfo3 <- val4
 #          x_value$lambdaCol <- colorNumeric(palette = pale_rev, x_value$info3)
 #        }
 #        #x_value$lambdaCol <- colorBin( "plasma", bins=quantile(x_value$info3,c(0,0.05,0.5,0.95,1)), na.color = "#aaff56")
 #        # density_info <- data[(ts()[x_value$sel]+1):ts()[x_value$sel+1],]
 #        density_info <- input_max[(ts()[x_value$sel]+1):ts()[x_value$sel+1],]
 #        a <- unique(density_info$C[density_info$col == 1])
 #        if(input$Model != "MultRank")
 #        {
 #            # density_info <- data[ts()[x_value$sel]:(ts()[x_value$sel+1]-1)]
 #            # a <- unique(density_info$C[density_info$col == 1])
 #            p <- c()
 #            for(i in 1:length(a))
 #            {
 #               p <- c(p,sum((density_info$C == a[i])&(density_info$col==1)))
 #            }
 #            p <- (p+sum(density_info$col==0)/length(a))/nrow(density_info)
 #            absi <- seq(0.001,max(input_max$C[(ts()[x_value$sel]+1):ts()[x_value$sel+1]]),length.out = max(input_max$C[(ts()[x_value$sel]+1):ts()[x_value$sel+1]])*1000)
 #            dens_mel <- rep(0,length(absi))
 #            for(i in 1:length(p))
 #            {
 #                p1 <- rep(0,length(absi))
 #                if(input$Model == "Fréchet")
 #                {
 #                    p1[absi == a[i]] <- exp(-a[i]^(-x_value$lambda))
 #                    p1[absi > a[i]] <- x_value$lambda*absi[absi > a[i]]^(-1-x_value$lambda)*exp(-absi[absi > a[i]]^(-x_value$lambda))
 #                    dens_mel <- dens_mel + p1 * p[i]
 #                }else if(input$Model == "Exponentiel")
 #                {
 #                    p1[absi == a[i]] <- (1-exp(-x_value$lambda*a[i]))
 #                    p1[absi > a[i]] <- x_value$lambda*exp(-x_value$lambda*absi[absi > a[i]])
 #                    dens_mel <- dens_mel + p1 * p[i]
 #                }else if(input$Model == "Weibull")
 #                {
 #                  p1[absi == a[i]] <- (1-exp(-(x_value$lambda*a[i])^(0.5)))
 #                  p1[absi > a[i]] <- (x_value$lambda*0.5)*(x_value$lambda*absi[absi > a[i]])^(-0.5)*exp(-(x_value$lambda*absi[absi > a[i]])^0.5)
 #                  dens_mel <- dens_mel + p1 * p[i]
 #                }
 #            }
 #        }else if(input$Model == "MultRank")
 #        {
 #            absi <- seq(0.001,max(input_max$C[(ts()[x_value$sel]+1):ts()[x_value$sel+1]]),length.out = max(input_max$C[(ts()[x_value$sel]+1):ts()[x_value$sel+1]])*1000)
 #            dens_mel <-  rep(0,length(absi))
 #        }
 #        x_value$a_plot <- a
 #        x_value$x_plot <- absi
 #        x_value$y_plot <- dens_mel
 #        n <- input_max[(ts()[x_value$sel]+1):ts()[x_value$sel+1],]
 #        span <- as.numeric(n$d[nrow(n)]-n$d[1])
 #        if(span >= 365)
 #        {
 #             x_value$sais_C <- n$C
 #             x_value$year_C <- rep(paste(as.character(unique(year(n$d))),collapse = "-"),nrow(n))
 #        }else
 #        {
 #             span2 <- seq(n$d[1],n$d[nrow(n)],by = "day")
 #             x_value$sais_C <- n$C
 #             x_value$year_C <- rep(paste(as.character(unique(year(n$d))),collapse = "-"),nrow(n))
 #             for(i in unique(year(data$d))[-length(unique(year(data$d)))])
 #             {
 #                add_sais <- span2[1]
 #                year(span2[1]) <- i
 #                add_sais <- seq(add_sais,by = "day",length.out = length(span2))
 #                x_value$sais_C <- c(x_value$sais_C,data$C[data$d %in% add_sais])
 #                x_value$year_C <- c(x_value$year_C,rep(paste(as.character(unique(year(add_sais))),collapse = "-"),length(data$C[data$d %in% add_sais])))
 #             }
 #        }
 #        # print(x_value$a_plot)
 #        # print(x_value$y_plot[x_value$x_plot %in% x_value$a_plot])
 #        leafletProxy("mymap") %>%
 #            clearMarkers() %>%
 #            clearControls() %>%
 #            addCircleMarkers(lng = st_coordinates(coord_sta)[-x_value$station,1],
 #                             lat = st_coordinates(coord_sta)[-x_value$station,2],
 #                             color = colorNumeric(palette = 'Greys',0),
 #                             fillOpacity = 0.25,
 #                             radius = 2) %>%
 #            addCircleMarkers(layerId = x_value$info1,
 #                             lng = st_coordinates(coord_sta)[x_value$station,1],
 #                             lat = st_coordinates(coord_sta)[x_value$station,2],
 #                             color = "black",
 #                             stroke = TRUE,
 #                             weight = 1,
 #                             opacity = 1,
 #                             fillColor = x_value$lambdaCol(x_value$info3),
 #                             fillOpacity = 1,
 #                             radius = sqrt(x_value$info2*7.5)*(!is.na(x_value$info3))+sqrt(x_value$info2*2)*(is.na(x_value$info3)),
 #                             popup = paste("Code station :",as.character(x_value$info1),"<br>",
 #                                           "n =",as.character(x_value$info2),"<br>",
 #                                           "n_quant :",as.character(x_value$info6),"<br>",
 #                                           "lambda =",as.character(x_value$popupinfo3),"<br>",
 #                                           "Seg_Hydro :",as.character(x_value$seg_hydro))) %>%
 #            addLegend('topright', pal = x_value$lambdaCol, values = x_value$info3,
 #                        title = x_value$title,
 #                        opacity = 1)
 #    })
 #  
 #   observeEvent(input$mymap_marker_click, {
 #        clickId <- input$mymap_marker_click$id
 #        n <- data[which(data$d %in% seq(tstrt()[x_value$sel],tstp()[x_value$sel],by = "day")),]  # (ts()[x_value$sel]+1):ts()[x_value$sel+1]
 #        pos <- which(n$sta == clickId)
 #        n <- n[pos,]
 #        X_sta <- n$d
 #        Y_sta <- n$C
 #        col_sta <- n$col
 #        map_clicker$d <- X_sta
 #        map_clicker$C <- Y_sta
 #        map_clicker$col <- as.factor(1-col_sta)
 #    })
 # 
 #   output$Plot_sta <- renderPlot({
 #        ggplot()+geom_point(aes(x = as.Date(map_clicker$d),
 #                                y = map_clicker$C,
 #                                col = map_clicker$col))+
 #            xlim(c(tstrt()[x_value$sel],tstp()[x_value$sel]))+
 #            xlab("Temps")+ylab("Concentration")+labs(color='Quantification')+
 #            ggtitle("Concentration de la station selectionnée ")
 #    })
 # 
 #   output$Plot_Seg <- renderPlot({
 #        ggplot()+geom_point(aes(x = input_max$d,y = input_max$C,col = Quantification))+
 #            geom_rect(aes(xmin = tstrt(),xmax = tstp(),ymin = -0.1,ymax = 3),col = "white",alpha = 0.3,show.legend = FALSE)+ #as.character(2:length(ts())%%2)
 #            geom_rect(aes(xmin = tstrt()[x_value$sel],xmax = tstp()[x_value$sel],ymin = -0.1,ymax = 3),alpha = 0.0,col = "black")+
 #            ylim(c(-0.1,3))+ylab("Concentrations (µg/L)")+xlab("Temps")+
 #            ggtitle(paste("Ruptures obtenues par la modèlisation",as.character(input$Model)))
 #    })
 #  
 #   output$Plot_Dens <- renderPlot({
 #        ggplot()+geom_line(aes(x = x_value$x_plot,
 #                               y = x_value$y_plot,
 #                               col = "red"))+
 #            geom_density(aes(x = input_max$C[(ts()[x_value$sel]+1):ts()[x_value$sel+1]],
 #                             y = ..density../length((ts()[x_value$sel]+1):ts()[x_value$sel+1])*15,
 #                             col = "blue"),stat = "density")+
 #            geom_rug(aes(x = input_max$C[(ts()[x_value$sel]+1):ts()[x_value$sel+1]]))+
 #            geom_vline(aes(xintercept = x_value$a_plot,col = "black"),show.legend = T)+
 #            scale_color_identity(name = "Densité",
 #                                 breaks = c("red", "blue","black"),
 #                                 labels = c("Théorique", "Empirique","Valeurs de LQ"),
 #                                 guide = "legend")+
 #            xlab("")+
 #            ylab("")+
 #            theme(legend.position = "bottom",legend.direction = "horizontal")
 #    })
 # 
 #   output$Plot_part <- renderPlot({
 #     ggplot()+geom_point(aes(x = x_value$part_abs,
 #                            y = x_value$part_ord))+
 #       geom_vline(xintercept = x_value$Kopt)+
 #       xlab("Nombre de ruptures")+
 #       ylab("Coût de la segmentation")
 #   })
 #   
 #   output$mymap <- renderLeaflet({
 #           m
 #        })
 # 
 #   output$Desc1_mod1 <- renderText(
 #    paste(
 #        paste("Dates des ruptures : \n",
 #              as.character(tstrt()[x_value$sel]),";",as.character(tstp()[x_value$sel])),
 #        paste("Valeur du lambda : \n",
 #              as.character(x_value$lambda)),
 #        paste("Nombre de données : \n",as.character(sum(x_value$info2))),
 #        paste("% de données quantifiées : \n",as.character(round(sum(1-data$col[ts()[x_value$sel]:(ts()[x_value$sel+1]-1)])/length(ts()[x_value$sel]:(ts()[x_value$sel+1]-1))*100,3))),
 #        sep = "\n"
 #    )
 # )
 #    output$Desc_mod1 <- renderText(
 #        paste(
 #            paste("Nombre de données (N) :", as.character(nrow(data))),
 #            paste("% de données quantifiées :",as.character(round(sum(1-(data$col))/nrow(data)*100,3))),
 #            paste("Taille de segment minimale :",x_value$info4),
 #            sep = "\n"
 #        )
 #    )
 # 
 #    output$Sais_mod1 <- renderPlot(
 #       ggplot()+ylim(c(-0.1,max(x_value$sais_C)+0.1))+
 #          geom_violin(aes(y = x_value$sais_C,
 #                                 x = x_value$year_C),
 #                             trim = FALSE)+
 #            stat_ydensity(aes(y = x_value$sais_C,
 #                              x = x_value$year_C,))+
 #            geom_violin(aes(y = x_value$sais_C[1:length((ts()[x_value$sel]+1):ts()[x_value$sel+1])],
 #                            x = x_value$year_C[1:length((ts()[x_value$sel]+1):ts()[x_value$sel+1])]),
 #                        col = "red",trim = FALSE)+
 #            stat_ydensity(aes(y = x_value$sais_C[1:length((ts()[x_value$sel]+1):ts()[x_value$sel+1])],
 #                              x = x_value$year_C[1:length((ts()[x_value$sel]+1):ts()[x_value$sel+1])]),col = "red")+
 #            xlab("Années")+
 #            ylab("Concentrations (µg/L)")
 #    )
 #    
 #    ##### Onglet 2 #####
 #    
 #    slidinfo_comp1 <- reactive({
 #      if(input$Model_comp1 == "Fréchet")
 #      {
 #        list("pen_comp1","Valeur de pénalité",
 #             crops_frechet$beta)
 #      }else if(input$Model_comp1 == "Weibull")
 #      {
 #        list("pen_comp1","Valeur de pénalité",
 #             crops_weibull$beta)
 #      }else if(input$Model_comp1 == "Exponentiel")
 #      {
 #        list("pen_comp1","Valeur de pénalité",
 #             crops_exp$beta)
 #      }else if(input$Model_comp1 == "MultRank")
 #      {
 #        list("pen_comp1","Nombre de ruptures",
 #             1:length(ts_Multrank))
 #      }
 #    })
 #    
 #    output$sliders_comp1 <- renderUI({
 #      inputID <- slidinfo_comp1()[[1]]
 #      inputName <- slidinfo_comp1()[[2]]
 #      sliderInput(inputID, inputName, min = min(slidinfo_comp1()[[3]]), 
 #                  max = max(slidinfo_comp1()[[3]])+10^-5, value = min(slidinfo_comp1()[[3]])+10^-5,round = -2,step = 0.05
 #      )
 #    })
 #    
 #    ts_comp1 <- reactive({
 #      if(input$Model_comp1 == "Fréchet")
 #      {
 #        if(length(which(crops_frechet$beta<=input$pen_comp1+10^-5)) == 0)
 #        {
 #          ts_mod_comp1 <- crops_frechet$segments[[1]]
 #        }else
 #        {
 #          ts_mod_comp1 <- crops_frechet$segments[[max(which(crops_frechet$beta<=input$pen_comp1+10^-5))]]
 #        }
 #      }else if(input$Model_comp1 == "Weibull")
 #      {
 #        if(length(which(crops_weibull$beta<=input$pen_comp1+10^-5)) == 0)
 #        {
 #          ts_mod_comp1 <- crops_weibull$segments[[1]]
 #        }else
 #        {
 #          ts_mod_comp1 <- crops_weibull$segments[[max(which(crops_weibull$beta<=input$pen_comp1+10^-5))]]
 #        }
 #      }else if(input$Model_comp1 == "Exponentiel")
 #      {
 #        if(length(which(crops_exp$beta<=input$pen_comp1+10^-5)) == 0)
 #        {
 #          ts_mod_comp1 <- crops_exp$segments[[1]]
 #        }else
 #        {
 #          ts_mod_comp1 <- crops_exp$segments[[max(which(crops_exp$beta<=input$pen_comp1+10^-5))]]
 #        }
 #      }else if(input$Model_comp1 == "MultRank")
 #      {
 #        if(length(which(1:length(ts_Multrank)<=input$pen_comp1+0.1)) == 0)
 #        {
 #          ts_mod_comp1 <- ts_Multrank[[1]]
 #        }else
 #        {
 #          ts_mod_comp1 <- ts_Multrank[[floor(input$pen_comp1+0.1)]] 
 #        }
 #      }
 #      ts_mod_comp1 <- c(0,ts_mod_comp1,nrow(data))
 #      return(ts_mod_comp1)
 #    })
 #    
 #    slidinfo_comp2 <- reactive({
 #      if(input$Model_comp2 == "Fréchet")
 #      {
 #        list("pen_comp2","Valeur de pénalité",
 #             crops_frechet$beta)
 #      }else if(input$Model_comp2 == "Weibull")
 #      {
 #        list("pen_comp2","Valeur de pénalité",
 #             crops_weibull$beta)
 #      }else if(input$Model_comp2 == "Exponentiel")
 #      {
 #        list("pen_comp2","Valeur de pénalité",
 #             crops_exp$beta)
 #      }else if(input$Model_comp2 == "MultRank")
 #      {
 #        list("pen_comp2","Nombre de ruptures",
 #             1:length(ts_Multrank))
 #      }
 #    })
 #    
 #    output$sliders_comp2 <- renderUI({
 #      inputID <- slidinfo_comp2()[[1]]
 #      inputName <- slidinfo_comp2()[[2]]
 #      sliderInput(inputID, inputName, min = min(slidinfo_comp2()[[3]]), 
 #                  max = max(slidinfo_comp2()[[3]])+10^-5, value = min(slidinfo_comp2()[[3]])+10^-5,round = -2,step = 0.05
 #      )
 #    })
 #    
 #    ts_comp2 <- reactive({
 #      if(input$Model_comp2 == "Fréchet")
 #      {
 #        if(length(which(crops_frechet$beta<=input$pen_comp2+10^-5)) == 0)
 #        {
 #          ts_mod_comp2 <- crops_frechet$segments[[1]]
 #        }else
 #        {
 #          ts_mod_comp2 <- crops_frechet$segments[[max(which(crops_frechet$beta<=input$pen_comp2+10^-5))]]
 #        }
 #      }else if(input$Model_comp2 == "Weibull")
 #      {
 #        if(length(which(crops_weibull$beta<=input$pen_comp2+10^-5)) == 0)
 #        {
 #          ts_mod_comp2 <- crops_weibull$segments[[1]]
 #        }else
 #        {
 #          ts_mod_comp2 <- crops_weibull$segments[[max(which(crops_weibull$beta<=input$pen_comp2+10^-5))]]
 #        }
 #      }else if(input$Model_comp2 == "Exponentiel")
 #      {
 #        if(length(which(crops_exp$beta<=input$pen_comp2+10^-5)) == 0)
 #        {
 #          ts_mod_comp2 <- crops_exp$segments[[1]]
 #        }else
 #        {
 #          ts_mod_comp2 <- crops_exp$segments[[max(which(crops_exp$beta<=input$pen_comp2+10^-5))]]
 #        }
 #      }else if(input$Model_comp2 == "MultRank")
 #      {
 #        if(length(which(1:length(ts_Multrank)<=input$pen_comp2+0.1)) == 0)
 #        {
 #          ts_mod_comp2 <- ts_Multrank[[1]]
 #        }else
 #        {
 #          ts_mod_comp2 <- ts_Multrank[[floor(input$pen_comp2+0.1)]] 
 #        }
 #      }
 #      ts_mod_comp2 <- c(0,ts_mod_comp2,nrow(data))
 #      return(ts_mod_comp2)
 #    })
 #    
 #    x_value$sel_comp1 <- 1
 #    x_value$sel_comp2 <- 1
 #    x_value$l_seg1  <- 200
 #    x_value$l_seg2  <- 200
 #    x_value$q_seg1 <- 0
 #    x_value$q_seg2 <- 0
 #    x_value$d_comp1 <- paste(as.character(data$d[1]),data$d[crops_weibull$segments[[1]][1]],sep=" : ")
 #    x_value$d_comp2 <- paste(as.character(data$d[1]),data$d[crops_weibull$segments[[1]][1]],sep=" : ")
 #    x_value$p_comp1 <- round(costfcpp(data = data,tstart = 0,tstop = 198,init = "mean",q_init = 1/2,k = 1/2,prec = 10^-8,Nmax = 100,lmax = 10)[1],2)
 #    x_value$p_comp2 <- round(costfcpp(data = data,tstart = 0,tstop = 198,init = "mean",q_init = 1/2,k = 1/2,prec = 10^-8,Nmax = 100,lmax = 10)[1],2)
 #      
 #    observeEvent({
 #      input$plot_click1
 #      #input$Model_comp1
 #      },{
 #      val <- as.Date(input$plot_click1$x,origin = "1970-01-01")
 #      x_value$sel_comp1 <- max(which(data$d[ts_comp1()[-length(ts_comp1())]+1]<val))
 #      x_value$d_comp1 <- paste(as.character(data$d[ts_comp1()[x_value$sel_comp1]+1]),
 #                               as.character(data$d[ts_comp1()[x_value$sel_comp1+1]]),
 #                               sep = " : ")
 #      x_value$l_seg1 <- length((ts_comp1()[x_value$sel_comp1]+1):ts_comp1()[x_value$sel_comp1+1])
 #      x_value$q_seg1 <- sum(1-data$col[(ts_comp1()[x_value$sel_comp1]+1):ts_comp1()[x_value$sel_comp1+1]])
 #      if(input$Model_comp1 == "Weibull")
 #      {
 #        sourceCpp("Weibull.cpp")
 #        val <- costfcpp(data = data,tstart = ts_comp1()[x_value$sel_comp1],tstop = ts_comp1()[x_value$sel_comp1+1]-1,init = "mean",q_init = 1/2,k = 1/2,prec = 10^-8,Nmax = 100,lmax = 10)[1]
 #        val <- round(val,2)
 #      }else if(input$Model_comp1 == "Fréchet")
 #      {
 #        sourceCpp("Frechetbis.cpp")
 #        val <- costfcpp(data = data,tstart = ts_comp1()[x_value$sel_comp1],tstop = ts_comp1()[x_value$sel_comp1+1]-1,q_init = 0.75,prec = 10^-8,Nmax = 100,lmax = 10)[1]
 #        val <- round(val,2)
 #      }else if(input$Model_comp1 == "Exponentiel")
 #      {
 #        sourceCpp("RCPPpropre3.cpp")
 #        val <- costfcpp(data = data,tstart = ts_comp1()[x_value$sel_comp1],tstop = ts_comp1()[x_value$sel_comp1+1]-1,init = "mean",q_init = 0.5,prec = 10^-8,Nmax = 100,lmax = 100)[1]
 #        val <- round(val,2)
 #      }else if(input$Model_comp1 == "MultRank")
 #      {
 #        # val <- UStatcens(data = c(data$C[(ts_comp1()[x_value$sel_comp1]+1):ts_comp1()[x_value$sel_comp1+1]],data$C[-c((ts_comp1()[x_value$sel_comp1]+1):ts_comp1()[x_value$sel_comp1+1])]),
 #        #                  cens = c(data$col[(ts_comp1()[x_value$sel_comp1]+1):ts_comp1()[x_value$sel_comp1+1]],data$col[-c((ts_comp1()[x_value$sel_comp1]+1):ts_comp1()[x_value$sel_comp1+1])]),
 #        #                  t = length((ts_comp1()[x_value$sel_comp1]+1):ts_comp1()[x_value$sel_comp1+1]))
 #        val <- "Non Applicable"
 #      }
 #      x_value$p_comp1 <- val
 #    }
 #  )
 # 
 #    observeEvent(input$plot_click2,{
 #                   val <- as.Date(input$plot_click2$x,origin = "1970-01-01")
 #                   x_value$sel_comp2 <- max(which(data$d[ts_comp2()[-length(ts_comp2())]+1]<val))
 #                   x_value$d_comp2 <- paste(as.character(data$d[ts_comp2()[x_value$sel_comp2]+1]),
 #                                            as.character(data$d[ts_comp2()[x_value$sel_comp2+1]]),
 #                                            sep = " : ")
 #                   x_value$l_seg2 <- length((ts_comp2()[x_value$sel_comp2]+1):ts_comp2()[x_value$sel_comp2+1])
 #                   x_value$q_seg2 <- sum(1-data$col[(ts_comp2()[x_value$sel_comp2]+1):ts_comp2()[x_value$sel_comp2+1]])
 #                   if(input$Model_comp2 == "Weibull")
 #                   {
 #                     sourceCpp("Weibull.cpp")
 #                     val <- costfcpp(data = data,tstart = ts_comp2()[x_value$sel_comp2],tstop = ts_comp2()[x_value$sel_comp2+1]-1,init = "mean",q_init = 1/2,k = 1/2,prec = 10^-8,Nmax = 100,lmax = 10)[1]
 #                     val <- round(val,2)
 #                   }else if(input$Model_comp2 == "Fréchet")
 #                   {
 #                     sourceCpp("Frechetbis.cpp")
 #                     val <- costfcpp(data = data,tstart = ts_comp2()[x_value$sel_comp2],tstop = ts_comp2()[x_value$sel_comp2+1]-1,q_init = 0.75,prec = 10^-8,Nmax = 100,lmax = 10)[1]
 #                     val <- round(val,2)
 #                   }else if(input$Model_comp2 == "Exponentiel")
 #                   {
 #                     sourceCpp("RCPPpropre3.cpp")
 #                     val <- costfcpp(data = data,tstart = ts_comp2()[x_value$sel_comp2],tstop = ts_comp2()[x_value$sel_comp2+1]-1,init = "mean",q_init = 0.5,prec = 10^-8,Nmax = 100,lmax = 100)[1]
 #                     val <- round(val,2)
 #                   }else if(input$Model_comp2 == "MultRank")
 #                   {
 #                     # val <- UStatcens(data = c(data$C[(ts_comp2()[x_value$sel_comp2]+1):ts_comp2()[x_value$sel_comp2+1]],data$C[-(ts_comp2()[x_value$sel_comp2]+1):ts_comp2()[x_value$sel_comp2+1]]),
 #                     #                  cens = c(data$col[(ts_comp2()[x_value$sel_comp2]+1):ts_comp2()[x_value$sel_comp2+1]],data$col[-(ts_comp2()[x_value$sel_comp2]+1):ts_comp2()[x_value$sel_comp2+1]]),
 #                     #                  t = length((ts_comp2()[x_value$sel_comp2]+1):ts_comp2()[x_value$sel_comp2+1])
 #                     # )
 #                     val <- "Non applicable"
 #                   }
 #                   x_value$p_comp2 <- val
 #                 })
 #    
 #    output$Plot_Seg_comp1 <- renderPlot({
 #      ggplot()+geom_point(aes(x = data$d,y = data$C,col = Quantification))+
 #        geom_rect(aes(xmin = data$d[ts_comp1()[-length(ts_comp1())]+1],xmax = data$d[ts_comp1()[-1]],ymin = -0.1,ymax = 3),col = "white",alpha = 0.3,show.legend = FALSE)+ #as.character(2:length(ts())%%2)
 #        geom_rect(aes(xmin = data$d[ts_comp1()[x_value$sel_comp1]+1],xmax = data$d[ts_comp1()[x_value$sel_comp1+1]],ymin = -0.1,ymax = 3),alpha = 0.0,col = "black")+
 #        ylim(c(-0.1,3))+ylab("Concentrations (µg/L)")+xlab("Temps")+
 #        ggtitle(paste("Ruptures obtenues par la modèlisation",as.character(input$Model)))
 #    })
 #    
 #    output$Plot_Seg_comp2 <- renderPlot({
 #      ggplot()+geom_point(aes(x = data$d,y = data$C,col = Quantification))+
 #        geom_rect(aes(xmin = data$d[ts_comp2()[-length(ts_comp2())]+1],xmax = data$d[ts_comp2()[-1]],ymin = -0.1,ymax = 3),col = "white",alpha = 0.3,show.legend = FALSE)+ #as.character(2:length(ts())%%2)
 #        geom_rect(aes(xmin = data$d[ts_comp2()[x_value$sel_comp2]+1],xmax = data$d[ts_comp2()[x_value$sel_comp2+1]],ymin = -0.1,ymax = 3),alpha = 0.0,col = "black")+
 #        ylim(c(-0.1,3))+ylab("Concentrations (µg/L)")+xlab("Temps")+
 #        ggtitle(paste("Ruptures obtenues par la modèlisation",as.character(input$Model)))
 #    })
 #    
 #    output$table_comp <- renderTable({
 #        tab <- data.table(matrix(c("Date",x_value$d_comp1,x_value$d_comp2,
 #                                   "Lambda",as.character(x_value$p_comp1),as.character(x_value$p_comp2),
 #                                   "N",as.character(x_value$l_seg1),as.character(x_value$l_seg2),
 #                                   "N_quant",as.character(x_value$q_seg1),as.character(x_value$q_seg2)),byrow = TRUE,ncol = 3,nrow = 4))
 #        colnames(tab) <- c("","Segment 1","Segment 2")
 #        tab
 #      }
 #    )
    
    # output$mymap_comp1 <- renderLeaflet({
    #   m
    # })
    # 
    # output$mymap_comp2 <- renderLeaflet({
    #   m
    # })
    
    
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)



