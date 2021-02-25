
#install.packages(c("shiny", "shinyWidgets", "sf", "stars", "gstat", "plotrix", "rhandsontable", "shinyalert", "ggplot2"))
library(shinyalert)
library(shiny)
library(shinyWidgets)
library(sf)
library(stars)
library(gstat)
library(plotrix)
library(rhandsontable)
library(ggplot2)



# Define UI for application
ui <- fixedPage(
  titlePanel("", windowTitle = "Krigzy"),
  
  tags$head(
    tags$style(
      HTML(
        ".selectize-input, .selectize-dropdown {
        padding: 2px 2px;
        min-height: 1px;
        min-width: 1px;
       }",
        "input[type=\"number\"] {
        height: 25px;
        padding: 2px 2px;
        min-height: 1px;
        min-width: 1px;
      }",
        "#gen_col {
        margin-top: 2em;
      }",
        "#nr {
        width: 30px;
        position: relative;
        right: 90%;
        bottom: -5px;
        }",
        "#n_p {
        margin-top: -30px;
        }",
        "#table {
        color: black;
        }",
        "#p_text {
        position: relative;
        bottom: -28px;
        }",
        "#grid_input_row {
        margin-bottom: -11px;
        }",
        "#ex, #ey {
        position: relative;
        width: 40px;
        right: -15px;
        bottom: 25px;
        }",
        "#e_text, #exey {
        position: relative;
        top: 15px;
        }",
        "#model {
        text-align: left;
        }",
        ".irs-max, .irs-min {
        color: black;
        background-color: transparent;
        margin-top: 4em;
        }",
        ".irs-single:after {
        transform: rotate(180deg);
        top: -5.5px;
        }",
        ".irs-single {
        top: 42px;
        }",
        ".irs {
        margin-top: -30px;
        }",
        ".irs-line {
        overflow: visible;
        }",
        ".shiny-input-container:not(.shiny-input-container-inline) {
        width: auto;
        }",
        "input#width.form-control.shiny-bound-input, input#height.form-control.shiny-bound-input {
        margin-bottom: -20px;
        }",
        "#estimation, #variance {
        background-color: white;
        padding: 7px;
        border: solid #e6e3ff 3px;
        font-size: 20px;
        }",
        "label.control-label {
        font-weight: normal;
        }",
        "#krigzy {
        position: absolute;
        bottom: 0;
        padding: 10%;
        }",
        "#krigzydiv {
        position: relative;
        height: 300px;
        }",
        "#author {
        font-size: 12px;
        font-family: monospace;
        font-style: italic;
        padding-left: 10%;
        margin-top: -10px;
        }",
        "#app_name {
        font-size: 12px;
        font-family: monospace;
        font-style: italic;
        margin-top: -8%;
        padding-left: 10%;
        font-weight: bold;
        }"
      ))),
  
        setBackgroundImage("bg.png"),
        chooseSliderSkin("Flat", color = "purple"),
        useShinyalert(),
    
        fixedRow(
            
            column(4, 
                   h4("Grid"),
                   div(plotOutput("siatka", height = "auto"), style = "border: solid #6800a9 3px;"),
                   br(),
                   fixedRow(id = "grid_input_row",
                       column(8,
                              fixedRow(
                                column(6,
                                       numericInput("width", "Width", 15, 4, Inf, 1)),
                                column(6,
                                       numericInput("height", "Height", 15, 4, Inf, 1))
                              ),
                              fixedRow(
                                column(7,
                                       selectInput("n", "", choices = list("nmax", "maxdist"), selected = "nmax")),
                                column(5,
                                       numericInput("nval", "", 5, 1, Inf, 1))
                              )
                       ),
                       column(4, align = "center", id = "gen_col",
                              actionButton("gen_grid", HTML("Update<br/>grid"), style = "padding: 10%; margin-left: -10px;")
                       )
                   ),
                   hr(),
                   fixedRow(id = "n_p",
                     column(6, 
                            p("Known points:", id = "p_text")
                            ),
                     column(3,
                            numericInput("nr", "", 5, 2, 7, 1)
                            )
                   ),
                   rHandsontableOutput("table"),
                   p("Estimated raster:", id = "e_text"),
                   fixedRow(id = "exey",
                     column(3,
                            numericInput("ex", "X", 8, 1, Inf, 1)),
                     column(3,
                            numericInput("ey", "Y", 7, 1, Inf, 1))
                     ),
                   br(),
                   textOutput("mv"),
                   textOutput("mv2")
                   ),
            
            column(4,
                   h4("Semivariogram", style = "text-align: left;"),
                   plotOutput("plot", height = "auto"),
                   br(),
                   fixedRow(
                     column(6, selectInput("model", "Model", choices = list("Spherical", "Exponential", "Gaussian", "Power"), selected = "Spherical"))
                   ),
                   sliderInput("nugget", "Nugget", min = 0, max = 60, value = 30, ticks = F),
                   br(),
                   sliderInput("psill", "Psill", min = 1, max = 60, value = 30, ticks = F),
                   br(),
                   sliderInput("range", "Range", min = 1, max = 11.2, value = 5.6, ticks = F)),
            
            column(4,
                   h4("Kriging weights"),
                   div(plotOutput("weights", height = "auto")),
                   br(),
                   fixedRow(
                     column(6, align = "right",
                            h5("Prediction value"),
                            textOutput("estimation")),
                     column(6, align = "right",
                            h5("Prediction variance"),
                            textOutput("variance"))
                            ),
                   br(),
                   div(imageOutput("krigzy", height = "auto"), id = "krigzydiv"),
                   p("Krigzy", id = "app_name"),
                   p("by Hubert Kwaśny (themercerus)", id = "author")
                   )
            )
)


server <- function(input, output, session) {
    
  output$krigzy = renderImage({ list(src = "www/krigzy_alpha.png", width = "70%") },
                              deleteFile = F)
  
    ###########creating function nmax/maxdist##########
  
    n_function = function(n_choice) {
      
      euc_dist = as.vector(st_distance(st_as_sf(r_values$data, coords = c("X", "Y")),
                                       st_as_sf(r_values$est_coords, coords = c("X", "Y"))))

      if (n_choice == "maxdist") {
        
        i = which(euc_dist <= input$nval)
        r_values$chosen = r_values$data[i,]

      } else

      if (n_choice == "nmax") {

        sorted_points = sort(euc_dist, index.return = T)
        r_values$chosen = r_values$data[sort(sorted_points$ix[1:input$nval]),]

      }

    }
    
    ##########tworzenie funkcji liczenia wag i wyników estymacji
    
    wev = function(known, unknown, model) {
      
      if ((model$model[-1] == "Pow" && model$range[-1] > 2) || (input$n == "nmax" && input$nval > input$nr)) 
        {}
      else {
        
          known_sf = st_as_sf(known, coords = c("X", "Y"))
          unknown_sf = st_as_sf(unknown, coords = c("X", "Y"))
          
          #  print(known_sf)
          #  print(unknown_sf)
          #  print(model)
          # print(known)
          
          dA = st_distance(known_sf)
          
          if (dim(dA)[1] == 0) {shinyalert("No points in the specified distance", "Change points coordinates or increase the distance.",
                                           closeOnClickOutside = T, confirmButtonCol = "#a565e0")
            r_values$ok_pred = " - "
            r_values$ok_var = " - "
            r_values$weights = 0}
          else {
          
          A = variogramLine(model, dist_vector = dA)
          A = cbind(A, 1)
          A = rbind(A, 1)
          diag(A) = 0
     
          dP = as.vector(st_distance(known_sf, unknown_sf))
          b = c(variogramLine(model, dist_vector = dP)$gamma, 1)
          
          if (det(A) == 0) {shinyalert("Cannot calculate weights", "Change the overlapping points coordinates.",
                                       closeOnClickOutside = T, confirmButtonCol = "#a565e0")
            r_values$ok_pred = " - "
            r_values$ok_var = " - "
            r_values$weights = 0
            }
          else {
            lambda = drop(solve(A) %*% b)
     
            r_values$weights = lambda[1:length(lambda)-1]
          
            r_values$ok_pred = round(as.vector(r_values$weights %*% known$Value),2)
            r_values$ok_var = round(as.vector(t(b) %*% lambda),2)
          }}
          
          # print(r_values$weights)
          # print(ok_pred)
          # print(ok_var)
          # print(dA)
          # print(A)
          # print(dP)
          # print(b)
          # print(lambda)
          # print(r_values$weights)
          
          
      }

      
    }

    #################################################
    

    
    
    
    ###########creating grid (with points coordinates and estimated cell index)#####################
    reactive_grid = eventReactive(input$gen_grid, ignoreNULL = F, {
        
        #create grid
        bb = c(xmin = 0, ymin = 0, xmax = input$width, ymax = input$height)       #bbox data
        class(bb) = "bbox"                                                        #bbox class
        siatka = st_as_stars(bb, deltax = 1, deltay = 1)                          #stars object
        grid = st_make_grid(st_as_sf(siatka), n = c(input$width, input$height))   #grid
        
        
        
        
        #coords
        punkty = st_as_sf(r_values$data, coords = c("X", "Y"))
        #print(punkty)
        
        

        #plot all: stars, grid, known points, estimated cell
        par(bg = "#e6e3ff")
        plot(siatka, main = NULL, reset = FALSE)
        plot(grid, add = TRUE)
        plot(grid[input$ex + (input$ey-1) * input$width], add = TRUE, col = "red")
        plot(st_geometry(punkty), pch = 21, col = "black", bg = "blue",cex = 1.5, add = T)
        
        if (is.na(r_values$chosen[1,1])) {} else {
        plot(st_geometry(st_as_sf(r_values$chosen, coords = c("X", "Y"))), pch = 21, col = "black", bg = "red",cex = 1.5, add = T)}
        
        text(r_values$data$X, r_values$data$Y, labels = rownames(punkty), cex = 1.3, font = 4, col = "white", adj = c(0,0))
        
        if (input$n == "maxdist") {
          draw.circle((input$ex - 0.5), (input$ey - 0.5), input$nval, border = "red", lwd = 1.5, col = rgb(1, 1, 1, 0.1))
        }

        
    })
    
    

    

    ##############################################
    
    
    
    ###########restrictions #1 (grid and estimated cell)  ######
    
    #grid values integer and positive
    observe({
        updateNumericInput(session, "width",
                          value = if (input$width < 4) 4 else round(input$width, 0))
    })
    
    observe({
        updateNumericInput(session, "height",
                           value = if (input$height < 4) 4 else round(input$height, 0))
    })

    #estimated cell positive and within grid
    observe({
      updateNumericInput(session, "ex",
                         value = if (input$ex > input$width) input$width else round(input$ex, 0))
    })
    
    observe({
      updateNumericInput(session, "ey",
                         value = if (input$ey > input$height) input$height else round(input$ey, 0))
    })
    
    
    #change of grid dimensions alters point coordinates
    observeEvent(input$width,
                 {
                   for (i in 1:length(r_values$data$X))
                   {
                     r_values$data$X[i] = if (r_values$data$X[i] > input$width) input$width else r_values$data$X[i]
                   }
                 })
    
    observeEvent(input$height,
                 {
                   for (i in 1:length(r_values$data$Y))
                   {
                     r_values$data$Y[i] = if (r_values$data$Y[i] > input$height) input$height else r_values$data$Y[i]
                   }
                 })
    ####################################
    
    
    
    ##############Create point dataframe############################
    
    #points values
    tx = c(3.5, 5.1, 7.8, 12.0, 11.1, 9.2, 5.1)
    ty = c(8.2, 5.0, 10.7, 8.2, 6.3, 4.3, 12.6)
    tv = c(22.64, 25.1, 33.2, 30.55, 29.7, 26.4, 42.1)
    df = data.frame(X = tx, Y = ty, Value = tv)
    
    r_values = reactiveValues(data = df, chosen = df)
    
    #rHOT table with reactive values
    observe({
      r_values$data = df[1:input$nr,]
      output$table = renderRHandsontable({
        rhandsontable(r_values$data) # converts the R dataframe to rhandsontable object
      })
    })
    
    #po każdym wprowadzeniu zmian r_values$data zostaje zaktualizowana i sprawdzona pod względem ograniczeń wartości
    #r_values$data updates after every change and checks the value restrictions
    observeEvent(
      input$table$changes$changes,
      {
        r_values$data = hot_to_r(input$table)
        
        xi = input$table$changes$changes[[1]][[1]]
        
        #restrictions for X column
        r_values$data[xi+1, 1] = if (is.na(r_values$data[xi+1, 1])) 1 else
                                 if (r_values$data[xi+1, 1] > input$width) input$width else
                                 if (r_values$data[xi+1, 1] <= 0) 1 else
                                   r_values$data[xi+1, 1]
        
        #restrictions for Y column
        r_values$data[xi+1, 2] = if (is.na(r_values$data[xi+1, 2])) 1 else
                                 if (r_values$data[xi+1, 2] > input$height) input$height else
                                 if (r_values$data[xi+1, 2] <= 0) 1 else
                                  r_values$data[xi+1, 2]
        
        #restrictions for Values column
        r_values$data[xi+1, 3] = if (is.na(r_values$data[xi+1, 3])) 1 else
                                  r_values$data[xi+1, 3]


      }     
    )

    
    ###############################################################

    
    ##############restrictions #2 (number of rows and nmax)    ###################
    
    
    #row number integer and between 2 and 7
    observe({
        updateNumericInput(session, "nr",
                           value = if (input$nr < 2) 2 else
                                   if (input$nr > 7) 7 else
                                   if (!is.integer(input$nr)) round(input$nr, 0))
    })
    
    #nmax value not greater than number of rows
    observe({
      updateNumericInput(session, "nval",
                         value = if (is.na(input$nval)) 1 else
                                 if (input$nval > input$nr && input$n == "nmax") input$nr else
                                 if (input$nval <= 0) 1 else
                                 if (!is.integer(input$nval) && input$n == "nmax") round(input$nval, 0)) 
                                 
    })
    
    
    #################################################

    ##########second panel#####################
    
    #semivariogram parameters
    #create fake semivariogram with reactive max gamma
    fake_sv = data.frame(np = 1,
                         dist = 150,
                         gamma = 200,
                         dir.hor = 0,
                         dir.ver = 0,
                         id = "var1")
    class(fake_sv) = c("gstatVariogram", "data.frame")
    
    r_values$fsv = fake_sv
    
    observe({

      r_values$vario = variogram(Value~1, locations = st_as_sf(r_values$data, coords = c("X", "Y")), width = 0.02, cutoff = 8)
      r_values$gamma = max(r_values$vario$gamma)
      
      r_values$fsv$dist = ceiling(sqrt(input$height * input$width) * 0.5)
      r_values$fsv$gamma = round(r_values$gamma*2, -1)+10
      
      if ((input$model == "Power" && input$range > 2))
      {}
      else {
      r_values$preds = variogramLine(vgm(psill = input$psill,
                                         model = substr(input$model, 1, 3),
                                         range = input$range,
                                         nugget = input$nugget),
                                     maxdist = r_values$fsv$dist)
      }
      
      
      # print(r_values$vario)
      # print(r_values$fsv$gamma)
      # print(r_values$fsv$dist)
      # print(max(r_values$preds$gamma))

    })
    
      
    output$plot = renderPlot({ ggplot(r_values$fsv, aes(dist, gamma)) + geom_point(color = "transparent") +
        geom_line(data = r_values$preds, color = "purple") + theme_classic() +
        ylim(0, r_values$fsv$gamma) + xlab("distance") + ylab("semivariance") +
        theme(panel.background = element_rect(fill = "transparent"),
              plot.background = element_rect(fill = "#e6e3ff", color = "#6800a9", size = 2),
              plot.margin = unit(c(7, 7, 10, 10), "pt"),
              axis.line = element_line(size = 0.6),
              axis.text = element_text(size = 10))

      },
                             height = function() { session$clientData$output_plot_width * 0.90 })
    
    
    
    #control of sliders values

    observe({
        updateSliderInput(session, "range",
                          min = if (input$model == "Power") 0.1 else 1,
                          max = if (input$model == "Power") 2 else r_values$fsv$dist + (2/5*r_values$fsv$dist),
                          value = if (input$model == "Power") 1 else 0.5*(r_values$fsv$dist + (2/5*r_values$fsv$dist)),
                          step = if (input$model == "Power") 0.05 else (r_values$fsv$dist + (2/5*r_values$fsv$dist))/100,
                          label = if (input$model == "Spherical") "Range" else "Practical range") })
    observeEvent(r_values$fsv$gamma, {
      updateSliderInput(session, "psill",
                        max = r_values$fsv$gamma*0.5,
                        value = r_values$fsv$gamma*0.25)
      })
    observe({
      updateSliderInput(session, "nugget",
                        max = r_values$fsv$gamma*0.5,
                        value = r_values$fsv$gamma*0.25)
    })
    

    

    
    
    ##########third panel############
    
    #weights plot#
    output$weights = renderPlot( {
      
      rn = if (is.na(r_values$chosen[1,1])) 0 else c(rownames(r_values$chosen))
      
      df = data.frame("point" = rn, "weight" = r_values$weights)
      
      ggplot(df, aes(x=factor(point), y=weight, label = round(weight, 2))) + 
        geom_bar(stat = "identity", fill = "purple", color = "purple") + 
        geom_text(aes(vjust = ifelse(weight > 0, 1.5, -0.5)), color = "white", size = 4.8) +
        theme_classic() + ylab(NULL) + xlab(NULL) + theme(axis.text = element_text(size = 15),
                                                          panel.grid.major.y = element_line(linetype = "dashed"),
                                                          panel.grid = element_line(color = "grey"),
                                                          panel.background = element_rect(fill = "transparent"),
                                                          plot.background = element_rect(fill = "#e6e3ff", color = "#6800a9", size = 2))
      
      },
      
      height = function() { session$clientData$output_plot_width * 0.90 } )
    
    output$estimation = renderText({r_values$ok_pred})
    output$variance = renderText({r_values$ok_var})
    
    observe({
      r_values$est_coords = data.frame(X = input$ex - 0.5, Y = input$ey - 0.5)
      r_values$model = vgm(psill = input$psill, model = substr(input$model, 1, 3), nugget = input$nugget, range = input$range)
      n_function(input$n)
      
      wev(r_values$chosen, r_values$est_coords, r_values$model)
    })
    
    output$siatka = renderPlot( {reactive_grid()}, height = function() { session$clientData$output_plot_width * 0.90 } )
    
}

# Run Krigzy v1.0
shinyApp(ui = ui, server = server)
