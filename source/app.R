#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
   # Application title
   titlePanel("Level lines"),
   
   # Sidebar with a slider and radio input 
   sidebarLayout(
      sidebarPanel(
         fluidRow(
           column(6, radioButtons("lines", "Lines: ", c("Off" = 0,"On" = 1), selected = 1, inline = TRUE),
                  style = "text-align: center"),
           column(6, radioButtons("map", "Maps: ", c("Off" = 0,"On" = 1), selected = 1, inline = TRUE),
                  style = "text-align: center")
         ),
         fluidRow(
           column(6, sliderInput("a11", "A11:", min = 0, max = 10, value = 1, step = 0.1)),
           column(6, sliderInput("a22", "A22:", min = 0, max = 10, value = 1, step = 0.1))
         ),
         sliderInput("a12", "A12, A21:", min = -10, max = 10, value = 0, step = 0.1),
         fluidRow(
           column(10, NULL),
           column(12, textOutput("message"), style = "color: red; text-align: center")
         ),
         
         sliderInput("step", "Step:", min = 0.001, max = 0.02, value = 0.005, step = 0.001),
         fluidRow(
           column(6, sliderInput("mX", "mX:", min = -3, max = 3, value = 0, step = 0.1)),
           column(6, sliderInput("mY", "mY:", min = -3, max = 3, value = 0, step = 0.1))
         )
      ),
      
      # Show a plot of the generated distribution
      mainPanel(
         plotOutput("distPlot", width = "800px", height = "600px")
      )
   )
)

point <- function(x, mu, sigma) { 
  k <- 1 / sqrt((2 * pi) ^ 2 * det(sigma))
  e <- exp(-0.5 * (x - mu) %*% solve(sigma) %*% t(x - mu))
  return (k * e)
}

# Define server logic required to draw a histogram
server <- function(input, output) {
   
   output$distPlot <- renderPlot({
      # generate bins based on input$bins from ui.R c(as.numeric(input$A11))
      m <- matrix(c(input$mX, input$mY), 1, 2)
      sigma <- matrix(c(input$a11, input$a12, input$a12, input$a22), 2, 2)
      
      # check discrim on <= 0 
      if (det(sigma) <= 0) {
        output$message <- renderText({"det < 0"})
        return()
      }
      output$message = renderText({"det > 0"})
      
      #par(bg = 'white', fg = 'black')
      plot(-5:5, -5:5, type = "n", asp = 1)
      
      minX <- -sigma[1, 1] - 4
      maxX <- sigma[1, 1] + 4
      minY <- -sigma[2, 2] - 3
      maxY <- sigma[2, 2] + 3
      
      x <- seq(minX, maxX, len = 100)
      y <- seq(minY, maxY, len = 100)
      
      z <- matrix(NA, length(x), length(y))
      
      for(i in 1:length(x)) {
        for(j in 1:length(y)) {
          A <- matrix(c(x[i], y[j]), 1, 2)
          z[i,j] <- point(A, m, sigma)
          
          
        }
      }
      
      add = F
      if(input$lines == 1) {
        for (level in seq(0, 0.2, input$step)) 
        {
          contour(x, y, z, levels = level, drawlabels = T, col = "black", add = add, asp = 1) 
          add = T
        }
      }
      
      if(input$map == 1) {
        for(i in 1:length(x)) {
          for(j in 1:length(y)) {
            p <- c(x[i], y[j])
            color <- adjustcolor("red", z[i, j]*10)
            points(p[1], p[2], pch = 16, col = color, bg = color)
          }
        }
      }
   })
}

# Run the application 
shinyApp(ui = ui, server = server)

