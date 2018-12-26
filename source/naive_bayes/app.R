library(shiny)

ui <- fluidPage(
  
  titlePanel("Наивный байесовский классификатор"),
  titlePanel("ВОВ"),
  
  sidebarLayout(
    sidebarPanel(
      fluidRow(
        column(12, "СССР", style = "color: red; text-align: center; font-size: 24px"),
        column(6, sliderInput("n", "Количество отрядов", 10, 100, 20, 1)),
        column(6, sliderInput("lmbd1", "Сила отрядов", 0.5, 2, 1, 0.01)),
        column(12, "Центр штаба", style = "color: red; text-align: center; font-size: 12px"),
        column(6, sliderInput("mu11", "Центр по OX", 0, 6, 2, 0.1)), 
        column(6, sliderInput("mu21", "Центр по OY", 0, 4, 1, 0.1)),
        column(12, "Распределение сил от центра", style = "color: red; text-align: center; font-size: 12px"),
        column(6, sliderInput("sigma11", "Среднее отдаление от штаба OX", 0.1, 1, 0.7, 0.1)), 
        column(6, sliderInput("sigma21", "Среднее отдаление от штаба OY", 0.1, 1, 1, 0.1)),
        
        column(12, sliderInput("P", "Вероятность появления отряда СССР|Germany", 0.01, 0.99, 0.5, 0.01)),
        
        column(12, "Germany", style = "color: blue; text-align: center; font-size: 24px"),
        column(6, sliderInput("g", "Количество отрядов", 10, 100, 20, 1)),
        column(6, sliderInput("lmbd2", "Сила отрядов", 0.5, 2, 1, 0.01)),
        column(12, "Центр штаба", style = "color: blue; text-align: center; font-size: 12px"),
        column(6, sliderInput("mu12", "Центр по OX", 0, 6, 4, 0.1)), 
        column(6, sliderInput("mu22", "Центр по OY", 0, 4, 3, 0.1)),
        column(12, "Распределение сил от центра", style = "color: blue; text-align: center; font-size: 12px"),
        column(6, sliderInput("sigma12", "Среднее отдаление от штаба OX", 0.1, 1, 0.6, 0.1)), 
        column(6, sliderInput("sigma22", "Среднее отдаление от штаба OY", 0.1, 1, 0.7, 0.1))
      )
    ),
    
    mainPanel(
      plotOutput("plot", width = "800px", height = "600px")
    )
  )
)

naive_bayes <- function(units_array, sssr, germ, P, l) {
  p <- function(ksi, mu, sigma) (1 / (sigma * sqrt(2 * pi))) * exp(-(ksi - mu) ^ 2 / (2 * sigma ^ 2)) # density
  
  classificator <- function(x, classes, mu, sigma, Py, lambda) { # lambda - class weight (potential)
    sum_class <- rep(0, length(classes))
    names(sum_class) <- classes
    
    for (i in 1:length(sum_class)) {
      sum <- 0
      for (j in 1:length(x)) sum <- sum + log(p(x[j], mu[i,j], sigma[i,j]))
      sum_class[i] <- log(lambda[i] * Py[i]) + sum
    }
    return (names(which.max(sum_class)))
  }
  
  strategic_map <- function(classes, mu, sigma, Py, units_dispersion, l) {
    classified_map <- c()
    for (i in seq(units_dispersion[1,1] - 6, units_dispersion[1,2] + 6, 0.1))
      for (j in seq(units_dispersion[2,1] - 6, units_dispersion[2,1] + 6, 0.1)) 
        classified_map <- rbind(classified_map, 
                                c(i, j, classificator(c(i, j), classes, mu, sigma, Py, l)))
    #print(classified_map)
    return (classified_map)
  }
  
  draw_strategic_map <- function(units_array, classified_territory) {
    n <- ncol(units_array)
    colors <- c("germ" = "blue", "sssr" = "red")
    
    plot(units_array[,1:(n-1)], 
         pch = 21, 
         bg = colors[units_array[,n]], 
         col = colors[units_array[,n]], 
         main = "Карта распределения сил", asp = 1)
    points(classified_territory[,1:(n-1)], 
           pch = 21, 
           col = colors[classified_territory[,n]])
  }
  
  get_hq <- function(units_array) sum(units_array) / length(units_array)
  
  get_sigma <- function(units_array, mu) 
    sum((units_array - mu)^2) / (length(units_array)-1)
  
  # main function for naive bayes algorithm
  main <- function(units_array, sssr, germ, P, l) {
    Py <- P
    s <- sssr
    g <- germ
    sum <- s + g
    sssr_x <- units_array[1:s,1]
    sssr_y <- units_array[1:s,2]
    germ_x <- units_array[(s+1):sum,1]
    germ_y <- units_array[(s+1):sum,2]
    #print(sssr_x)
    #print(sssr_y)
    
    hq_center <- rbind(c(get_hq(sssr_x), 
                         get_hq(sssr_y)), 
                       c(get_hq(germ_x), 
                         get_hq(germ_y))) # search center
    
    sigma <- rbind(c(get_sigma(sssr_x, hq_center[1,1]), 
                     get_sigma(sssr_y, hq_center[1,2])), 
                   c(get_sigma(germ_x, hq_center[2,1]), 
                     get_sigma(germ_y, hq_center[2,2]))) # search sigma
    
    classes <- unique(units_array[,ncol(units_array)])
    
    units_dispersion <- matrix(c(min(hq_center[,1]), 
                                 min(hq_center[,2]), 
                                 max(hq_center[,1]), 
                                 max(hq_center[,2])), 
                               2, 2)
    #print(units_dispersion)
    
    classified_territory <- strategic_map(classes, hq_center, sigma, Py, units_dispersion, l)
    draw_strategic_map(units_array, classified_territory)
  }
  
  main(units_array, sssr, germ, P, l)
}

# server
server <- function(input, output) {
  
  output$plot <- renderPlot({
    s <- input$n #sssr
    g <- input$g #germ
    Ps <- input$P
    Pg <- 1 - Ps
    Py <- c(Ps, Pg)
    #print(Py)
    l1 <- input$lmbd1
    l2 <- input$lmbd2
    l <- c(l1, l2)
    #print(P)
    #print(l)
    sum <- s + g
    
    sssr_x <- rnorm(s, input$mu11, input$sigma11)
    sssr_y <- rnorm(s, input$mu21, input$sigma21)
    germ_x <- rnorm(g, input$mu12, input$sigma12)
    germ_y <- rnorm(g, input$mu22, input$sigma22)
    
    sssr <- cbind(sssr_x, sssr_y)
    germ <- cbind(germ_x, germ_y)
    
    colnames(sssr) <- c()
    colnames(germ) <- c()
    
    units_array <- data.frame()
    units_array <- rbind(units_array, sssr)
    units_array <- rbind(units_array, germ)
    
    classes <- 1:sum
    classes[1:s] <- "sssr"
    classes[(s+1):sum] <- "germ"
    units_array <- cbind(units_array, classes)
    
    colnames(units_array) <- c("X координаты", "Y координаты", "Сторона конфликта")
    naive_bayes(units_array, s, g, Py, l)
  })
}

shinyApp(ui = ui, server = server)