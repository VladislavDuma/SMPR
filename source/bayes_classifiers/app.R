library(shiny)

ui <- fluidPage(
  
  titlePanel("Байесовский классификатор"),
  titlePanel("ВОВ"),
  
  sidebarLayout(
    sidebarPanel(
      fluidRow(
        column(12, radioButtons("classifiers", "Классификатор", 
                                c("Plug-in" = 1, "ЛДФ" = 2), selected = 1, inline = T)),
        column(12, "СССР", style = "color: red; text-align: center; font-size: 24px"),
        column(12, sliderInput("s", "Количество отрядов", 10, 300, 20, 1)),
        column(12, "Центр штаба", style = "color: red; text-align: center; font-size: 12px"),
        column(6, sliderInput("mu11", "Центр по OX", 0, 6, 2, 0.1)), 
        column(6, sliderInput("mu21", "Центр по OY", 0, 4, 1, 0.1)),
        column(12, "Распределение сил от центра", style = "color: red; text-align: center; font-size: 12px"),
        column(6, sliderInput("sigma11", "Среднее отдаление от штаба OX", 0.1, 1, 0.5, 0.1)), 
        column(6, sliderInput("sigma21", "Среднее отдаление от штаба OY", 0.1, 1, 1, 0.1)),
        
        column(12, sliderInput("P", "Вероятность появления отряда СССР|Germany", 0.01, 0.99, 0.5, 0.01)),
        
        column(12, "Germany", style = "color: blue; text-align: center; font-size: 24px"),
        column(12, sliderInput("g", "Количество отрядов", 10, 300, 20, 1)),
        column(12, "Центр штаба", style = "color: blue; text-align: center; font-size: 12px"),
        column(6, sliderInput("mu12", "Центр по OX", 0, 6, 4, 0.1)), 
        column(6, sliderInput("mu22", "Центр по OY", 0, 4, 3, 0.1)),
        column(12, "Распределение сил от центра", style = "color: blue; text-align: center; font-size: 12px"),
        column(6, sliderInput("sigma12", "Среднее отдаление от штаба OX", 0.1, 1, 0.5, 0.1)), 
        column(6, sliderInput("sigma22", "Среднее отдаление от штаба OY", 0.1, 1, 1, 0.1))
      )
    ),
    
    mainPanel(
      plotOutput("plot", width = "800px", height = "600px"),
      #plotOutput("plot1", width = "800px", height = "600px"),
      textOutput("error")
    )
  )
)

# PLUG_IN
plug_in <- function(units_array, sssr, germ, Py) {
  
  get_hq <- function(units_array) colMeans(units_array)
  
  get_sigma <- function(units_array, hq_center) {
    sum <- 0
    for (i in 1:nrow(units_array)) {
      xi <- matrix(c(units_array[i,1], units_array[i,2]), 1, 2)
      sum <- sum + t(xi - hq_center) %*% (xi - hq_center)
    }
    sum / (nrow(units_array)-1)
  }
  
  get_discriminant_coeffs <- function(hq_center, sigma, Py) {
    hq_center1 <- matrix(c(hq_center[1,1], hq_center[1,2]),1,2)
    hq_center2 <- matrix(c(hq_center[2,1], hq_center[2,2]),1,2)
    sigma1 <- sigma[1:2,]
    sigma2 <- sigma[3:4,]
    revSigma1 <- solve(sigma1)
    revSigma2 <- solve(sigma2)
    a <- revSigma1[1,1]
    b <- revSigma1[1,2]
    c <- revSigma1[2,2]
    k <- revSigma2[1,1]
    l <- revSigma2[1,2]
    n <- revSigma2[2,2]
    
    # q^2
    x2 <- a - k
    y2 <- c - n
    xy <- 2 * b - 2 * l
    # q^1
    x <- 2 * l * hq_center2[2] - 2 * b * hq_center1[2] - 2 * a * hq_center1[1] + 2 * k * hq_center2[1]
    y <- 2 * l * hq_center2[1] + 2 * n * hq_center2[2] - 2 * b * hq_center1[1] - 2 * c * hq_center1[2]
    # c
    f <- -k * hq_center2[1]^2 - 2 * l * hq_center2[1] * hq_center2[2] - n * hq_center2[2]^2 + a * hq_center1[1]^2 +
      2 * b * hq_center1[1] * hq_center1[2] + c * hq_center1[2]^2 + log(det(sigma1)) - log(det(sigma2)) - Py[1]/Py[2]
    #print(c("x^2" = x2, "y^2" = y2, "xy" = xy, "x" = x, "y" = y, "1" = f))
    return (c("x^2" = x2, "y^2" = y2, "xy" = xy, "x" = x, "y" = y, "1" = f))
  }
  
  draw_force_balance_line <- function(units_array, hq_center, sigma, Py, units_dispersion) {
    x <- seq(units_dispersion[1,1] - 5, units_dispersion[1,2] + 5, length.out = 100)
    y <- seq(units_dispersion[2,1] - 5, units_dispersion[2,1] + 5, length.out = 100)
    
    coeffs <- get_discriminant_coeffs(hq_center, sigma, Py)
    #print(coeffs)
    z <- outer(x, y, function(x, y) coeffs["x^2"]*x^2 + coeffs["xy"]*x*y + 
                 coeffs["y^2"]*y^2 + coeffs["x"]*x + coeffs["y"]*y + coeffs["1"])
    n <- ncol(units_array)
    colors <- c("germ" = "blue", "sssr" = "red")
    plot(units_array[,1:(n-1)], 
         pch = 21, 
         bg = colors[units_array[,n]], 
         col = colors[units_array[,n]], 
         main = "Карта классификации нормального распределения\n Зелёная линия - баланс сил между конфликтующими сторонами", 
         asp = 1)
    contour(x, y, z, 
            lwd = 4, 
            col = "green", 
            levels = 0, 
            drawlabels = F,
            add = T)
  }
  
  main <- function(units_array, sssr, germ, Py) {
    s <- sssr
    g <- germ
    sum <- s + g
    sssr_m <- units_array[1:s,1:2]
    germ_m <- units_array[(s+1):sum,1:2]
    
    hq_center <- rbind(get_hq(sssr_m), get_hq(germ_m))
    hq_center <- matrix(c(hq_center[1,1], 
                          hq_center[2,1], 
                          hq_center[1,2], 
                          hq_center[2,2]), 
                        2, 2)
    
    #print(hq_center)
    
    sigma <- rbind(get_sigma(sssr_m, hq_center[1,]), 
                   get_sigma(germ_m, hq_center[2,]))
    
    units_dispersion <- matrix(c(min(hq_center[,1]), 
                                 min(hq_center[,2]), 
                                 max(hq_center[,1]), 
                                 max(hq_center[,2])), 
                               2, 2)
    
    classes <- unique(units_array[,ncol(units_array)])
    
    draw_force_balance_line(units_array, hq_center, sigma, Py, units_dispersion)
  }
  
  main(units_array, sssr, germ, Py)
}

#LDF
LDF <- function(units_array, s, g, Py) {
  
  get_hq <- function(units_array) colMeans(units_array)
  
  get_sigma <- function(units_array, hq_center, s, g) {
    sum <- 0
    #print(hq_center)
    m <- nrow(units_array)
    for (i in 1:s) {
      xi <- matrix(c(units_array[i,1], units_array[i,2]), 1, 2)
      sum <- sum + t(xi - hq_center[1,]) %*% (xi - hq_center[1,])
    }
    for (i in (s+1):(s+g)) {
      xi <- matrix(c(units_array[i,1], units_array[i,2]), 1, 2)
      sum <- sum + t(xi - hq_center[2,]) %*% (xi - hq_center[2,])
    }
    sum / (m-2)
  }
  
  get_discriminant_coeffs <- function(hq_center, sigma) {
    hq_center1 <- matrix(c(hq_center[1,1], hq_center[1,2]),1,2)
    hq_center2 <- matrix(c(hq_center[2,1], hq_center[2,2]),1,2)
    revSigma <- solve(sigma)
    a <- revSigma[1,1]
    b <- revSigma[1,2]
    c <- revSigma[2,2]
    x <- 2 * b * hq_center2[2] - 2 * b * hq_center1[2] - 2 * a * hq_center1[1] + 2 * a * hq_center2[1]
    y <- 2 * b * hq_center2[1] + 2 * c * hq_center2[2] - 2 * b * hq_center1[1] - 2 * c * hq_center1[2]
    f <- -a * hq_center2[1]^2 - 2 * b * hq_center2[1] * hq_center2[2] - c * hq_center2[2]^2 + 
      a * hq_center1[1]^2 + 2 * b * hq_center1[1] * hq_center1[2] + c * hq_center1[2]^2 - 2 * Py[1]/Py[2]
    return (c("x" = x, "y" = y, "1" = f))
  }
  
  draw_force_balance_line <- function(units_array, hq_center, sigma, units_dispersion) {
    x <- seq(units_dispersion[1,1] - 5, units_dispersion[1,2] + 5, length.out = 100)
    y <- seq(units_dispersion[2,1] - 5, units_dispersion[2,1] + 5, length.out = 100)
    coeffs <- get_discriminant_coeffs(hq_center, sigma)
    
    z <- outer(x, y, function(x, y) coeffs["x"]*x + coeffs["y"]*y + coeffs["1"])
    n <- ncol(units_array)
    
    colors <- c("germ" = "blue", "sssr" = "red")
    plot(units_array[,1:(n-1)], 
         pch = 21, 
         bg = colors[units_array[,n]], 
         col = colors[units_array[,n]], 
         main = "Карта классификации нормального распределения\n Зелёная линия - баланс сил между конфликтующими сторонами", 
         asp = 1)
    contour(x, y, z, 
            lwd = 3, 
            levels = 0, 
            col = "green", 
            drawlabels = F,
            add = T)
  }
  
  main <- function(units_array, sssr, germ, Py) {
    s <- sssr
    g <- germ
    sum <- s + g
    sssr_m <- units_array[1:s,1:2]
    germ_m <- units_array[(s+1):sum,1:2]
    
    hq_center <- rbind(get_hq(sssr_m), get_hq(germ_m))
    print(hq_center)
    hq_center <- matrix(c(hq_center[1,1], hq_center[2,1], hq_center[1,2], hq_center[2,2]), 2, 2)
    print(hq_center)
    
    sigma <- (get_sigma(units_array, hq_center, s, g))
    
    units_dispersion <- matrix(c(min(hq_center[,1]), 
                                 min(hq_center[,2]), 
                                 max(hq_center[,1]), 
                                 max(hq_center[,2])), 
                               2, 2)
    
    classes <- unique(units_array[,ncol(units_array)])
    
    draw_force_balance_line(units_array, hq_center, sigma, units_dispersion)
  }
  
  main(units_array, s, g, Py)
}

# server
server <- function(input, output) {
  
  output$plot <- renderPlot({
      if(input$classifiers == 1) {
      s <- input$s #sssr
      g <- input$g #germ
      Ps <- input$P
      Pg <- 1 - Ps
      Py <- c(Ps, Pg)
      
      sum <- s + g
      
      sigma1 <- matrix(c(input$sigma11, 0, 0, input$sigma21), 2, 2)
      sigma2 <- matrix(c(input$sigma12, 0, 0, input$sigma22), 2, 2)
      
      hq_center1 <- c(input$mu11, input$mu21)
      hq_center2 <- c(input$mu12, input$mu22)
      
      sssr <- mvrnorm(s, hq_center1, sigma1)
      germ <- mvrnorm(g, hq_center2, sigma2)
      
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
      
      plug_in(units_array, s, g, Py)
    }
    else if (input$classifiers == 2) {
      s <- input$s #sssr
      g <- input$g #germ
      Ps <- input$P
      Pg <- 1 - Ps
      Py <- c(Ps, Pg)
      
      sum <- s + g
      
      sigma1 <- matrix(c(input$sigma11, 0, 0, input$sigma21), 2, 2)
      sigma2 <- matrix(c(input$sigma12, 0, 0, input$sigma22), 2, 2)
      if (sigma1[1,1] != sigma2[1,1] || sigma1[2,2] != sigma2[2,2]) {
        output$error = renderText("Отклонения соответствующих признаков Sigma должны быть равны")
        return ()
      }
      hq_center1 <- c(input$mu11, input$mu21)
      hq_center2 <- c(input$mu12, input$mu22)
      
      sssr <- mvrnorm(s, hq_center1, sigma1)
      germ <- mvrnorm(g, hq_center2, sigma2)
      
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
      print(nrow(units_array))
      
      LDF(units_array, s, g, Py)
    }
  })
  #output$plot1 <- renderPlot
}

shinyApp(ui = ui, server = server)