library(shiny)

ui <- fluidPage(
  
  titlePanel("Линейные классификаторы"),
  
  sidebarLayout(
    sidebarPanel(
      fluidRow(
        column(12, radioButtons("classifiers", "Классификаторы", c("ADALINE" = 0, 
                                                                  "Персептрон Розенблатта" = 1, 
                                                                  "Логистич. регрессия" = 2, 
                                                                  "Сравнение" = 3),
                                selected = 0, inline = T), style = "text-align: center"),
        column(12, "Дополнительные свойства:", style = "text-align: center; font-size: 16px"),
        column(6, checkboxInput("drow_iter", "Отображать итерации", F)), 
        column(6, checkboxInput("drow_q", "Отображать график Q", T)),
        column(6, selectInput("eps_q", "Погрешность Q", c(0.00001, 0.0001, 0.001, 0.01, 0.1), 0.001)),
        column(6, selectInput("max_cnt", "MAX количество итераций", c(1000, 2000, 5000, 10000, 20000), 5000)),
        
        column(12, "Особые критерии остановки (не учитывается Q):", style = "text-align: center; font-size: 16px"),
        column(6, checkboxInput("adaline_stop", "Для Adaline", F)), 
        column(6, checkboxInput("logistic_stop", "Для логистич. регрессии", F)),
        column(12, sliderInput("n", "Количество объектов", 100, 600, 200, 1)),
        
        column(12, "Класс 1", style = "color: red; text-align: center; font-size: 24px"),
        column(6, sliderInput("mu11", "Мат. ожид. X1", 0, 20, 0, 1)), 
        column(6, sliderInput("mu21", "Мат. ожид. X2", 0, 20, 4, 1)),
        column(6, sliderInput("sigma11", "Отклонение X1", 0.1, 1, 0.5, 0.1)), 
        column(6, sliderInput("sigma21", "Отклонение X2", 0.1, 1, 0.5, 0.1)),
        
        column(12, "Класс 2", style = "color: blue; text-align: center; font-size: 24px"),
        column(6, sliderInput("mu12", "Мат. ожид. X1", 0, 20, 4, 1)), 
        column(6, sliderInput("mu22", "Мат. ожид. X2", 0, 20, 0, 1)),
        column(6, sliderInput("sigma12", "Отклонение X1", 0.1, 1, 0.5, 0.1)), 
        column(6, sliderInput("sigma22", "Отклонение X2", 0.1, 1, 0.5, 0.1)),
        
        column(12, uiOutput("ui"))
      )
    ),
    
    mainPanel(
      plotOutput("plot", width = "1000px", height = "600px")
    )
  )
)

normalize_sample <- function(xl) {
  for (i in 1:(ncol(xl)-1)) xl[,i] <- (xl[,i] - mean(xl[,i])) / sd(xl[,i])
  return (xl)
}

margin <- function(w, object, class) w %*% object * class

loss.Quad <- function(m) (1-m)^2
loss.Pers <- function(m) max(-m,0)
loss.Log <- function(m) log2(1 + exp(-m))
sigmoid <- function(z) 1 / (1 + exp(-z))

adaline.get_w <- function(w, object, class, eta) w - c(eta) * (w %*% object - class) %*% object
hebb.get_w <- function(w, object, class, eta) w + c(eta) * object * class
logistic.get_w <- function(w, object, class, eta) w + c(eta) * object * class * sigmoid(c(w %*% object) * class)

gradient <- function(xl, eta, lambda, rule, loss_function, max_cnt, eps_q, stop_criterion) { # eta - темп обучения, lambda - параметр сглаживания
  l <- nrow(xl)
  n <- ncol(xl)
  w <- matrix(c(runif(n-1, -1/(2*(n-1)), 1/(2*(n-1)))), 1, 3) # генерация w случайными весами
  objects <- xl[,-n]
  classes <- xl[, n]
  q <- sum(sapply(1:l, function(i) loss_function(margin(w, objects[i,], classes[i])))) # инициализация q
  q_full <- matrix(q, 1, 1)
  
  cnt <- 0
  while (T) {
    cnt <- cnt + 1 # Счётчик итераций
    margins <- sapply(1:l, function(i) margin(w[cnt,], objects[i,], classes[i]))
    errors <- which(margins < 0)
    
    if (length(errors) == 0 && stop_criterion == T) break;
    
    if (length(errors) > 0) rand <- sample(errors, 1)
    else rand <- sample(1:l, 1)
    
    eps <- loss_function(margin(w[cnt,], objects[rand,], classes[rand])) # Ошибка алгоритма на объекте
    
    eta <- 1 / (objects[rand,] %*% objects[rand,])^2 # изменение темпа обучения
    
    # Обновление весов
    w <- rbind(w, rule(w[cnt,], objects[rand,], classes[rand], eta))
    
    # Перерасчёт q
    q_prev <- q
    q <- (1 - lambda) * q + lambda * eps
    q_full <- rbind(q_full, q)
    
    if (abs(q_prev - q) / max(q_prev, q) <= eps_q) # выход из цикла
      { break; }
    else if (cnt == max_cnt) 
      { break; }
    
  }
  w <- cbind(w, q_full)
  return (w)
}

adaline <- function(xl, drow_iter, drow_q, max_cnt, eps_q, stop_criterion) {
  w <- gradient(xl, 1, 1/6, adaline.get_w, loss.Quad, max_cnt, eps_q, stop_criterion)
  q <- w[,ncol(w)]
  n <- ncol(xl)
  l <- nrow(w)
  if (drow_q == T) {
    par(mfrow=c(1, 2))
    plot(q, type = "l", bg = "red", col = "red", main = "График изменения Q", xlab = "Итерации", ylab = "Значения Q")
  }
  colors <- c("blue", "white", "red")
  plot(xl[,1:(n-2)], type="n", asp = 1, main = "Классификатор: ADALINE")
  if (drow_iter == T) 
    for (i in 1:(l-1)) 
      abline(a = w[i,3]/w[i,2], b = -w[i,1]/w[i,2], lwd = 1, col = "black")
  abline(a = w[l,3]/w[l,2], b = -w[l,1]/w[l,2], lwd = 3, col = "green")
  points(xl[,1:(n-2)], pch = 21, col = colors[xl[,n]+2], bg = colors[xl[,n]+2])
}

perceptron <- function(xl, drow_iter, drow_q, max_cnt, eps_q, stop_criterion) {
  w <- gradient(xl, 1, 1/6, hebb.get_w, loss.Pers, max_cnt, eps_q, stop_criterion)
  q <- w[,ncol(w)]
  n <- ncol(xl)
  l <- nrow(w)
  if (drow_q == T) {
    par(mfrow=c(1, 2))
    plot(q, type = "l", bg = "red", col = "red", main = "График изменения Q", xlab = "Итерации", ylab = "Значения Q")
  }
  colors <- c("blue", "white", "red")
  plot(xl[,1:(n-2)], type="n", asp = 1, main = "Классификатор: Персептрон Розенблатта")
  if (drow_iter == T) 
    for (i in 1:(l-1)) 
      abline(a = w[i,3]/w[i,2], b = -w[i,1]/w[i,2], lwd = 1, col = "black")
  abline(a = w[l,3]/w[l,2], b = -w[l,1]/w[l,2], lwd = 3, col = "gold")
  points(xl[,1:(n-2)], pch = 21, col = colors[xl[,n]+2], bg = colors[xl[,n]+2])
}

logistic <- function(xl, drow_iter, drow_q, max_cnt, eps_q, stop_criterion, map) {
  w <- gradient(xl, 1, 1/6, logistic.get_w, loss.Log, max_cnt, eps_q, stop_criterion)
  q <- w[,ncol(w)]
  n <- ncol(xl)
  l <- nrow(w)
  if (drow_q == T) {
    par(mfrow=c(1, 2))
    plot(q, type = "l", bg = "red", col = "red", main = "График изменения Q", xlab = "Итерации", ylab = "Значения Q")
  }
  colors <- c("blue", "white", "red")
  plot(xl[,1:(n-2)], type="n", asp = 1, main = "Классификатор: логистическая регрессия")
  if (drow_iter == T) 
    for (i in 1:(l-1)) 
      abline(a = w[i,3]/w[i,2], b = -w[i,1]/w[i,2], lwd = 1, col = "black")
  abline(a = w[l,3]/w[l,2], b = -w[l,1]/w[l,2], lwd = 3, col = "violet")
  points(xl[,1:(n-2)], pch = 21, col = colors[xl[,n]+2], bg = colors[xl[,n]+2])
}

compare <- function(xl, drow_iter, max_cnt, eps_q, stop_criterion) {
  n <- ncol(xl)
  w1 <- gradient(xl, 1, 1/6, adaline.get_w, loss.Quad, max_cnt, eps_q, stop_criterion["Adaline"])
  w2 <- gradient(xl, 1, 1/6, hebb.get_w, loss.Pers, max_cnt, eps_q, stop_criterion["Hebb"])
  w3 <- gradient(xl, 1, 1/6, logistic.get_w, loss.Log, max_cnt, eps_q, stop_criterion["Logistic"])
  l1 <- nrow(w1)
  l2 <- nrow(w2)
  l3 <- nrow(w3)
  colors <- c("blue", "white", "red")
  plot(xl[,1:(n-2)], type="n", asp = 1, main = "Сравнение линейных классификаторов")
  if (drow_iter == T) {
    for (i in 1:(l1-1)) 
      abline(a = w1[i,3]/w1[i,2], b = -w1[i,1]/w1[i,2], lwd = 1, col = "black")
    for (i in 1:(l2-1)) 
      abline(a = w2[i,3]/w2[i,2], b = -w2[i,1]/w2[i,2], lwd = 1, col = "black")
    for (i in 1:(l3-1)) 
      abline(a = w3[i,3]/w3[i,2], b = -w3[i,1]/w3[i,2], lwd = 1, col = "black")
  }
  abline(a = w1[l1,3]/w1[l1,2], b = -w1[l1,1]/w1[l1,2], lwd = 3, col = "green")
  abline(a = w2[l2,3]/w2[l2,2], b = -w2[l2,1]/w2[l2,2], lwd = 3, col = "gold")
  abline(a = w3[l3,3]/w3[l3,2], b = -w3[l3,1]/w3[l3,2], lwd = 3, col = "violet")
  points(xl[,1:(n-2)], pch = 21, col = colors[xl[,n]+2], bg = colors[xl[,n]+2])
  legend("bottomright", c("ADALINE", "Персептрон", "Логистич. регрессия"), pch = c("l","l","l"), col = c("green", "gold", "violet"))
}

server <- function(input, output) {
  
  output$plot = renderPlot({
    n <- input$n
    xl11 <- rnorm(n/2, input$mu11, input$sigma11)
    xl21 <- rnorm(n/2, input$mu21, input$sigma21)
    xl12 <- rnorm(n/2, input$mu12, input$sigma12)
    xl22 <- rnorm(n/2, input$mu22, input$sigma22)
    tmp1 <- cbind(xl11, xl21)
    tmp2 <- cbind(xl12, xl22)
    xl <- rbind(cbind(tmp1, 1), cbind(tmp2, -1)) # первый класс с 1, второй с -1
    colnames(xl) <- c("X1", "X2", "class")
    
    norm_xl <- normalize_sample(xl)
    norm_xl <- cbind(norm_xl[,1:(ncol(xl)-1)], -1, norm_xl[, ncol(norm_xl)])
    
    drow_iter <- input$drow_iter
    drow_q <- input$drow_q
    
    stop_criterion <- c("Adaline" = input$adaline_stop, "Hebb" = T, "Logistic" = input$logistic_stop)
    eps_q <- as.numeric(input$eps_q)
    max_cnt <- as.numeric(input$max_cnt)
    
    if (input$classifiers == 0) 
      adaline(norm_xl, drow_iter, drow_q, max_cnt, eps_q, stop_criterion["Adaline"])
    else if (input$classifiers == 1) 
      perceptron(norm_xl, drow_iter, drow_q, max_cnt,eps_q, stop_criterion["Hebb"])
    else if (input$classifiers == 2) 
      logistic(norm_xl, drow_iter, drow_q, max_cnt, eps_q, stop_criterion["Logistic"])
    else if (input$classifiers == 3) 
      compare(norm_xl, drow_iter, max_cnt, eps_q, stop_criterion)
  })
  
}

shinyApp(ui, server)