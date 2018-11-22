require("plotrix")

euclide <- function(u, v) sqrt(sum((u-v)^2)) # Евклидова метрика

core.E <- function(r) (3/4)*(1-r^2)*(abs(r) <= 1) # Ядро Епанечникова
core.Q <- function(r) (15/16)*((1 - r^2)^2)*(abs(r) <= 1) # Квартическое ядро
core.T <- function(r) (1 - abs(r))*(abs(r) <= 1) # Треугольное ядро
core.G <- function(r) (2*pi)^(-0.5)*exp(-0.5*(r^2)) # Гауссовское ядро
core.P <- function(r) (0.5)*(abs(r) <= 1) # Прямоугольное ядро

getNewSetByPotentials <- function(xl, potentials, res) {
  # Получить новый массив со значениями, соответствующими ненулевым значениям потенциалов
  if(!is.null(ncol(res))) return (res[as.numeric(rownames(xl[which(potentials != 0),])), ]) # Если не вектор
  return (res[as.numeric(rownames(xl[which(potentials != 0),]))]) # Если вектор
}

getDistances <- function(xl, z, metricFunction = euclide) {
  # Посчитать расстояния от объекта z до каждого объекта выборки
  l <- nrow(xl)
  n <- ncol(xl)
  distances <- rep(0, l)
  for (i in 1:l) 
    distances[i] <- metricFunction(xl[i, 1:(n-1)], z)
  return (distances)
}

getHWidth <- function(xl) {
  # Задать ширину окна для каждого объекта выборки
  l <- nrow(xl)
  h <- rep(0, l)
  for(i in 1:l) { 
    h[i] <- 0.4 # ширину окна задаём как 0.4, так это значение является оптимальным
  }
  return (h)
}

getPotentials <- function(xl, h, eps, type_core)
{
  l <- nrow(xl) # строки
  n <- ncol(xl) # столбцы (размерность)
  potentials <- rep(0, l)
  distances_to_points <- matrix(0, l, l)
  err <- eps + 1 # будем считать количество ошибок на выборке 
  # err = (eps + 1) для предотвращения раннего выхода из цикла
  # матрица дистанций до других точек выборки
  for (i in 1:l)
    distances_to_points[i,] <- getDistances(xl, c(xl[i, 1], xl[i, 2])) 
  # xl - выборка, xl[i] - текущая точка и её координата
  while(err > eps){
    while (TRUE) {
      # Продолжаем, пока не встретим ошибочное определение к классу
      cur <- sample(1:l, 1) # выбираем случайную точку из выборки
      class <- potentialFunction(distances_to_points[cur, ], potentials, h, xl, type_core)
      
      if (class != xl[cur, n]) { # если не соответсвует классу
        potentials[cur] = potentials[cur] + 1 # увеличиваем потенциал
        break
      } #if
    } #while(true)
    
    # считаем количество ошибок
    err <- 0
    for (i in 1:l) {
      class <- potentialFunction(distances_to_points[i, ], potentials, h, xl, type_core)
      err <- err + (class != xl[i, n])
    }
  }
  return (potentials)
}

potentialFunction <- function(distances, potentials, h, xl, type_core) {
  l <- nrow(xl)
  n <- ncol(xl)
  classes <- xl[, n]
  weights <- table(classes) # Таблица для весов классов
  weights[1:length(weights)] <- 0 # По умолчанию все веса равны нулю
  for (i in 1:l) { # Для каждого объекта выборки
    class <- xl[i, n] # Берется его класс
    r <- distances[i] / h[i]
    weights[class] <- weights[class] + potentials[i] * type_core(r) # Считается его вес, и прибавляется к общему ввесу его класса
  }
  if (max(weights) != 0) return (names(which.max(weights))) # Если есть веса больше нуля, то вернуть класс с наибольшим весом
  return (0) # Если точка не проклассифицировалась вернуть 0
}

ClassMap <- function(xl, h, potentials, type_core) {
  # Проклассифицируем объекты на основе обучающей выборки, и запишем их в матрицу
  ox <- seq(0, 7, 0.1)
  oy <- seq(0, 2.5, 0.1)
  classMatrix <- matrix(NA, length(ox)*length(oy), ncol(xl))
  cnt <- 1
  for (i in ox)
    for (j in oy) {
      z <- c(i, j)
      distances <- getDistances(xl, z)
      class <- potentialFunction(distances, potentials, h, xl, type_core)
      if (class != 0) {
        classMatrix[cnt, ] <- c(z[1], z[2], class)
        cnt <- cnt + 1
      }
    }
  return (classMatrix)
}

drawPlots <- function(xl, classMatrix, potentials, h) {
  l <- nrow(xl)
  n <- ncol(xl)
  colors <- c("setosa" = "red", "versicolor" = "green3", "virginica" = "blue")
  
  # Полупрозрачные цвета для потенциалов
  transp_red <- col2rgb("red")
  transp_red <- rgb(transp_red[1], transp_red[2], transp_red[3], alpha = 60, max = 255)
  transp_green <- col2rgb("green3")
  transp_green <- rgb(transp_green[1], transp_green[2], transp_green[3], alpha = 60, max = 255)
  transp_blue <- col2rgb("blue")
  transp_blue <- rgb(transp_blue[1], transp_blue[2], transp_blue[3], alpha = 60, max = 255)
  transp_colors <- c("setosa" = transp_red, "versicolor" = transp_green, "virginica" = transp_blue)
  
  par(mfrow=c(1,2))
  # Карта потенциалов
  plot(xl[, 1:(n-1)], pch = 21, bg = colors[xl[,n]], col = colors[xl[,n]], main = "Карта потенциалов", xlab = "Длина лепестка", ylab = "Ширина лепестка", asp = 1)
  for (i in 1:l)
    if (potentials[i] != 0)
      draw.circle(xl[i, 1], xl[i, 2], radius = potentials[i], border = transp_colors[xl[i, n]], col = transp_colors[xl[i, n]])
  # Карта классификации
  plot(xl[, 1:(n-1)], pch = 21, bg = colors[xl[,n]], col = colors[xl[,n]], main = "Метод потенциальных функций", xlab = "Длина лепестка", ylab = "Ширина лепестка", asp = 1)
  points(classMatrix[, 1:(n-1)], pch = 1, col = colors[classMatrix[, n]])
}

main <- function(type_core) {
  xl <- iris[, 3:5]
  h <- getHWidth(xl)
  potentials <- getPotentials(xl, h, 5, type_core)
  
  new_xl <- getNewSetByPotentials(xl, potentials, xl) # новая выборка с теми элементами, у которых потенциалы ненулевые
  new_h <- getNewSetByPotentials(xl, potentials, h) # Соответствующий вектор с ширинами окон этих элементов
  new_potentials <- getNewSetByPotentials(xl, potentials, potentials) # Соответствующие ненулевые потенциалы
  
  classMatrix <- ClassMap(new_xl, new_h, new_potentials, type_core)
  drawPlots(xl, classMatrix, potentials, h)
}

main(core.E)