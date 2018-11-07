euclide <- function(u, v) sqrt(sum((u-v)^2)) # метрика

core.E <- function(r) (3/4)*(1-r^2)*(abs(r) <= 1) # Ядро Епанечникова
core.Q <- function(r) (15/16)*((1 - r^2)^2)*(abs(r) <= 1) # Квартическое ядро
core.T <- function(r) (1 - abs(r))*(abs(r) <= 1) # Треугольное ядро
core.G <- function(r) (2*pi)^(-0.5)*exp(-0.5*(r^2)) # Гауссовское ядро
core.P <- function(r) (0.5)*(abs(r) <= 1) # Прямоугольное ядро


Distances <- function(xl, z, metricFunction = euclide) {
  # Получить вектор расстояний от объекта z до каждого объекта выборки
  l <- nrow(xl)
  n <- ncol(xl)
  distances <- rep(0, l)
  for (i in 1:l)
    distances[i] <- metricFunction(xl[i, 1:(n-1)], z)
  return (distances)
}

parzen <- function(xl, h, distances, type_core) {
  # Оценка весовой функции по расстоянию, а не по рангу 
  # h - ширина окна
  # distances - расстояния от точки z до каждого объекта из выборки xl 
  # type_core - текущая функция ядра
  l <- nrow(xl) # строки
  n <- ncol(xl) # столбцы (размерность)
  
  classes <- xl[1:l, n] # Классы объектов выборки
  weights <- table(classes) # Таблица весов классов
  weights[1:length(weights)] <- 0
  
  for (i in 1:l) { # Для всех объектов выборки
    class <- xl[i, n] # Берём его класс
    r <- distances[i] / h
    weights[class] <- weights[class] + type_core(r) # И прибавляем его вес к общему весу его класса
  }
  
  if (max(weights) != 0) # Если веса точки по классам не равны 0 (точка попала в окно)
    return (names(which.max(weights))) # Вернуть класс с максимальным весом
  else
    return (0) # Иначе - вернуть 0
}

LOO <- function(xl, type_core) {
  l <- nrow(xl)
  n <- ncol(xl)
  h_temp <- seq(0.1, 2, 0.1)
  sum <- rep(0, length(h_temp))
  
  for (i in 1:l) {
    cnt <- 1
    xi <- xl[i, 1:(n-1)]
    xt <- xl[-i, ]
    
    distances <- Distances(xt, xi)
    for (h in h_temp) {
      class <- parzen(xt, h, distances, type_core)
      if (class != xl[i, n] || class == 0) {
        sum[cnt] = sum[cnt] + 1/l
      }
      cnt <- cnt + 1
    }
  }
  return (sum)
}

Opt_H <- function(LOO_H) which.min(LOO_H) / 10

ClassMap <- function(xl, h, type_core) {
  l <- nrow(xl)
  n <- ncol(xl)
  ox <- seq(0, 7, 0.1)
  oy <- seq(0, 2.5, 0.1)
  classMatrix <- matrix(NA, length(ox)*length(oy), n)
  cnt <- 1
  for (i in ox)
    for (j in oy) {
      z <- c(i, j)
      distances <- Distances(xl, z)
      class <- parzen(xl, h, distances, type_core)
      if (class != 0) {
        classMatrix[cnt, ] <- c(i, j, class)
        cnt <- cnt + 1
      }
    }
  return (classMatrix)
}

DrowPlots <- function(xl, classMatrix, LOO_H, h) {
  l <- nrow(classMatrix)
  n <- ncol(classMatrix)
  colors <- c("setosa" = "red", "versicolor" = "green3", "virginica" = "blue")
  h10 <- h*10
  par(mfrow=c(1, 2))
  
  # Карта классификации
  plot(xl[,1:(n-1)], pch = 21, bg = colors[xl[, n]], col = colors[xl[, n]], main = "Метод парзеновского окна", xlab = "Длина лепестка", ylab = "Ширина лепестка", asp = 1)
  points(classMatrix[, 1:(n-1)], pch = 1, col = colors[classMatrix[, n]])
  
  # График lOO
  plot(seq(0.1, 2, 0.1), LOO_H[1:length(LOO_H)], type = "l", bg = "red", col = "red", main = "Оценка оптимальности h по LOO", xlab = "h", ylab = "LOO")
  points(h, LOO_H[h10], pch = 21, bg = "blue", col = "blue")
  label <- paste("h = ", h, "\n", "LOO = ", round(LOO_H[h10], 3))
  text(h, LOO_H[h10], labels = label, pos = 3)
}

main <- function(type_core) {
  xl <- iris[, 3:5]
  LOO_H <- LOO(xl, type_core)
  h <- Opt_H(LOO_H)
  classMatrix <- ClassMap(xl, h, type_core)
  DrowPlots(xl, classMatrix, LOO_H, h)
}

main(core.P)