
euclide <- function(u, v)
{
  #return (1)
  sqrt(sum(u - v)^2)
}


sortObjectByDist <- function(xl, z, metricFunction = euclide)
{
  l <- dim(xl)[1]
  n <- dim(xl)[2] - 1
  
  distances <- rep(0, l) #вектор на L элементов, заполненный 0
  for (i in 1:l)
  {
    distances[i] <- c(metricFunction(xl[i, 1:n], z)) # вектор с расстояниями
  }
  orderedXl <- xl[order(distances), ] # сортирует вектор дистанции по индексам
  return (orderedXl)
}


kNN <- function(xl, z, k)
{
  orderedXl <- sortObjectByDist(xl, z)
  n <- dim(orderedXl)[2]
  
  classes <- orderedXl[1:k, n] # убрал -1 & +1, так как даёт доступ к тем же данным
  
  counts <- table(classes)
  class <- names(which.max(counts))
  return (class)
}
q <- 15 # количество элементов выборки из каждого класса

res1 <- sample(c(1:50), q, replace = FALSE)
res2 <- sample(c(51:100), q, replace = FALSE)
res3 <- sample(c(101:150), q, replace = FALSE)
x1 <- iris[res1, 3:5]
x2 <- iris[res2, 3:5]
x3 <- iris[res3, 3:5]
xl <- rbind(x1,x2,x3)
colors <- c("setosa" = "red", "versicolor" = "green3", "virginica" = "blue")
#plot(iris[, 3:4], pch = 21, bg = colors[iris$Species], col = colors[iris$Species], asp = 1)
plot(xl[1:2], pch = 20, bg = colors[xl$Species], col = colors[xl$Species], asp = 1)


#for(i in seq(0, 7, 0.1)) {
#  for(j in seq(0, 2.5, 0.1)){
#    z <- c(i, j)
#    class <- kNN(xl, z, k = 5)
#    points(z[1], z[2], pch = 1, bg = colors[class], col = colors[class])
#  }
#}

#for(i in 1:30){  points(xl[1:2], pch = 22, bg = colors[xl$Species], col = colors[xl$Species])}
len <- dim(xl)[1]
res <- rep(0, 150)

LOO <- function(xl, k){
  
  for(i in 1:len){
    
    xt <- xl[-i, ]
    z <- xl[i, 1:2]
    q <- xl[i, 3]
    #print(xt)
    #print(z)
    class <- kNN(xt, z, k)
    if (class != q)
      res[i] <- res[i]+1
    #points(z[1], z[2], pch = 1, bg = colors[class], col = colors[class])
    
  }# завершение цикла в LOO
  
  #print(len)
  #print(xl)
  #print(xt)
} #конец функции



for(k in 1:len)
{
  LOO(xl, k)
}
result <- res/150

