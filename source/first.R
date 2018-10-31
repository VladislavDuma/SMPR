

euclide <- function(u, v)
{
  sqrt(sum((u - v)^2))
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


oneNN <- function(xl, z)
{
  orderedXl <- sortObjectByDist(xl, z)
  n <- dim(orderedXl)[2]
  
  classes <- orderedXl[1, n]
  
  return (classes)
}

xl_arr <- function(q)
{
  if (q > 50) 
    q = 50
    
  res1 <- sample(c(1:50), q, replace = FALSE)
  res2 <- sample(c(51:100), q, replace = FALSE)
  res3 <- sample(c(101:150), q, replace = FALSE)
  x1 <- iris[res1, 3:5]
  x2 <- iris[res2, 3:5]
  x3 <- iris[res3, 3:5]
  arr <- rbind(x1,x2,x3)
  return(arr)
}

drowPoints <- function(xl)
{
  colors <- c("setosa" = "red", "versicolor" = "green3", "virginica" = "blue")
  plot(xl[1:2], pch = 21, bg = colors[xl$Species], col = colors[xl$Species], asp = 1)
  
  for(i in seq(0, 7, 0.1)){
    for(j in seq(0, 2.5, 0.1))
    {
      z <- c(i, j)
      class <- oneNN(xl, z)
      points(z[1], z[2], pch = 1, bg = colors[class], col = colors[class])
    }
  }
}

main <- function()
{
  xl <- xl_arr(50)
  drowPoints(xl)
} # end main function

main()


