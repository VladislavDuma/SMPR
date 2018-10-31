
euclide <- function(u, v)
{
  sqrt(sum((u - v)^2))
}

sortObjectByDist <- function(xl, z, metricFunction = euclide) # Sort
{
  l <- dim(xl)[1]
  n <- dim(xl)[2] - 1
  
  distances <- rep(0, l)
  for (i in 1:l)
  {
    distances[i] <- c(metricFunction(xl[i, 1:n], z))
  }
  orderedXl <- xl[order(distances), ] 
  return (orderedXl)
}

kNN <- function(xl, z, k) #kNN method
{
  orderedXl <- sortObjectByDist(xl, z)
  n <- dim(xl)[2]
  
  classes <- orderedXl[1:k, n]
  
  counts <- table(classes)
  class <- names(which.max(counts))
  return (class)
}

drowPoints <- function(xl)
{
  colors <- c("setosa" = "red", "versicolor" = "green3", "virginica" = "blue")
  plot(xl[1:2], pch = 21, bg = colors[xl$Species], col = colors[xl$Species],main = "Ìåòîä êëàññèôèêàöèè kNN ïðè k = 10", asp = 1)
  
  for(i in seq(0, 7, 0.1)){
    for(j in seq(0, 2.5, 0.1))
    {
      z <- c(i, j)
      class <- kNN(xl, z, 10)
      points(z[1], z[2], pch = 1, bg = colors[class], col = colors[class])
    }
  }
}

main <- function()
{
  xl <- iris[, 3:5]
  drowPoints(xl)
} # end main function

main()






