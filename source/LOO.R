
euclide <- function(u, v)
{
  sqrt(sum((u - v)^2))
}

sortObjectByDist <- function(xl, z, metricFunction = euclide)
{
  l <- dim(xl)[1]
  n <- dim(xl)[2] - 1
  
  distances <- rep(0, l)
  for (i in 1:l)
  {
    distances[i] <- c(metricFunction(xl[i, 1:n], z)) # array distances
  }
  orderedXl <- xl[order(distances), ] # sorted by index
  return (orderedXl)
}

kNN <- function(xl, z, k)
{
  orderedXl <- sortObjectByDist(xl, z)
  n <- dim(xl)[2]
  
  classes <- orderedXl[1:k, n]
  
  counts <- table(classes)
  class <- names(which.max(counts))
  return (class)
}

kNN_for_sorted_array <- function(xl, z, k)
{
  n <- dim(xl)[2]
  
  classes <- xl[1:k, n]
  
  counts <- table(classes)
  class <- names(which.max(counts))
  return (class)
}

LOO <- function(xl, len){ # LOO function
  
  res <- rep(0, len)
  
  for(i in 1:len){
    
    xt <- xl[-i, ]
    z <- xl[i, 1:2]
    q <- xl[i, 3]
    
    x_sort <- sortObjectByDist(xt, z)
    for(k in 1:len)
    {
      class <- kNN_for_sorted_array(x_sort, z, k)
      if (class != q)
        res[k] <- res[k]+1
    }
  }
  
  result <- res/len
  
  plot(1:len, result, type = "l", col = "red", main = "LOO äëÿ kNN",xlab = "Çíà÷åíèå k",ylab = "Çíà÷åíèå LOO")
  
  min = which.min(result)
  label = paste("k = ", min, "\n", "LOO = ", result[min])
  text(min, result[min], labels = label, pos = 4)
  
} # end LOO

main <- function()
{
  xl <- iris[, 3:5]
  k <- dim(xl)[1]
  LOO(xl, k)
  
} # end main function

main()






