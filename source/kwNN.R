euclide <- function(u, v) sqrt(sum((u - v)^2))

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

kwNN <- function(xl, z, k, q)
{
  #print("kwNN")
  orderedXl <- sortObjectByDist(xl, z)
  n <- dim(xl)[2]
  
  classes <- orderedXl[1:k, n]
  
  classes <- table(classes)        
  classes[1:length(classes)] <- 0
  for (i in names(classes)) {
    for (j in 1:k) {
      if (orderedXl[j, n] == i)
        classes[i] = classes[i] + q^j
    }
  }
  class <- names(which.max(classes))
  return (class)
}

kwNN_for_sorted_array <- function(xl, k, q)
{
  #print("kwNN_sort")
  n <- dim(xl)[2]
  
  classes <- xl[1:k, n]
  
  classes <- table(classes)        
  classes[1:length(classes)] <- 0
  for (i in names(classes)) {
    for (j in 1:k) {
      if (xl[j, n] == i)
        classes[i] = classes[i] + q^j
    }
  }
  class <- names(which.max(classes))
  return (class)
}

LOO <- function(xl)                             
{
  print("LOO")
  l <- nrow(xl)
  n <- ncol(xl)
  qRange <- seq(0.1, 1, 0.1)                    
  LOO_K <- matrix(0, l-1, length(qRange))
  for (i in 1:l) {
    xt <- xl[i, 1:2]                            
    orderedXl <- sortObjectByDist(xl[-i, ], xt) 
    for (k in 1:(l-1)) {
      q_cnt <- 1
      for (q in qRange) {
        class <- kwNN_for_sorted_array(orderedXl, k, q)
        if (class != xl[i, n])
          LOO_K[k, q_cnt] <- LOO_K[k, q_cnt] + 1 / l
        q_cnt <- q_cnt + 1
      }
    }
  }
  return (LOO_K) 
}

Opt_K <- function(LOO_K) 
{ 
  print("Opt_K")
  optIndex <- 1
  optVal <- LOO_K[1, 1]
  for (i in 1:ncol(LOO_K)) {
    minIndex <- which.min(LOO_K[, i])
    minVal <- LOO_K[minIndex, i]
    
    if (optVal > minVal) {
      optIndex <- minIndex
      optVal <- minVal
    }
  }
  return (optIndex)
}

Opt_Q <- function(k, LOO_K) 
{
  print("Opt_Q")
  optVal <- LOO_K[k, 1]
  optIndex <- 1
  for (i in 2:ncol(LOO_K)) {
    minValue <- LOO_K[k, i]
    if (optVal > minValue) {
      optVal <- minValue
      optIndex <- i
    }
  }
  return (optIndex / 10) 
}

drowPoints <- function(xl, LOO_K, k, q)
{
  print("Field")
  l <- nrow(xl)
  n <- ncol(xl)
  cnt <- 1
  
  classField <- matrix(NA, length(seq(0, 7, 0.1)) * length(seq(0, 2.5, 0.1)), n)
  
  for(i in seq(0, 7, 0.1)){
    for(j in seq(0, 2.5, 0.1))
    {
      z <- c(i, j)
      class <- kwNN(xl, z, k, q)
      classField[cnt, ] <- c(i, j, class)
      cnt <- cnt + 1
    }
  }
  print("Drow")
  
  colors <- c("setosa" = "red", "versicolor" = "green3", "virginica" = "blue")
  q10 = q * 10
  par(mfrow=c(1, 2))
  
  plot(iris[, 3:4], pch = 21, bg = colors[iris$Species], col = colors[iris$Species], main="Êëàññèôèêàöèÿ èðèñîâ Ôèøåðà ìåòîäîì kwNN", xlab = "Äëèíà ëåïåñòêà", ylab = "Øèðèíà ëåïåñòêà", asp = 1)
  points(classField[, 1:(n-1)], pch = 1, bg = colors[classField[, n]], col = colors[classField[, n]])
  
  plot(LOO_K[1:nrow(LOO_K), q10], type = "l", bg = "red", col = "red", main = "Îöåíêà îïòèìàëüíîñòè  k ïî LOO", xlab = "Çíà÷åíèÿ k", ylab = "Çíà÷åíèÿ LOO")
  points(k, LOO_K[k, q10], pch = 21, bg = "blue", col = "blue")
  label = paste("k = ", k, "\n", "LOO = ", round(LOO_K[k, q10], 3))
  text(k, LOO_K[k, q10], labels = label, pos = 3)
  lines(LOO_K, col = "red")  
}

main <- function()
{
  print("main")
  xl <- iris[, 3:5]
  LOO_K <- LOO(xl)
  k <- Opt_K(LOO_K)
  q <- Opt_Q(k, LOO_K)
  drowPoints(xl, LOO_K, k, q)
}

main()
