norm <- function(x){
  r <- apply(x^2, 1, sum)
  sqrt(r)
}

mc.knn <- function(xl, u, k, sorted=FALSE, metric=norm){
  cols <- ncol(xl)
  rows <- nrow(xl)
  if(sorted != TRUE){
     umat <- matrix(rep(u, rows), rows, cols-1, byrow=TRUE)
     xl <- xl[order(metric(umat - xl[,1:(cols-1)] )),]
  }
  xl[which.max(table(xl[1:k,cols])), cols]
}

mc.kwnn <- function(xl, u, k, wf, sorted=FALSE, metric=norm){
  cols <- ncol(xl)
  rows <- nrow(xl)
  
  if(sorted != TRUE){
    umat <- matrix(rep(u, rows), rows, cols-1, byrow=TRUE)
    w <- apply(data.frame(k, 1:k), 1, wf)
    w <- matrix(w, length(w), 1)
    xl <- xl[order(metric(umat - xl[,1:(cols-1)] )),]
  }
  xl <- xl[1:k,cols]
  
  classes <- names(table(xl))
  score <- rep(0, length(classes))
  i <- 1
  
  for(el in xl){
    class <- xl[i]
    score[class] <- score[class]+w[i]
    i <- i+1
  }

  classes[which.max(score)]
}

mc.window <- function(xl, u, h, K, sorted=FALSE, metric=norm){
  cols <- ncol(xl)
  rows <- nrow(xl)
  
  umat <- matrix(rep(u, rows), rows, cols-1, byrow=TRUE)
  distances <- metric(umat - xl[,1:(cols-1)])
  if(sorted != TRUE){
    orderedIndexes <- order(distances)
    xl <- xl[orderedIndexes,]
    distances <- distances[orderedIndexes]
  }
  xl <- xl[,cols]
  wp <- distances/h
  
  classes <- c(names(table(xl)), "none")
  score <- rep(0, length(classes))
  score[length(classes)]<-0.00001
  i <- 1
  
  for(el in xl){
    class <- xl[i]
    score[class] <- score[class]+K(wp[i])
    i <- i+1
  }
  
  classes[which.max(score)]
}

mc.window.auto <- function(xl, u, k, K, sorted=FALSE, metric=norm){
  cols <- ncol(xl)
  rows <- nrow(xl)
  
  umat <- matrix(rep(u, rows), rows, cols-1, byrow=TRUE)
  distances <- metric(umat - xl[,1:(cols-1)])
  if(sorted != TRUE){
    orderedIndexes <- order(distances)
    xl <- xl[orderedIndexes,]
    distances <- distances[orderedIndexes]
  }
  xl <- xl[,cols]
  wp <- distances[1:k]/distances[k+1]
  
  classes <- c(names(table(xl)), "none")
  score <- rep(0, length(classes))
  score[length(classes)]<-0.00001
  i <- 1
  
  for(el in xl){
    class <- xl[i]
    score[class] <- score[class]+K(wp[i])
    i <- i+1
  }
  
  classes[which.max(score)]
}

setka <- function(f, xl, color, mi, ma, acc, ...){
  for(x in seq(mi[1], ma[1], acc)){
    for(y in seq(mi[2], ma[2], acc)){
      u <- c(x, y)
      class <- f(xl, u, ...)
      points(x, y, col=color[class], pch="+")
    }
  }
}

loo <-function(f, xl, ...){
  cols <- ncol(xl)
  rows <- nrow(xl)
  acc <- 0
  
  for(i in 1:rows){
    xli <- xl[-i,]
    u <- unname(unlist(xl[i, 1:(cols-1)]))
    class <- f(xli, u, ...)
    if(xl[i,cols] == class){
      acc <- acc+1
    }
  }
  
  1-acc/rows
}

loo.list <- function(f, xl, lstOfX,...){
  p <- c()
  for(x in lstOfX){
      y <- loo(f, xl, x, ...)
      p <- rbind(p, c(x,y))
  }
  p
}

p <- loo.list(mc.window, sel, seq(from=0.1, to=5, by=0.1), mc.K.T)
plot(p, type="l", xlab="k in knn", ylab="errors using loo")


mc.wlin <- function(el){
  k <- el[1]
  i <- el[2]
  (k+1-i)/k
}

mc.K.T <- function(r){
  if(abs(r)<=1)
    1-abs(r)
  else
    0
}

mc.K.P <- function(r){
  if(abs(r)<=1)
    1/2
  else
    0
}

mc.K.G <- function(r){
    (2*pi)^(1/2 * exp(-1/2 * r^2))
}

mc.K.E <- function(r){
  if(abs(r)<=1)
    3/4 * (1-r^2)
  else
    0
}

mc.K.Q <- function(r){
  if(abs(r)<=1)
    15/16*(1-r^2)^2
  else
    0
}

#init
sel<-iris[,3:5]
pl<-sel[,1:2]
mi <- c(1, 0.1)
ma <- c(7, 2.5)
acc <- 3
colors <- c(setosa="red", versicolor="green", virginica="blue")
par("lwd")

#1NN setka with LOO
plot(pl, col=colors[sel$Species], main="1NN", pch=19)
k <- 1
setka(mc.knn, sel, colors, mi, ma, 0.1, k)
legend(mi[1], ma[2], legend=paste("LOO =", round(loo(mc.knn, sel, 1), digits=3)))

#LOO(k) in KNN
p<-loo.list(mc.knn, sel, 1:50)
plot(p, type="l", xlab="k", ylab="error", main="LOO KNN")
legend(0, 0.066, legend=paste("LOO(6) =", round(loo(mc.knn, sel, 6), digits=3)))

#6NN setka with LOO
plot(pl, col=colors[sel$Species], main="6NN", pch=19)
k <- 6
setka(mc.knn, sel, colors, mi, ma, 0.1, k)
legend(mi[1], ma[2], legend=paste("LOO =", round(loo(mc.knn, sel, k), digits=3)))

#LOO in kwnn
p<-loo.list(mc.kwnn, sel, 1:50, mc.wlin)
plot(p, type="l", xlab="k", ylab="error", main="LOO KWNN")
legend(0, 0.066, legend=paste("LOO(6) =", round(loo(mc.knn, sel, 6), digits=3)))

#6WNN setka with LOO
plot(pl, col=colors[sel$Species], main="6WNN", pch=19)
k <- 6
setka(mc.kwnn, sel, colors, mi, ma, 0.1, k, mc.wlin)
legend(mi[1], ma[2], legend=paste("LOO =", round(loo(mc.knn, sel, k), digits=3)))

#LOO in kwnn
hs<-seq(from=0.1, to=2, by=0.1)
p<-loo.list(mc.window, sel, hs, mc.K.T)
plot(p, type="l", xlab="k", ylab="error", main="LOO PARZEN WINDOW", ylim=c(0, 0.2))
lines(loo.list(mc.window, sel, hs, mc.K.P),  col = "red")
lines(loo.list(mc.window, sel, hs, mc.K.E),  col = "green")
lines(loo.list(mc.window, sel, hs, mc.K.G),  col = "blue")
lines(loo.list(mc.window, sel, hs, mc.K.Q),  col = "grey")
legend(1.8, 0.2, legend=c("K=T", "K=P", "K=E", "K=G", "K=Q"), lty=5, col=c("black", "red", "green", "blue", "grey") )

hs<-1:20
p<-loo.list(mc.window.auto, sel, hs, mc.K.T)
plot(p, type="l", xlab="k", ylab="error", main="LOO PARZEN WINDOW", ylim=c(0, 0.2))
lines(loo.list(mc.window.auto, sel, hs, mc.K.P),  col = "red")
lines(loo.list(mc.window.auto, sel, hs, mc.K.E),  col = "green")
lines(loo.list(mc.window.auto, sel, hs, mc.K.G),  col = "blue")
lines(loo.list(mc.window.auto, sel, hs, mc.K.Q),  col = "grey")
legend(1.8, 0.2, legend=c("K=T", "K=P", "K=E", "K=G", "K=Q"), lty=5, col=c("black", "red", "green", "blue", "grey") )

#6WNN setka with LOO
plot(pl, col=colors[sel$Species], main="6WNN", pch=19)
k <- 6
setka(mc.kwnn, sel, colors, mi, ma, 0.1, k, mc.wlin)
legend(mi[1], ma[2], legend=paste("LOO =", round(loo(mc.knn, sel, k), digits=3)))


plot(pl, col=colors[sel$Species], pch=19)
k <- 6
w <- apply(data.frame(k, 1:k), 1, mc.w )
loo(mc.knn, sel, k)
setka(mc.window.auto, sel, colors, c(1,0.1), c(7,2.5), 0.1, k, mc.K.T)
setka(mc.kwnn, sel, colors, c(1,0.1), c(7,2.5), 0.1, k, mc.w)

