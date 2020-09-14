nma.networkplot <- function(s.id, t.id, n, data, title = "", title.cex = 1, trtname,
  alphabetic = TRUE, multi.show = FALSE, multi.col, weight.edge = TRUE, adjust.thick = 5,
  weight.node = TRUE, weight.node.ss = FALSE, adjust.node.size = 10,
  node.col = "gray", edge.col = "black", text.cex = 1,
  adjust.figsizex = 1.1, adjust.figsizey = 1.1){
  if(missing(s.id) | missing(t.id)){
    stop("both study id and treatment id must be specified.")
  }
  if(!missing(data)){
    s.id <- eval(substitute(s.id), data, parent.frame())
    t.id <- eval(substitute(t.id), data, parent.frame())
    if(weight.node & weight.node.ss){
      if(missing(n)){
        stop("sample sizes are needed for weights of node sizes.")
      }
      n <- eval(substitute(n), data, parent.frame())
    }
  }

  unique.sid <- unique(s.id)
  nstudy <- length(unique.sid)
  sid.ori <- s.id
  for(s in 1:nstudy){
    s.id[sid.ori == unique.sid[s]] <- s
  }

  unique.tid <- sort(unique(t.id))
  ntrt <- length(unique.tid)
  if(ntrt <= 2) stop("there are less than 3 treatments, no need for network plot.")
  if(!missing(trtname) & alphabetic){
    trtname.order <- order(trtname)
    unique.tid <- unique.tid[trtname.order]
  }
  if(!missing(trtname) & !alphabetic){
    trtname.order <- 1:ntrt
  }
  if(missing(trtname)){
    trtname <- unique.tid
    trtname.order <- 1:ntrt
  }
  if(length(trtname) != ntrt){
    stop("the length of trtname do not match the data.")
  }
  ## make treatment id to be from 1 to ntrt
  tid.ori <- t.id
  for(t in 1:ntrt){
    t.id[tid.ori == unique.tid[t]] <- t
  }

  pi <- asin(1)*2
  polar <- pi/2 - 2*pi/ntrt*(0:(ntrt - 1))
  x <- cos(polar)
  y <- sin(polar)
  plot(0, 0, axes = FALSE, xlab="", ylab="", type = "n",
    xlim = c(-adjust.figsizex, adjust.figsizex),
    ylim = c(-adjust.figsizey, adjust.figsizey),
    main = title, cex.main = title.cex)

  if(multi.show){
    multi <- NULL
    multi.sid <- table(s.id)
    if(any(multi.sid > 2)){
      multi <- list(NULL)
      multi.sid <- names(multi.sid)[multi.sid > 2]
      for(j in 1:length(multi.sid)){
        multi.temp <- t.id[s.id == multi.sid[j]]
        if(j == 1){
          multi[[j]] <- multi.temp
        }else{
          check <- 1
          for(k in 1:length(multi)){
            if(length(multi.temp) == length(multi[[k]])){
              if(all(sort(multi.temp) - sort(multi[[k]]) == 0)) check <- 0
            }
          }
          if(check == 1){
            multi[[length(multi) + 1]] <- multi.temp
          }
        }
      }
    }
    if(!is.null(multi)){
      if(missing(multi.col)){
        if(length(multi) <= 13){
          multi.col <- c("red", "yellow", "blue", "green4", "orange",
            "purple", "palegreen4", "saddlebrown", "tan", "khaki",
            "bisque", "deeppink", "lightsalmon")
        }else{
          multi.col <- rainbow(length(multi))
        }
        multi.col <- adjustcolor(multi.col, alpha.f = 0.2)
      }else{
        if(length(multi.col) < length(multi)){
          stop("more colors are needed for visualizing multi-arm studies.")
        }
      }
      all.multi.cols <- multi.col
      multi.col <- all.multi.cols[1:length(multi)]
      for(i in 1:length(multi)){
        multi.temp <- multi[[i]]
        multi.temp <- sort(multi.temp)
        multi.id <- which(is.element(1:ntrt, multi.temp))
        polygon(x = x[multi.id], y = y[multi.id], border = NA, col = multi.col[i])
      }
    }
  }

  wt <- matrix(0, ntrt, ntrt)
  for(t1 in 2:ntrt){
    for(t2 in 1:(t1 - 1)){
      study.t1 <- s.id[t.id == t1]
      study.t2 <- s.id[t.id == t2]
      study.t1.t2 <- intersect(study.t1, study.t2)
      wt[t1, t2] <- length(study.t1.t2)
    }
  }
  wt <- c(wt)
  if(weight.edge == TRUE){
    wt.unique <- unique(wt[wt > 0])
    wtmin <- min(wt.unique)
    wtmax <- max(wt.unique)
    if(wtmin < wtmax){
      wt[wt > 0] <- round(1 + adjust.thick*(wt[wt > 0] - wtmin)/(wtmax - wtmin))
    }else{
      wt[wt > 0] <- 2
    }
  }else{
    wt[wt > 0] <- 2
  }
  wt <- matrix(wt, ntrt, ntrt)

  for(t1 in 2:ntrt){
    for(t2 in 1:(t1 - 1)){
      if(t1 != t2 & wt[t1, t2] > 0){
        lines(x = x[c(t1, t2)], y = y[c(t1, t2)], lwd = wt[t1, t2], col = edge.col)
      }
    }
  }
  if(weight.node){
    if(weight.node.ss){
      wt <- rep(NA, ntrt)
      for(i in 1:ntrt){
        wt[i] <- sum(n[t.id == i])
      }
    }else{
      wt <- matrix(0, ntrt, ntrt)
      for(t1 in 2:ntrt){
        for(t2 in 1:(t1 - 1)){
          study.t1 <- s.id[t.id == t1]
          study.t2 <- s.id[t.id == t2]
          study.t1.t2 <- intersect(study.t1, study.t2)
          wt[t1, t2] <- length(study.t1.t2)
        }
      }
      wt <- wt + t(wt)
      wt <- colSums(wt)
    }
    node.sizes <- 3 + (wt - min(wt))/(max(wt) - min(wt))*adjust.node.size
    points(x, y, pch = 20, cex = node.sizes, col = node.col)
  }else{
    points(x, y, pch = 20, cex = 3, col = node.col)
  }

  sides <- numeric(ntrt)
  eps <- 10^(-4)
  for(t in 1:ntrt){
    if((polar[t] <= pi/2 & polar[t] > pi/4) |
      (polar[t] < -5*pi/4 & polar[t] >= -3*pi/2)){
      sides[t] <- 3
    }
    if(polar[t] <= pi/4 & polar[t] >= -pi/4){
      sides[t] <- 4
    }
    if(polar[t] < -pi/4 & polar[t] > -3*pi/4){
      sides[t] <- 1
    }
    if(polar[t] <= -3*pi/4 & polar[t] >= -5*pi/4){
      sides[t] <- 2
    }
  }
  for(t in 1:ntrt){
    if(weight.node){
      text(x = x[t], y = y[t], labels = trtname[trtname.order[t]], cex = text.cex)
    }else{
      text(x = x[t], y = y[t], labels = trtname[trtname.order[t]], pos = sides[t], cex = text.cex)
    }
  }
}