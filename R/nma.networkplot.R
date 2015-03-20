nma.networkplot <- function(s.id, t.id, title = "", trtname,
	alphabetic = TRUE, weight = TRUE, adjust.thick = 5,
	node.size = 3, node.col = "blue", edge.col = "black", text.cex = 1,
	adjust.figsizex = 1.1, adjust.figsizey = 1.1){
	if(missing(s.id) | missing(t.id)){
		stop("both study id and treatment id must be specified.")
	}

	unique.sid <- unique(s.id)
	nstudy <- length(unique.sid)
	sid.ori <- s.id
	for(s in 1:nstudy){
		s.id[sid.ori == unique.sid[s]] <- s
	}

	unique.tid <- sort(unique(t.id))
	ntrt <- length(unique.tid)
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
	## make treatment id to be 1 to ntrt
	tid.ori <- t.id
	for(t in 1:ntrt){
		t.id[tid.ori == unique.tid[t]] <- t
	}

	polar <- pi/2 - 2*pi/ntrt*(0:(ntrt - 1))
	x <- cos(polar)
	y <- sin(polar)
	plot(x, y, axes = FALSE, xlab="", ylab="", cex = 0.1,
		xlim = c(-adjust.figsizex, adjust.figsizex),
		ylim = c(-adjust.figsizey, adjust.figsizey),
		main = title)

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
	if(weight == TRUE){
		wt.unique <- unique(wt[wt > 0])
		wtmin <- min(wt.unique)
		wtmax <- max(wt.unique)
		wt[wt > 0] <- round(1 + adjust.thick*(wt[wt > 0] - wtmin)/
			(wtmax - wtmin))
	}else{
		wt[wt > 0] <- 2
	}
	wt <- matrix(wt, ntrt, ntrt)
	
	if(length(unique(c(wt)[c(wt) > 0])) <= 10) 

	for(t1 in 2:ntrt){
		for(t2 in 1:(t1 - 1)){
			if(t1 != t2 & wt[t1, t2] > 0){
				lines(x = x[c(t1, t2)], y = y[c(t1, t2)],
					lwd = wt[t1, t2], col = edge.col)
			}
		}
	}

	points(x, y, pch = 20, cex = node.size, col = node.col)

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
		text(x = x[t], y = y[t], labels = trtname[trtname.order[t]],
			pos = sides[t], cex = text.cex)
	}
}