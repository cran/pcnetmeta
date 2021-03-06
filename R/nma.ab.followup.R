nma.ab.followup <-
function(s.id, t.id, event.n, total.n, followup, data, trtname, param = c("lograte", "logratio", "rank.prob"),
  model = "het_cor", prior.type, a = 0.001, b = 0.001, c = 10, higher.better = FALSE, digits = 4,
  n.adapt = 5000, n.iter = 100000, n.burnin = floor(n.iter/2), n.chains = 3,
  n.thin = max(1,floor((n.iter - n.burnin)/100000)), conv.diag = FALSE, trace = NULL,
  dic = FALSE, postdens = FALSE, mcmc.samples = FALSE){
  ## check the input parameters
  options(warn = 1)
  if(missing(s.id)) stop("please specify study id.")
  if(missing(t.id)) stop("please specify treatment.")
  if(missing(event.n)) stop("please specify event number.")
  if(missing(total.n)) stop("please specify total number of participants.")
  if(missing(followup)) stop("please specify the follow-up times for different studies.")
  if(!missing(data)){
    s.id <- eval(substitute(s.id), data, parent.frame())
    t.id <- eval(substitute(t.id), data, parent.frame())
    event.n <- eval(substitute(event.n), data, parent.frame())
    total.n <- eval(substitute(total.n), data, parent.frame())
    followup <- eval(substitute(followup), data, parent.frame())
  }
  if(length(s.id) != length(t.id) | length(t.id) != length(event.n) | length(event.n) != length(total.n) | length(total.n) != length(followup)){
    stop("s.id, t.id, event.n, total.n, and followup have different lengths.")
  }
  if(!all(followup > 0)) stop("follow-up times must be positive.")
  if(!is.element(model, c("hom_eqcor", "het_eqcor", "het_cor"))) stop("model should be specified as \"hom_eqcor\", \"het_eqcor\", or \"het_cor\".")

  if(any(is.na(s.id)) | any(is.na(t.id)) | any(is.na(event.n)) | any(is.na(total.n)) | any(is.na(followup))){
    dat <- cbind(s.id, t.id, event.n, total.n, followup)
    s.id <- s.id[complete.cases(dat)]
    t.id <- t.id[complete.cases(dat)]
    event.n <- event.n[complete.cases(dat)]
    total.n <- total.n[complete.cases(dat)]
    followup <- followup[complete.cases(dat)]
    cat("NA is not allowed in the input data set;\n")
    cat("the rows containing NA are removed.\n")
  }

  ## make ids consecutive
  s.id.o <- s.id
  t.id.o <- t.id
  s.label <- sort(unique(s.id.o))
  t.label <- sort(unique(t.id.o))
  nstudy <- length(s.label) # total number of studies
  ntrt <- length(t.label) # total number of treatments
  len <- length(s.id)
  s.id <- numeric(nstudy)
  for(i in 1:nstudy){
    s.id[which(s.id.o == s.label[i])] <- i
  }
  t.id <- numeric(ntrt)
  for(i in 1:ntrt){
    t.id[which(t.id.o == t.label[i])] <- i
  }

  if(missing(trtname)){
    if(is.numeric(t.id.o)){
      trtname <- paste("Trt", t.label, sep = "")
    }
    if(is.character(t.id.o)){
      trtname <- t.label
    }
  }
  if(length(trtname) != length(unique(t.id))) stop("the number of treatment names does not match for specified treatment id.")
  if(missing(prior.type)) prior.type <- ifelse(model == "het_cor", "invwishart", "unif")

  ## JAGS model
  if(model == "hom_eqcor"){
    modelstring <- model.followup.hom.eqcor(prior.type, is.element("rank.prob", param))
  }
  if(model == "het_eqcor"){
    modelstring <- model.followup.het.eqcor(prior.type, is.element("rank.prob", param))
  }
  if(model == "het_cor"){
    modelstring <- model.followup.het.cor(prior.type, is.element("rank.prob", param))
  }

  ## JAGS data
  if(model == "hom_eqcor" | model == "het_eqcor"){
    if(prior.type == "unif"){
      if(is.element("rank.prob", param)) data.jags <- list(s = s.id, t = t.id, y = event.n, totaln = total.n, f = followup, len = len, nstudy = nstudy, ntrt = ntrt, zeros = rep(0, ntrt), c = c, higher.better = higher.better)
      if(!is.element("rank.prob", param)) data.jags <- list(s = s.id, t = t.id, y = event.n, totaln = total.n, f = followup, len = len, nstudy = nstudy, ntrt = ntrt, zeros = rep(0, ntrt), c = c)
    }
    if(prior.type == "invgamma"){
      if(is.element("rank.prob", param)) data.jags <- list(s = s.id, t = t.id, y = event.n, totaln = total.n, f = followup, len = len, nstudy = nstudy, ntrt = ntrt, zeros = rep(0, ntrt), a = a, b = b, higher.better = higher.better)
      if(!is.element("rank.prob", param)) data.jags <- list(s = s.id, t = t.id, y = event.n, totaln = total.n, f = followup, len = len, nstudy = nstudy, ntrt = ntrt, zeros = rep(0, ntrt), a = a, b = b)
    }
  }
  if(model == "het_cor"){
    if(prior.type == "invwishart"){
      I <- diag(ntrt)
      if(is.element("rank.prob", param)) data.jags <- list(s = s.id, t = t.id, y = event.n, totaln = total.n, f = followup, len = len, nstudy = nstudy, ntrt = ntrt, zeros = rep(0, ntrt), I = I, higher.better = higher.better)
      if(!is.element("rank.prob", param)) data.jags <- list(s = s.id, t = t.id, y = event.n, totaln = total.n, f = followup, len = len, nstudy = nstudy, ntrt = ntrt, zeros = rep(0, ntrt), I = I)
    }
    if(prior.type == "chol"){
      if(is.element("rank.prob", param)) data.jags <- list(s = s.id, t = t.id, y = event.n, totaln = total.n, f = followup, len = len, nstudy = nstudy, ntrt = ntrt, zeros = rep(0, ntrt), c = c, higher.better = higher.better)
      if(!is.element("rank.prob", param)) data.jags <- list(s = s.id, t = t.id, y = event.n, totaln = total.n, f = followup, len = len, nstudy = nstudy, ntrt = ntrt, zeros = rep(0,ntrt), c = c)
    }
  }

  ## JAGS initial value
  rng.seeds <- sample(1000000, n.chains)
  mu.init <- numeric(ntrt)
  for(i in 1:ntrt){
    mu.init[i] <- 0
  }
  init.jags <- list(NULL)
  if(model == "hom_eqcor"){
    if(prior.type == "unif"){
      for(ii in 1:n.chains){
        init.jags[[ii]] <- list(mu = mu.init, vi = matrix(0, nstudy, ntrt), sigma = c/2, rho = 0.5, .RNG.name = "base::Wichmann-Hill", .RNG.seed = rng.seeds[ii])
      }
    }
    if(prior.type == "invgamma"){
      for(ii in 1:n.chains){
        init.jags[[ii]] <- list(mu = mu.init, vi = matrix(0, nstudy, ntrt), inv.sig.sq = a/b, rho = 0.5, .RNG.name = "base::Wichmann-Hill", .RNG.seed = rng.seeds[ii])
      }
    }
  }
  if(model == "het_eqcor"){
    if(prior.type == "unif"){
      for(ii in 1:n.chains){
        init.jags[[ii]] <- list(mu = mu.init, vi = matrix(0, nstudy, ntrt), sigma = rep(c/2, ntrt), rho = 0.5, .RNG.name = "base::Wichmann-Hill", .RNG.seed = rng.seeds[ii])
      }
    }
    if(prior.type == "invgamma"){
      for(ii in 1:n.chains){
        init.jags[[ii]] <- list(mu = mu.init, vi = matrix(0, nstudy, ntrt), inv.sig.sq = rep(a/b, ntrt), rho = 0.5, .RNG.name = "base::Wichmann-Hill", .RNG.seed = rng.seeds[ii])
      }
    }
  }
  if(model == "het_cor"){
    if(prior.type == "invwishart"){
      for(ii in 1:n.chains){
        init.jags[[ii]] <- list(mu = mu.init, vi = matrix(0, nstudy, ntrt), T = (ntrt + 1)*I, .RNG.name = "base::Wichmann-Hill", .RNG.seed = rng.seeds[ii])
      }
    }
    if(prior.type == "chol"){
      for(ii in 1:n.chains){
        init.jags[[ii]] <- list(mu = mu.init, vi = matrix(0, nstudy, ntrt), sigma = rep(c/2, ntrt), psi = matrix(3.1415926/2, ntrt - 1, ntrt - 1), .RNG.name = "base::Wichmann-Hill", .RNG.seed = rng.seeds[ii])
      }
    }
  }

  ## parameters to be monitored in JAGS
  if(!is.element("lograte", param)) param <- c("lograte", param)
  if(!is.null(trace)){
    if(!any(is.element(trace, param))) stop("at least one effect measure in argument trace is not specified in argument param.")
  }
  monitor <- param[!is.element(param, c("ratio", "logratio"))]
  if(is.element("ratio", param)){
    for(ii in 1:ntrt){
      for(jj in 1:ntrt){
        if(ii != jj) monitor <- c(monitor, paste("ratio[", ii, ",", jj, "]", sep = ""))
      }
    }
  }
  if(is.element("logratio", param)){
    for(ii in 1:ntrt){
      for(jj in 1:ntrt){
        if(ii < jj) monitor <- c(monitor, paste("logratio[", ii, ",", jj, "]", sep = ""))
      }
    }
  }

  ## run JAGS
  cat("Start running MCMC...\n")
  jags.m <- jags.model(file = textConnection(modelstring), data = data.jags, inits = init.jags, n.chains = n.chains, n.adapt = n.adapt)
  update(jags.m, n.iter = n.burnin)
  jags.out <- coda.samples(model = jags.m, variable.names = monitor, n.iter = n.iter, thin = n.thin)
  smry <- summary(jags.out)
  smry <- cbind(smry$statistics[,c("Mean", "SD")], smry$quantiles[,c("2.5%", "50%", "97.5%")])
  smry <- signif(smry, digits = digits)

  out <- NULL
  out$model <- "Binomial likelihood with cloglog link."
  lograte.id <- which(is.element(rownames(smry), paste("lograte[", 1:ntrt, "]", sep = "")))
  lograte.stat <- array(paste(format(round(smry[lograte.id, "Mean"], digits = digits), nsmall = digits), " (", format(round(smry[lograte.id, "SD"], digits = digits), nsmall = digits), ")", sep = ""), dim = c(ntrt, 1))
  colnames(lograte.stat) <- "Mean (SD)"
  rownames(lograte.stat) <- trtname
  lograte.quan <- array(paste(format(round(smry[lograte.id, "50%"], digits = digits), nsmall = digits), " (", format(round(smry[lograte.id, "2.5%"], digits = digits), nsmall = digits),
    ", ", format(round(smry[lograte.id, "97.5%"], digits = digits), nsmall = digits), ")", sep = ""), dim = c(ntrt, 1))
  colnames(lograte.quan) <- "Median (95% CI)"
  rownames(lograte.quan) <- trtname
  out$LogRate <- list(Mean_SD = noquote(lograte.stat), Median_CI = noquote(lograte.quan))

  if(is.element("rate", param)){
    rate.id <- which(is.element(rownames(smry), paste("rate[", 1:ntrt, "]", sep = "")))
    rate.stat <- array(paste(format(round(smry[rate.id, "Mean"], digits = digits), nsmall = digits), " (", format(round(smry[rate.id, "SD"], digits = digits), nsmall = digits), ")", sep = ""), dim = c(ntrt, 1))
    colnames(rate.stat) <- "Mean (SD)"
    rownames(rate.stat) <- trtname
    rate.quan <- array(paste(format(round(smry[rate.id, "50%"], digits = digits), nsmall = digits), " (", format(round(smry[rate.id, "2.5%"], digits = digits), nsmall = digits),
      ", ", format(round(smry[rate.id, "97.5%"], digits = digits), nsmall = digits), ")", sep = ""), dim = c(ntrt, 1))
    colnames(rate.quan) <- "Median (95% CI)"
    rownames(rate.quan) <- trtname
    out$Rate <- list(Mean_SD = noquote(rate.stat), Median_CI = noquote(rate.quan))
  }

  if(is.element("ratio", param)){
    ratio.stat <- ratio.quan <- array(NA, dim = c(ntrt, ntrt))
    colnames(ratio.stat) <- colnames(ratio.quan) <- rownames(ratio.stat) <- rownames(ratio.quan) <- trtname
    for(i in 1:ntrt){
      ratio.stat[i,i] <- ratio.quan[i,i] <- "--"
      for(j in 1:ntrt){
        if(i != j){
          ratio.ij <- paste("ratio[", i, ",", j, "]", sep = "")
          ratio.stat[i,j] <- paste(format(round(smry[ratio.ij, "Mean"], digits = digits), nsmall = digits), " (", format(round(smry[ratio.ij, "SD"], digits = digits), nsmall = digits), ")", sep = "")
          ratio.quan[i,j] <- paste(format(round(smry[ratio.ij, "50%"], digits = digits), nsmall = digits), " (", format(round(smry[ratio.ij, "2.5%"], digits = digits), nsmall = digits),
            ", ", format(round(smry[ratio.ij, "97.5%"], digits = digits), nsmall = digits), ")", sep = "")
        }
      }
    }
    out$RateRatio <- list(Mean_SD = noquote(ratio.stat), Median_CI = noquote(ratio.quan))
  }

  if(is.element("logratio", param)){
    logratio.stat <- logratio.quan <- array(NA, dim = c(ntrt, ntrt))
    colnames(logratio.stat) <- colnames(logratio.quan) <- rownames(logratio.stat) <- rownames(logratio.quan) <- trtname
    for(i in 1:ntrt){
      logratio.stat[i,i] <- logratio.quan[i,i] <- "--"
      for(j in 1:ntrt){
        if(i < j){
          logratio.ij <- paste("logratio[", i, ",", j, "]", sep = "")
          logratio.stat[i,j] <- paste(format(round(smry[logratio.ij, "Mean"], digits = digits), nsmall = digits), " (", format(round(smry[logratio.ij, "SD"], digits = digits), nsmall = digits), ")", sep = "")
          logratio.stat[j,i] <- paste(format(round(-smry[logratio.ij, "Mean"], digits = digits), nsmall = digits), " (", format(round(smry[logratio.ij, "SD"], digits = digits), nsmall = digits), ")", sep = "")
          logratio.quan[i,j] <- paste(format(round(smry[logratio.ij, "50%"], digits = digits), nsmall = digits), " (", format(round(smry[logratio.ij, "2.5%"], digits = digits), nsmall = digits),
            ", ", format(round(smry[logratio.ij, "97.5%"], digits = digits), nsmall = digits), ")", sep = "")
          logratio.quan[j,i] <- paste(format(round(-smry[logratio.ij, "50%"], digits = digits), nsmall = digits), " (", format(round(-smry[logratio.ij, "97.5%"], digits = digits), nsmall = digits),
            ", ", format(round(-smry[logratio.ij, "2.5%"], digits = digits), nsmall = digits), ")", sep = "")
        }
      }
    }
    out$LogRateRatio <- list(Mean_SD = noquote(logratio.stat), Median_CI = noquote(logratio.quan))
  }

  if(is.element("rank.prob", param)){
    rank.prob.id <- grep("rank.prob", rownames(smry))
    rank.prob.stat <- array(format(round(smry[rank.prob.id, "Mean"], digits = 4), nsmall = 4), dim = c(ntrt, ntrt))
    colnames(rank.prob.stat) <- paste("rank", 1:ntrt, sep = "")
    rownames(rank.prob.stat) <- trtname
    out$TrtRankProb <- noquote(rank.prob.stat)
  }

  if(conv.diag){
    cat("Start calculating MCMC convergence diagnostic statistics...\n")
    conv.out <- gelman.diag(jags.out, multivariate = FALSE)
    conv.out <- conv.out$psrf
    if(is.element("rank.prob", param)){
      rank.prob.id <- grep("rank.prob", rownames(conv.out))
      conv.out <- conv.out[-rank.prob.id,]
    }
    write.table(conv.out, "ConvergenceDiagnostic.txt", row.names = rownames(conv.out), col.names = TRUE)
  }

  if(dic){
    cat("Start calculating deviance information criterion statistics...\n")
    dic.out <- dic.samples(model = jags.m, n.iter = n.iter, thin = n.thin)
    dev <- sum(dic.out$deviance)
    pen <- sum(dic.out$penalty)
    pen.dev <- dev + pen
    dic.stat <- rbind(dev, pen, pen.dev)
    rownames(dic.stat) <- c("D.bar", "pD", "DIC")
    colnames(dic.stat) <- ""
    out$DIC <- dic.stat
  }

  if(mcmc.samples){
    out$mcmc.samples <- jags.out
  }

  if(!is.null(trace)){
    cat("Start saving trace plots...\n")
  }

  if(is.element("rate", trace)){
    for(i in 1:ntrt){
      png(paste("TracePlot_rate_", trtname[i], ".png", sep = ""), res = 600, height = 8.5, width = 11, units = "in")
      par(mfrow = c(n.chains, 1))
      for(j in 1:n.chains){
        temp <- as.vector(jags.out[[j]][,paste("rate[", i, "]", sep = "")])
        plot(temp, type = "l", col = "red", ylab = "Rate", xlab = "Iteration", main = paste("Chain", j))
      }
      dev.off()
    }
  }
  if(is.element("lograte", trace)){
    for(i in 1:ntrt){
      png(paste("TracePlot_lograte_", trtname[i], ".png", sep = ""), res = 600, height = 8.5, width = 11, units = "in")
      par(mfrow = c(n.chains, 1))
      for(j in 1:n.chains){
        temp <- as.vector(jags.out[[j]][,paste("lograte[", i, "]", sep = "")])
        plot(temp, type = "l", col = "red", ylab = "Log Rate", xlab = "Iteration", main = paste("Chain", j))
      }
      dev.off()
    }
  }
  if(is.element("ratio", trace)){
    for(i in 1:ntrt){
      for(k in 1:ntrt){
        if(i < k){
          png(paste("TracePlot_ratio_", trtname[i], "_", trtname[k], ".png", sep = ""), res = 600, height = 8.5, width = 11, units = "in")
          par(mfrow = c(n.chains, 1))
          for(j in 1:n.chains){
            temp <- as.vector(jags.out[[j]][,paste("ratio[", i, ",", k, "]", sep = "")])
            plot(temp, type = "l", col = "red", ylab = "Rate Ratio", xlab = "Iteration", main = paste("Chain", j))
          }
          dev.off()
        }
      }
    }
  }
  if(is.element("logratio", trace)){
    for(i in 1:ntrt){
      for(k in 1:ntrt){
        if(i < k){
          png(paste("TracePlot_logratio_", trtname[i], "_", trtname[k], ".png", sep = ""), res = 600, height = 8.5, width = 11, units = "in")
          par(mfrow = c(n.chains, 1))
          for(j in 1:n.chains){
            temp <- as.vector(jags.out[[j]][,paste("logratio[", i, ",", k, "]", sep = "")])
            plot(temp, type = "l", col = "red", ylab = "Log Rate Ratio", xlab = "Iteration", main = paste("Chain", j))
          }
          dev.off()
        }
      }
    }
  }

  if(postdens){
    cat("Start saving posterior density plot for log rates...\n")
    mcmc <- NULL
    dens <- matrix(0, ntrt, 3)
    colnames(dens) <- c("ymax", "xmin", "xmax")
    for(i in 1:ntrt){
      temp <- NULL
      for(j in 1:n.chains){
        temp <- c(temp, as.vector(jags.out[[j]][,paste("lograte[", i, "]", sep = "")]))
      }
      mcmc[[i]] <- temp
      tempdens <- density(temp)
      dens[i,] <- c(max(tempdens$y), quantile(temp, 0.001), quantile(temp, 0.999))
    }
    ymax <- max(dens[,"ymax"])
    xmin <- min(dens[,"xmin"])
    xmax <- max(dens[,"xmax"])
    cols <- rainbow(ntrt, s = 1, v = 0.6)
    pdf("LogRateDensityPlot.pdf")
    par(mfrow = c(1, 1), mar = c(5.5, 5.5, 2, 2) + 0.1)
    plot(density(mcmc[[1]]), xlim = c(xmin, xmax), ylim = c(0, ymax), xlab = "Log Rate", ylab = "Density", main = "", col = cols[1], lty = 1, lwd = 2, cex.axis = 2, cex.lab = 2)
    for(i in 2:ntrt){
      lines(density(mcmc[[i]]), col = cols[i], lty = i, lwd = 2)
    }
    legend("topright", legend = trtname, col = cols, lty = 1:ntrt, lwd = 2, cex = 1.5)
    dev.off()
  }
  class(out) <- "nma.ab"
  return(out)
  options(warn = 0)
}