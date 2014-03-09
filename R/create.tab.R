create.tab <-
function(summary.stat,type="CI",compare=1,f.name="",trtname,order=TRUE,o.path=getwd()){
  ## check arguments
  if(type!="CI" & type!="SE") stop("type must be set as \"CI\" or \"SE\".")
  if(compare!=1 & compare!=2 & compare!=3 & compare!=4 & compare!=5 & compare!=6) stop("compare must be set as 1, 2, 3, 4, 5, or 6.")
    # 1 for upper RR lower RD
    # 2 for upper RD lower RR
    # 3 for upper RR lower OR
    # 4 for upper OR lower RR
    # 5 for upper RD lower OR
    # 6 for upper OR lower RD
  
  ## set variables
  cm<-matrix(c("RR","RD","RD","RR","RR","OR","OR","RR","RD","OR","OR","RD"),6,2,byrow=TRUE)
  r.probt<-summary.stat[grep("probt",row.names(summary.stat)),]
  r.RR<-summary.stat[grep("RR",row.names(summary.stat)),]
  r.RD<-summary.stat[grep("RD",row.names(summary.stat)),]
  r.OR<-summary.stat[grep("OR",row.names(summary.stat)),]
  n.tr<-nrow(r.probt)
  o.table<-matrix(0,nrow=n.tr,ncol=n.tr)
  probtname<-rownames(r.probt)
  rrname<-rownames(r.RR)
  rdname<-rownames(r.RD)
  orname<-rownames(r.OR)
  nm<-list(list(rrname,rdname),list(rdname,rrname),list(rrname,orname),list(orname,rrname),list(rdname,orname),list(orname,rdname))
  if(type=="CI"){
    t.probt<-sprintf(paste("%.2f(%.2f,%.2f)",sep=""),r.probt[,1],r.probt[,3],r.probt[,7])
    t.RR<-sprintf(paste("%.2f(%.2f,%.2f)",sep=""),r.RR[,1],r.RR[,3],r.RR[,7])
    t.RD<-sprintf(paste("%.2f(%.2f,%.2f)",sep=""),r.RD[,1],r.RD[,3],r.RD[,7])
    t.OR<-sprintf(paste("%.2f(%.2f,%.2f)",sep=""),r.OR[,1],r.OR[,3],r.OR[,7])
  }
  if(type=="SE"){
    t.probt<-sprintf(paste("%.2f(%.2f)",sep=""),r.probt[,1],r.probt[,2])
    t.RR<-sprintf(paste("%.2f(%.2f)",sep=""),r.RR[,1],r.RR[,2])
    t.RD<-sprintf(paste("%.2f(%.2f)",sep=""),r.RD[,1],r.RD[,2])
    t.OR<-sprintf(paste("%.2f(%.2f)",sep=""),r.OR[,1],r.OR[,2])
  }
  sm<-list(list(t.RR,t.RD),list(t.RD,t.RR),list(t.RR,t.OR),list(t.OR,t.RR),list(t.RD,t.OR),list(t.OR,t.RD))

  ## create table
  setwd(o.path)
  trtids<-row.names(summary.stat)[grep("probt",row.names(summary.stat))]
  trtids<-gsub("probt\\[","",trtids)
  trtids<-as.numeric(gsub("\\]","",trtids))
  orders<-order(trtids)
  sorted<-sort(trtids)
  if(missing(trtname)){
    trtname<-paste("treat",sorted,sep="")}else{
      if(length(trtname)!=n.tr) stop("the length of trtname is not equal to treatment number.")
    }
  for(i in 1:n.tr){
    for(j in 1:n.tr){
      if(i==j) o.table[i,j]<-t.probt[orders[i]]
      if(i<j) o.table[i,j]<-sm[[compare]][[1]][which(nm[[compare]][[1]]==paste(cm[compare,1],"[",ifelse(order,min(c(sorted[i],sorted[j])),max(c(sorted[i],sorted[j]))),",",ifelse(order,max(c(sorted[i],sorted[j])),min(c(sorted[i],sorted[j]))),"]",sep=""))]
      if(i>j) o.table[i,j]<-sm[[compare]][[2]][which(nm[[compare]][[2]]==paste(cm[compare,2],"[",ifelse(order,min(c(sorted[i],sorted[j])),max(c(sorted[i],sorted[j]))),",",ifelse(order,max(c(sorted[i],sorted[j])),min(c(sorted[i],sorted[j]))),"]",sep=""))]
    }
  }
  write.table(o.table,file=paste(f.name,cm[compare,1],"-",cm[compare,2],"-",type,"-table.txt",sep=""),row.names=trtname,col.names=trtname)
}
