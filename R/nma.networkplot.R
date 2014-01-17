nma.networkplot <-
function(c1,c2,percomparison,trtname,weight=FALSE,VAR1,graphtitle,thickness,nodetextsize,nodesize){
  ## check arguments
  if(missing(c1)) stop("need to specify the first argument.")
  if(missing(c2)) stop("need to specify the second argument.")
  if(length(c1)!=length(c2)) stop("the lengths of the first two arguments are not equal.")
  if(missing(percomparison)) stop("need to specify the type of the first two arguments: \n TRUE for both are treatments; FALSE for c1 is study id and c2 is treatment.")
  if(!missing(trtname)) {if(length(trtname)!=length(unique(c2))) stop("input trtname can not correspond to treatments.")}

  ## settings for compared treatments
  if(percomparison){
    t1<-c1
    t2<-c2}else{
      ID<-c1
      Treat<-c2
      treat<-as.double(as.factor(Treat))
      id<-ID
      nt<-length(unique(treat))
      ns<-length(unique(ID))
      na<-table(match(id,unique(id)))
      ## check whether the treatments are 1 to nt
      if(max(sort(unique(treat))-c(1:nt))>0){
        treat<-as.factor(treat)
        levels(treat)<-c(1:nt)
        treat<-as.double(treat)
      }
      TT<-matrix(NA,nrow=ns,ncol=max(na))
      for(i in 1:ns){
        TT[i,1:na[i]]<-treat[id==unique(id)[i]]
      }
      u<-TT
      new.id<-1:length(table(id))
      TT<-apply(TT,1,sort)
      na<-na
      u<-c()
      torepeat<-na[na>2]
      ## if there are more than 3 arms
      if(sum(na>2)){
        for(i in 1:sum(na>2)){
          u<-rbind(u,t(combn(unlist(TT[na>2][i]),2)))
        }
      }
      TT<-rbind(matrix(unlist(TT[na==2]),ncol=2,byrow=TRUE),u)
      t1<-TT[,1]
      t2<-TT[,2]
    }

  numoftreatments<-max(cbind(t1,t2))
  ## default adjacency-style matrix initialization
  mat_treat<-matrix(0,numoftreatments,numoftreatments)

  if(missing(trtname)){
    nam=""
    for(i in 1:numoftreatments){
      nam[i]<-paste("treat",i,sep=".")
      assign(nam[i],1:i)
    }
    colnames(mat_treat)<-nam
    rownames(mat_treat)<-nam
  }else{
    colnames(mat_treat)<-trtname
    rownames(mat_treat)<-trtname
  }

  ## based on the frequency of treatments we construct the node's thickness
  ## divisor has fixed value equal to 5 for aesthetical reasons
  if(missing(nodesize)){
    divisor<-5}else{
      divisor<-5/nodesize
    }

  if(missing(VAR1)){
    nodethickness<-table(c(t1,t2))/divisor}else{
      nodethickness<-VAR1
    }

  if(missing(thickness)) thickness<-10

  ## combined treatments
  tr<-cbind(t1,t2)

  for(i in 1:length(t1)){
    mat_treat[tr[i,1],tr[i,2]]<-mat_treat[tr[i,1],tr[i,2]]+1
  }
  mat_treatb<-t(mat_treat)
  mat_treat<-mat_treat+mat_treatb

  ## Based on package "network" we construct the netdata
  netdata<-network(mat_treat,directed=FALSE)

  ## network plotting
  mat_treat_new<-mat_treat*thickness/max(mat_treat)

  if(missing(nodetextsize)){
    par(cex=1)}else{
      par(cex=nodetextsize)
    }

  if(missing(graphtitle)){
    par(mai=c(0,0,0,0))
    par(cex=1)
    pltitle<-""
  }else{
      par(mai=c(0.2,0.2,0.2,0.2))
      par(cex=1)
      pltitle<-graphtitle
    }

  if(weight) wt_label<-mat_treat else wt_label<-NULL
  ## plot of the network
  plot(netdata,mode="circle",displaylabels=TRUE,vertex.cex=nodethickness,edge.col="skyblue2",boxed.label=FALSE,edge.lwd=mat_treat_new,edge.label=wt_label)
  title(main=list(pltitle,col="blue",font=2))
}
