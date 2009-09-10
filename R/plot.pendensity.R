#plotting the estimated density/densities, using the estimated density. Here, it's called 'obj'.
plot.pendensity <- function(x,plot.val=1,latt=FALSE,kernel=FALSE,confi=TRUE,main=NULL,sub=NULL,xlab=NULL,ylab=NULL,plot.base=FALSE,lwd=NULL,legend.txt=NULL,...) {
  library(lattice)
  library(fda)
  obj <- x
  if(plot.val==1) {
    weight <- obj$results$ck
    q <- obj$splines$q
    x.factor <- obj$values$covariate$x.factor
    help.env <- new.env()
    if(!is.null(x.factor)) {
      len.x.fac <- length(x.factor[1,])
      all.x <- 0
      all.x2 <- len.x.fac
    }
    max.bsp <-0
    max.kern <-0
    base <- obj$splines$base
    base.den <- obj$splines$base.den
    MeanW <- obj$splines$MeanW
    knots.spline <- obj$splines$knots.val$all
    knots.val <- obj$splines$knots.val
    Stand.abw <- obj$splines$Stand.abw
    help.degree <- obj$splines$help.degree
    K <- obj$splines$K
    N <- obj$splines$N
    Z <- obj$values$Z
    x <- obj$values$x
    m <- obj$splines$m
    Dm <- obj$splines$Dm
    h <- obj$splines$h
    sort <- obj$values$sort
    my.AIC <- obj$results$AIC$my.AIC
    lambda0 <- obj$results$lambda0
    var.par <- obj$results$variance.par
    eps <- 1e-4
    levels <- obj$values$covariate$levels
    how.combi <- obj$values$covariate$how.combi
    how.levels <- obj$values$covariate$how.levels
    
    lev1 <- c()
    for(i in 1:how.levels) lev1 <- c(lev1,as.numeric(as.vector(levels[[i]])))
    lev <- lev1[-1]
    if(base=="bspline") {
      h.help <- abs(knots.spline[1]-knots.spline[2])
      m <- length(knots.spline)
    }
    if(base=="gaussian") {
      h.help <- abs(MeanW[1]-MeanW[2])
      m <- length(MeanW)
    }

    max.bw.all <- c()
    
    if (is.null(x)) {
      y <- obj$values$y
      if(!sort) {
        o <- order(y)
        base.den <- base.den[,o]
        y <- y[o]
      }
      y.r <- range(knots.val$val)
      y.help <- seq(y.r[1],y.r[2],length=200)
      list <- c()
      
      if(base=="bspline") {
        if(q>2) K.help <- K-q+2 else K.help <- 0
        base.val<- my.bspline(h,q,knots.val,y.help,K.help,plot.bsp=FALSE)
        help.base.den <- base.val$base.den
        assign(paste("base1"),help.base.den,env=help.env)
        assign("y.help1",y.help,env=help.env)
      }
      if(base=="gaussian") {
        for(i in 1:m) {
          if(MeanW[i] >=y.r[1] & MeanW[i]<=y.r[2]) list <- c(list,i)
        }
        nn <- matrix(1:length(y.help))
        help.base.den <- apply(nn,1,function(i,MeanW,Stand.abw,y.help) dnorm(y.help[i],MeanW,Stand.abw),MeanW,Stand.abw,y.help)[list,]
        assign(paste("y.help1",sep=""),y.help,env=help.env)
        assign(paste("base1"),help.base.den,env=help.env)
      }   
      sd.cal <- new.env()
      sd.cal1 <- variance.val(help.base.den,var.par,weight,K,x,Z,x.factor) #confidence intervals
      assign("sd.cal1",sd.cal1,sd.cal)
      
      assign("weight1",weight,env=help.env)
      assign("y.later1",y,env=help.env)
      assign("len.x.fac1",1,env=help.env)
      assign("b.w1",bw <- get("weight1",env=help.env)%*%get("base1",env=help.env),help.env)

      max.bw <- max(bw)
      cut <- 0.0005*max.bw

      y.help <- get("y.help1",help.env)
      ind <- which(bw<cut)
      assign("ind1",ind,help.env)

      if(length(ind)!=0){
        assign("b.w1",bw[-ind],env=help.env)
        assign("y.help1",get("y.help1",help.env)[-ind],help.env)
      }
      max.bsp <- max.bw
      if(kernel) assign("kern.den1",density(get("y.later1",env=help.env),kernel="gaussian",bw="ucv"),env=help.env)
      if(kernel) max.kern <- max(get("kern.den1",env=help.env)$y) else max.kern <- 0
      all.x2 <- 1
      list.len <- 1
    }
    
    if (!is.null(x)) {
      y <- obj$values$y
      help.knots <- c()
      list <- c()
      list.len <- length(x.factor[,1])
      for(i in 1:list.len) {
        name <- paste("weight",i,sep="")
        obj2 <- weight[i,]
        assign(name,obj2,env=help.env)
      }
      if(!sort) {
        o <- order(y)
        y.help <- y[o]
        base.den <- base.den[,o]
        Z <- Z[o,]
      }
      else y.help <- y
      
      x.factor.len <- length(x.factor)
      for(j in 1:list.len) {
        set <- c()
        y <- c()
        for (i in 1:length(obj$values$y)) {
          if (all.equal(as.vector(Z[i,]),as.vector(x.factor[j,]))==TRUE) {
            y <- c(y,y.help[i])
            set <- c(set,i)
          }
        }
        assign(paste("y.later",j,sep=""),y,env=help.env)
      }
   
      for(j in 1:list.len) { #which covariate?
        knot.left <- 1
        knot.right <- K-help.degree
        
        y.help <- seq(knots.val$val[knot.left],knots.val$val[knot.right],length=200)
        assign(paste("y.help",j,sep=""),y.help,env=help.env)
        if(base=="bspline") {
          if(q>2) K.help <- K-q+2 else K.help <- 0
          base.val<- my.bspline(h,q,knots.val,y.help,K.help,plot.bsp=FALSE)
          help.base.den <- base.val$base.den
          assign(paste("base",j,sep=""),help.base.den,env=help.env)
        }
 
        if(base=="gaussian") {
          for(i in 1:m) if(MeanW[i] >=y.r[1] & MeanW[i]<=y.r[2]) list <- c(list,i)
          nn <- matrix(1:length(y.help))
          help.base.den <- apply(nn,1,function(i,k,MeanW,Stand.abw,y.help) dnorm(y.help[i],MeanW,Stand.abw),MeanW,Stand.abw,y.help)[list,]
          assign(paste("base",j,sep=""),help.base.den,env=help.env)
        }
        
        bw <- get(paste("weight",j,sep=""),env=help.env)%*%get(paste("base",j,sep=""),env=help.env)

        max.bw.all <- c(max.bw.all,max(bw))
        assign(paste("b.w",j,sep=""),bw,env=help.env)
        
        if(plot.base) max.bsp <- max(max(get(paste("b.w",j,sep=""),env=help.env)),max.bsp)
        if(kernel) assign(paste("kern.den",j,sep=""),density(get(paste("y.later",j,sep=""),env=help.env),kernel="epanechnikov"),env=help.env)
        if(kernel) max.kern <- max(max(get(paste("kern.den",j,sep=""),env=help.env)$y),max.kern) else max.kern <- 0
        
      }
      sd.cal <- variance.val(help.env,var.par,weight,K,x,list.len,Z,x.factor)

      max.bw <- max(max.bw.all)
      cut <- 0.0005*max.bw

      max.x.all <- c()
      min.x.all <- c()

      for(i in 1:list.len) {
        bw <- get(paste("b.w",i,sep=""),help.env)
        y.help <- get(paste("y.help",i,sep=""),help.env)
        ind <- which(bw<cut)
        if(length(ind)!=0) {
          assign(paste("b.w",i,sep=""),bw[-ind],help.env)
          y.help <- y.help[-ind]
        }
        assign(paste("y.help",i,sep=""),y.help,help.env)
        assign(paste("ind",i,sep=""),ind,help.env)

        max.x.all <- c(max(y.help),max.x.all)
        min.x.all <- c(min(y.help),min.x.all)
      }

      max.x <- max(max.x.all)
      min.x <- min(min.x.all)
    }
    AIC <- round(my.AIC,2)
    
    if(is.null(lambda0)) lambda0 <- 0
    lam <- round(lambda0,2)
    help.lam <- substitute(lambda ==s,list(s=lam))
    
    if(is.null(x)) {
      ind <- get("ind1",help.env)
      if(length(ind)!=0) {
        assign("conf.plus1",get("b.w1",help.env)+2*get("sd.cal1",sd.cal)[-ind],help.env)
        assign("conf.minus1",get("b.w1",help.env)-2*get("sd.cal1",sd.cal)[-ind],help.env)
      }
      else {
        assign("conf.plus1",get("b.w1",help.env)+2*get("sd.cal1",sd.cal),help.env)
        assign("conf.minus1",get("b.w1",help.env)-2*get("sd.cal1",sd.cal),help.env)
      }
        
      list.len <- 1
    }
    else {
      list.len <- length(x.factor[,1])
      c1 <- c()
      for(i in 1:list.len) {#which grouping of covariates?
        ind <- get(paste("ind",i,sep=""),help.env)
        if(length(ind)!=0) {
          assign(paste("conf.plus",i,sep=""),conf.plus <- get(paste("b.w",i,sep=""),help.env)+2*get(paste("sd.cal",i,sep=""),sd.cal)[-ind])
          c1 <- c(c1,conf.plus)
          assign(paste("conf.minus",i,sep=""),get(paste("b.w",i,sep=""),help.env)-2*get(paste("sd.cal",i,sep=""),sd.cal)[-ind])
        }
        else {
          assign(paste("conf.plus",i,sep=""),conf.plus <- get(paste("b.w",i,sep=""),help.env)+2*get(paste("sd.cal",i,sep=""),sd.cal))
          c1 <- c(c1,conf.plus)
          assign(paste("conf.minus",i,sep=""),get(paste("b.w",i,sep=""),help.env)-2*get(paste("sd.cal",i,sep=""),sd.cal))
        } 
      }
    }

    if(is.null(main)) main.title <- substitute("K= "*a*", AIC= "*b*", "*c*"="*d, list(a=K,b=AIC,c=parse(text="lambda")[[1]],d=lam))
    if(is.null(sub) & base=="bspline") sub.title <- paste("base is", base, " order is", m) else sub.title <- sub
    if(is.null(sub) & base=="gaussian") sub.title <- paste("base is", base) else sub.title <- sub
    if(is.null(xlab)) xlab.title <- "y" else xlab.title <- xlab
    if(is.null(ylab)) ylab.title <- "density" else ylab.title <- ylab
    if(is.null(lwd)) lwd.value <- 3 else lwd.value <- lwd
    if(!is.null(legend.txt) & length(legend.txt)!=list.len) stop("Input length not equal to length of groupings")
    if(!is.null(legend.txt) & kernel) legend.txt <- c(legend.txt,paste("kernel density estimation of ", legend.txt,sep=""))

    if(is.null(x)) x.range <- range(get("y.help1",help.env))
    else  x.range <- c(min.x,max.x)
    x.range.add <- 0.05*diff(x.range)
    x.lim <- c(x.range[1]-x.range.add,x.range[2]+x.range.add)

    if(plot.val==1 & latt==FALSE) {
      y.temp <- obj$values$y
      if(is.null(x)) plot(y.temp, dnorm(y.temp,mean=0,sd=1),xlab=xlab.title,ylab=ylab.title,xlim=x.lim,ylim=c(0,max(max(get("conf.plus1",help.env)),max.kern,max.bsp)),type="n", main=main.title, sub=sub.title,cex.axis=1.2,cex.lab=1.5,cex.main=1.8)
      
      if(!is.null(x)) plot(y.temp, dnorm(y.temp,mean=0,sd=1),xlab=xlab.title,ylab=ylab.title,xlim=x.lim,ylim=c(0,max(c1,max.kern,max.bsp)),type="n", main=main.title, sub=sub.title,cex.axis=1.2,cex.lab=1.5,cex.main=1.8)
      
      if(is.null(x)) {
        lines(get(paste("y.help1",sep=""),env=help.env),get(paste("b.w1",sep=""),env=help.env),type="l",lwd=lwd.value)
        if(confi)lines(get("y.help1",help.env),get("conf.plus1",help.env),type="l",lwd=lwd.value-1,lty=3)
        if(confi)lines(get("y.help1",help.env),get("conf.minus1",help.env),type="l",lwd=lwd.value-1,lty=3)
        if(kernel) lines(get(paste("kern.den1",sep=""),env=help.env),col=3,lwd=lwd.value,lty=2)
        if(!is.null(legend.txt)&!kernel) legend("topright",legend.txt,lty=1,lwd=lwd.value,cex=1.2)
        if(!is.null(legend.txt)&kernel) legend("topright",legend.txt,lty=c(1,2),col=c(1,3),lwd=lwd.value,cex=1.2)
        if(plot.base) {
          help <- get(paste("base",i,sep=""),env=help.env)
          help.w <- get(paste("weight",i,sep=""),env=help.env)
          help.l <- length(help[,1])
          for(j in 1:help.l) lines(get(paste("y.help",i,sep=""),env=help.env),help[j,]*help.w[j],col=1)
        }
      }
      else {
        for(i in 1:list.len) {#which grouping of covariates?
          if(kernel) lines(get(paste("kern.den",i,sep=""),env=help.env),col=2,lwd=lwd.value,lty=1+i-1)
          conf.plus <- get(paste("conf.plus",i,sep=""),help.env)
          conf.minus <- get(paste("conf.minus",i,sep=""),help.env)
          lines(get(paste("y.help",i,sep=""),env=help.env),t(get(paste("b.w",i,sep=""),env=help.env)),type="l",col=1,lwd=lwd.value,lty=1+i-1)
          if(confi) lines(get(paste("y.help",i,sep=""),help.env),conf.plus,type="l",col=1,lty=1+i-1,lwd=lwd.value-1)
          if(confi) lines(get(paste("y.help",i,sep=""),help.env),conf.minus,type="l",col=1,lty=1+i-1,lwd=lwd.value-1)
          if(!is.null(legend.txt)&!kernel) legend("topright",legend.txt,lty=seq(1,list.len),lwd=lwd.value,cex=1.2)
          if(!is.null(legend.txt)&kernel) legend("topright",legend.txt,lty=c(kronecker(seq(1,1+list.len-1),matrix(1,list.len,1))),col=c(kronecker(matrix(1,1,list.len),seq(1,1+list.len-1))),lwd=lwd.value)
          if(plot.base) {
            help <- get(paste("base",i,sep=""),env=help.env)
            help.w <- get(paste("weight",i,sep=""),env=help.env)
            help.l <- length(help[,1])
            for(j in 1:help.l) lines(get(paste("y.help",i,sep=""),env=help.env),help[j,]*help.w[j],col=1)
          }
        }
      }
    }
    
    if(plot.val==1 & latt==TRUE) {
      y <- c()
      fy <- c()
      fac <- c()
      confi.plus <- c()
      confi.minus <- c()
      for(i in 1:list.len) {
        y <- c(y,get(paste("y.help",i,sep=""),env=help.env))
        fy <- c(fy,get(paste("b.w",i,sep=""),env=help.env))
        if(confi) {
          y <- c(y,rep(get(paste("y.help",i,sep=""),help.env),2))
          fy <- c(fy,get(paste("conf.plus",i,sep=""),help.env))
          fy <- c(fy,get(paste("conf.minus",i,sep=""),help.env))
        }
        if(!is.null(x)) {
          fac <- c(fac,rep(paste(i,".x=",paste(x.factor[i,],collapse=",",sep=""),sep=""),length(get(paste("b.w",i,sep=""),env=help.env))))
          if(confi) fac <- c(fac,rep(paste(i,".confi+ x=",paste(x.factor[i,],collapse=",",sep=""),sep=""),length(get(paste("conf.plus",i,sep=""),help.env))),rep(paste(i,".confi- x=",paste(x.factor[i,],collapse=",",sep=""),sep=""),length(get(paste("conf.plus",i,sep=""),help.env))))
        }
        if(is.null(x)) {
          fac <- c(fac,rep("y~1",length(get(paste("b.w",i,sep=""),env=help.env))))
          if(confi) fac <- c(fac,rep("confi+ y~1",length(get(paste("conf.plus",i,sep=""),help.env))),rep("confi- y~1",length(get(paste("conf.plus",i,sep=""),help.env))))
        }
      }
      graph.sets <-list(superpose.line=list(col=sort(rep(1:list.len,3))),superpose.symbol = list(col = sort(rep(1:list.len,3))))
      datafr <- data.frame(y,fy,fac)
      help55 <- xyplot(fy~y,groups=fac,type="l",auto.key=list(space="right",title="grouping",sort=FALSE),main=main.title,sub=sub.title,data=datafr,xlab=xlab.title,ylab=ylab.title,par.settings=graph.sets,lty=rep(c(2,2,1),list.len))
      print(help55)
    }
  }
  
if(plot.val==2) {
  sort <- obj$values$sort
  y.list <- list()
  sum.list <- list()
  if(is.null(obj$x)) {
    knots <- obj$splines$knots
    beta <- obj$results$beta.val
    x.factor <- obj$values$covariate$x.factor
    N <- obj$splines$N
    p <- obj$splines$p
    if(sort) {
      y <- obj$values$y
      base.den2 <- obj$splines$base.den2
    }
    else {
      y.order <- order(obj$values$y)
      y <- obj$values$y[y.order]
      base.den2 <- obj$splines$base.den2[,y.order]
    }
    K <- length(base.den2[,1])-1
    weight <- c(obj$results$ck)
    row.help <- rep(0,length(base.den2[1,]))
    sum <- c(0)
    for(k in 1:K) {
      sum <- weight[k]*colSums(base.den2[(k:(K+1)),]) +sum
    }
    y <- sort(y)
    y.list[[1]] <- y
    sum.list[[1]] <- sum
    plot(y,sum,xlab="y",ylab="F(y)",main="Distribution of f(y)")
    return(list(y=y.list,sum=sum.list))
  }
  if(!is.null(obj$values$x)) {
    base.den2 <- obj$splines$base.den2
    K <- length(base.den2[,1])-1
    help.env <- new.env()
    base <- obj$splines$base   
    knots.spline <- obj$splines$knots.spline$all
    beta <- obj$results$beta.val
    x.factor <- obj$values$covariate$x.factor
    N <- obj$splines$N
    y <- obj$values$y
    x <- obj$values$x
    len.x.fac <- length(x.factor[,1])
    all.x <- 0
    all.x2 <- 1
    Z <- obj$values$Z
    weight <- obj$results$ck

    for(i in 1:len.x.fac) {
        name <- paste("weight",i,sep="")
        assign(name,weight[i,],env=help.env)
      }
    y.help <- c()
    y.call <- c()
    if(!sort) {
      o <- order(y)
      y.help <- y[o]
      base.den2 <- base.den2[,o]
    }
    else {
      y.help <- y
    }
    x.factor.len <- length(x.factor[1,])
    for(j in 1:len.x.fac) {
      set <- c()
      y <- c()
      for (i in 1:length(obj$values$y)) {
          if (all.equal(as.vector(Z[i,]),as.vector(x.factor[j,]))==TRUE) {
          y <- c(y,y.help[i])
          set <- c(set,i)
         }
      }
      assign(paste("y.later",j,sep=""),y,env=help.env)
      assign(paste("y.list",j,sep=""),set,env=help.env)
    }
    par(mfrow=c(2,ceiling(len.x.fac/2)))
    for(j in 1:len.x.fac) { #which covariate?
      y.r <- range(get(paste("y.later",j,sep=""),env=help.env))
      y.list <- get(paste("y.list",j,sep=""),env=help.env)
      if(base=="bspline") {
        row.help <- rep(0,length(base.den2[1,]))
        base.den2 <- rbind(base.den2,row.help)
        sum <- c(0)
        weight <- get(paste("weight",j,sep=""),help.env)
	  for(k in 1:K) {
        	sum <- weight[k]*colSums(base.den2[(k:(K+1)),y.list]) +sum
        }
        assign(paste("distr",j,sep=""),sum,env=help.env)
        y.call[[j]] <- get(paste("y.later",j,sep=""),env=help.env)
        sum.list[[j]] <- sum
        plot(get(paste("y.later",j,sep=""),env=help.env),sum,xlab="y",ylab=paste("F(y|x",j,")",sep=""),main=paste("Distribution of f(y|x",j,")",sep=""))
      }
      }
  return(list(y=y.call,sum=sum.list))
  }
}
  if(plot.val==3) {
    help.env <- distr.func.help(obj)
    func <- distr.func(yi=NULL,obj,help.env)
    len.b <- length(obj$splines$base.den[,1])-obj$splines$q+1
    x.factor <- get("x.factor",env=func)
    if(!is.null(obj$values$x)) x.factor.len <- length(x.factor[1,]) else x.factor.len <- 1
    knots.val <- obj$splines$knots.val
    all.x2 <- get("allx",env=func)
    par(mfrow=c(ceiling(sqrt(all.x2)),1))
    y <- obj$values$y
    Z <- obj$values$Z
    x <- obj$values$x
    eps <- 1e-10
    y.help <- c()
    for(i in 1:all.x2) {
      w <- c()
      tt <- c()
      if(!is.null(x)) com.h <- x.factor[i,]
      
      if(!is.null(x)) {
        for(j in 1:length(Z[,1])) {
          if(identical(as.vector(Z[j,]),as.vector(com.h))==TRUE) y.help <- c(y.help,y[j])
        }
      }
      else y.help <- y
      
      min.y <- min(y.help)
      max.y <- max(y.help)
      val.min <- c()
      val.max <- c()
      for(k in 1:(length(knots.val$val)-1)) {
        if(knots.val$val[k]-eps <= min.y & min.y <= knots.val$val[k+1]+eps) val.min <- k
        if(knots.val$val[k]-eps <= max.y & max.y <= knots.val$val[k+1]+eps) val.max <- k
      }
      for(j in val.min:val.max) {
        funcy <- get(paste("distr.func",i,".",j,sep=""),env=func)
        eval(parse(text=funcy))
        if(j==val.min) xi <- seq(min.y,knots.val$val[j+1],length=100)
        if(j!=val.min | j!=val.max) xi <- seq(knots.val$val[j],knots.val$val[j+1],length=100)
        if(j==val.max) xi <- seq(knots.val$val[j],max.y,length=100)
        ti <- obj(xi)
        w <- c(w,xi)
        tt <- c(tt,ti)
      }
      assign(paste("x",i,sep=""),w,env=func)
      assign(paste("F(x)",i,sep=""),tt,env=func)
      plot(w,tt,xlab="y",ylab="F(y)",main=paste("Distribution function of f(y|x",i,")",sep=""))
    }
  }
}