
multiplot = function(x,y,types="s",lwds=1,cols=1,ltys=1,xlab="",ylab="",ylim=1.05*range(y,finite=TRUE),...)
{
	y = as.matrix(y)
	nc = ncol(y)
	wh = which(is.finite(y[,1]))
	plot(x[wh],y[wh,1],type=types[1],lwd=lwds[1],col=cols[1],lty=ltys[1],xlab=xlab,ylab=ylab,ylim=ylim,...)
	if (length(types)==1) types = rep(types,nc) else if(length(types)!=nc) stop("wrong length of 'types'")
	if (length(lwds)==1) lwds = rep(lwds,nc) else if(length(lwds)!=nc) stop("wrong length of 'lwds'")
	if (length(cols)==1) cols = rep(cols,nc) else if(length(cols)!=nc) stop("wrong length of 'cols'")
	if (length(ltys)==1) ltys = rep(ltys,nc) else if(length(ltys)!=nc) stop("wrong length of 'ltys'")
	if (nc>1) {
		for (j in 2:nc) {
			wh = which(is.finite(y[,j]))
			lines(x[wh],y[wh,j],type=types[j],lwd=lwds[j],col=cols[j],lty=ltys[j],...)
		}
	}
	invisible()
}

plot.lwy.test = function(x,lwds=c(3,1),cols=c("black","grey"),ltys=c(1,1),...)
{
	if (!inherits(x,"lwy.test"))
		stop("must be an object of class 'lwy.test'")
	
	xx = c(0,sort(x$time))
	yy = if (x$nsim.plot>0) rbind(rep(0,1+x$nsim.plot),cbind(x$test.process,x$test.process.sim)) else c(0,x$test.process)
	
	lwds=c(lwds[1],rep(lwds[2],x$nsim.plot))
	cols=c(cols[1],rep(cols[2],x$nsim.plot))
	ltys=c(ltys[1],rep(ltys[2],x$nsim.plot))
	
	multiplot(xx,yy,types="l",lwds=lwds,cols=cols,ltys=ltys,...)
	lines(xx,yy[,1],lwd=lwds[1],col=cols[1],lty=ltys[1],...)
	
	invisible()
}

