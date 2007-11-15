### Cumulative incidence functions

cif = function(x, group)
{
	if (!is.Survcomp(x))
		stop("x must be a Survcomp object")
	
	n = nrow(x)
	if (n!=length(group)) stop("lengths of 'x' and 'group' differ")
	
	if (missing(group))
		group = rep(1,n)
	K = max(group)
	if ( !identical(sort(as.integer(unique(group))),1:K) )
		stop("values of group must be integers 1,2,...")
	J = attr(x,"ncauses")
	
	out = list()
	
	for (k in 1:K) {
		group.k = which(group==k)
		n.k = length(group.k)
		event.k = x[group.k,2]
		time.k = x[group.k,1]
		r = order(time.k)
		event.k = event.k[r]
		km.k = c(1,cumprod(1-(event.k>0)/(n.k-(1:n.k)+1))[-n.k])
		out[[k]] = list()
		out[[k]]$time = time.k
		out[[k]]$surv = km.k
		f = matrix(0,n.k,J)
		for (i in 1:n.k) {
			for (j in 1:J) {
				f[i,j] = km.k[i]*(event.k[i]==j)/(n.k-i+1)
			}
		}
		out[[k]]$f = apply(f,2,cumsum)
	}
	
	out$time = x[,1]
	out$event = x[,2]
	out$group = group
	out$ncauses = J
	out$ngroups = K
	
	class(out) = "cif"
	out
}

## print method for the "cif" class

print.cif = function(x,...)
{
	if (!inherits(x,"cif"))
	stop("must be an object of class 'cif'")
	print(table(group=x$group,event=x$event))
	invisible()
}

## plot method for the "cif" class

plot.cif = function(x,by="group",aggreg.cif=TRUE,orient="land",lwds=1,cols=1,
	ltys=if ((by=="cause")||(by=="c")) rep(1:6,len=x$ngroups) else 1,xlab="",ylab="",
	ylim,mfrow,mfcol,mains,...)
{
	if (!inherits(x,"cif"))
		stop("must be an object of class 'cif'")
	
	if ((b<-match(by,c("group","g","cause","c"),nomatch=0))==0)
		stop("unknown value of the 'by' parameter")
	b = c(1,1,2,2)[b]
	
	np = ifelse(b==1, x$ngroups, x$ncauses)
	if (orient=="land") {
		mfr = if (np == 1) c(1, 1)
			else if (np <= 2) c(1, 2)
				else if (np <= 4) c(2, 2)
					else if (np <= 6) c(2, 3)
						else if (np <= 8) c(2, 4)
							else if (np <= 9) c(3, 3)
								else if (np <= 12) c(3, 4)
									else if (np <= 16) c(4, 4)
										else if (np <= 20) c(4, 5)
											else stop("number of plots too high")
	} else {
		mfr = if (np == 1) c(1, 1)
			else if (np <= 2) c(2, 1)
				else if (np <= 4) c(2, 2)
					else if (np <= 6) c(3, 2)
						else if (np <= 8) c(4, 2)
							else if (np <= 9) c(3, 3)
								else if (np <= 12) c(4, 3)
									else if (np <= 16) c(4, 4)
										else if (np <= 20) c(5, 4)
											else stop("number of plots too high")
	}

	if (!missing(mfrow))
		par(mfrow=mfrow)
	else
		if (!missing(mfcol))
			par(mfcol=mfcol)
		else
			par(mfrow=mfr)
	
	if (!missing(mains)) {
		if (length(mains)==1)
			mains = rep(mains,np) # replicate the same title for all plots
		else
			if (length(mains)!=np) stop("wrong length of 'mains'")
	}
	if (b==1) { # by group
		for (k in 1:x$ngroups) {
			xx = c(0,x[[k]]$time)
			yy = rbind(rep(0,x$ncauses),as.matrix(x[[k]]$f))
			if (aggreg.cif) {
				for (i in 1:nrow(yy)) yy[i,] = cumsum(yy[i,])
			}
			# ylim automatically determined for each plot, or all plots will have the same ylim
			if (missing(ylim)) yl = 1.07*range(yy) else yl = ylim
			multiplot(xx,yy,ylim=yl,main=ifelse(missing(mains),paste("Group",k),mains[k]),
				lwds=lwds,cols=cols,ltys=ltys,...)
		}
	} else { # by cause
		# times are different for each plot, multiplot can't be used
		if (length(lwds)==1) lwds = rep(lwds,x$ngroups) else if(length(lwds)!=x$ngroups) stop("wrong length of 'lwds'")
		if (length(cols)==1) cols = rep(cols,x$ngroups) else if(length(cols)!=x$ngroups) stop("wrong length of 'cols'")
		if (length(ltys)==1) ltys = rep(ltys,x$ngroups) else if(length(ltys)!=x$ngroups) stop("wrong length of 'ltys'")
		for (j in 1:x$ncauses) {
			ylim = c(0,0)
			for (k in 1:x$ngroups) ylim = range(c(ylim,x[[k]]$f[,j]))
			ylim = 1.07*ylim
			plot(x[[1]]$time,x[[1]]$f[,j],type="s",lwd=lwds[1],col=cols[1],lty=ltys[1],
				xlab=xlab,ylab=ylab,ylim=ylim,main=ifelse(missing(mains),paste("Cause",j),mains[j]),...)
			if (x$ngroups>1) {
				for (k in 2:x$ngroups)
					lines(x[[k]]$time,x[[k]]$f[,j],type="s",lwd=lwds[k],col=cols[k],lty=ltys[k],...)
			}
		}
	}
	
	par(mfrow=c(1,1))
	invisible()
}

