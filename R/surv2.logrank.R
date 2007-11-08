### Weighted logrank two-sample tests of equality of survival distributions
### and combinations (sum, max) of weighted logrank tests

surv2.logrank = function(x, group, rho.gamma=c(0,0), comb, sum.weights,
	approx="perm", nsim=2000, choltol=1e-07)
{
	if (!is.Surv(x))
		stop("x must be a Surv object")
	if (attr(x,"type")!="right")
		stop("right censored data allowed only")
	if (missing(group))
		stop("group missing")
	n = nrow(x)
	if (n!=length(group)) stop("lengths of 'x' and 'group' differ")
	if ( !identical(sort(as.integer(unique(group))),1:2) )
		stop("values of group must be 1,2")
	
	if (is.vector(rho.gamma)&&(length(rho.gamma)==2)) {
		rho.gamma = list(rho.gamma)
	} else {
		if (!is.list(rho.gamma)) stop("'rho.gamma' must be a vector of length 2 or a list of vectors of length 2")
	}
	rho.gamma.mat = matrix(unlist(rho.gamma),ncol=2,byrow=TRUE) # each row one pair (rho,gamma)
	m = nrow(rho.gamma.mat)
	
	if (m==1)
		comb = "none"
	else {
		if (missing(comb)) {
			comb = "max"
		} else {
		comb = tolower(comb)
			if (is.na(match(comb,c("max","sum")))) stop("unknown combination method (must be 'max' or 'sum')")
		}
	}
	if ((comb=="sum")) {
		if (missing(sum.weights)) sum.weights = rep(1/m,m)
		if (length(sum.weights)!=m) stop("incorrect length of 'sum.weights'")
	} else {
		sum.weights = rep(1/m,m)
	}
	
	aprx = tolower(approx)
	if (is.na(match(aprx,c("asympt","perm","boot"))))
		stop("unknown approximation method (possible values 'asympt', 'perm', 'boot')")
	if ( ((m>1)||(aprx!="asympt")) && (nsim<=0) ) # for perm, boot, or for comb. stat nsim must be >0
		stop("incorrect nsim (nsim>0 needed)")
	
	tim = x[,1]
	r = order(tim)
	tim = tim[r]
	evt = x[r,2]
	grp = group[r]
	
	if (aprx=="asympt") {
		nsim.asympt = ifelse(m>1,nsim,0) # simulation from asympt. distribution of comb. stat.
		nperm = 0
		nboot = 0
	}
	if (aprx=="perm") {
		nperm = nsim
		nsim.asympt = 0
		nboot = 0
	}
	if (aprx=="boot") {
		nboot = nsim
		nsim.asympt = 0
		nperm = 0
	}
	
	temp = .C("twosample_grhogamma",
		as.double(tim),
		as.integer(evt),
		as.integer(grp),
		as.integer(n),
		as.double(rho.gamma.mat[,1]),
		as.double(rho.gamma.mat[,2]),
		as.double(sum.weights), # weights for the sum of statistics
		as.integer(m),
		as.integer(nsim.asympt),
		as.integer(nperm),
		as.integer(nboot),
		as.double(choltol),
		stats = double(m),
		sigma = matrix(double(1),m,m),
		pvals.asympt = double(m),
		pvals.perm = double(m),
		pvals.boot = double(m),
		stat.sum = double(1),
		pval.sum.asympt = double(1),
		pval.sum.perm = double(1),
		pval.sum.boot = double(1),
		stat.max = double(1),
		pval.max.asympt = double(1),
		pval.max.perm = double(1),
		pval.max.boot = double(1),
		PACKAGE = "surv2sample")
	
	if (m==1) { # weighted logrank
		out = list(stat = temp$stats[1],
			pval = temp[[paste("pvals",aprx,sep=".")]][1],
			rho.gamma = rho.gamma.mat[1,],
			approx = aprx)
		if (aprx!="asympt") out$nsim = nsim
	} else { # combination test
		out = list(stat = temp[[paste("stat",comb,sep=".")]],
			pval = temp[[paste("pval",comb,aprx,sep=".")]],
			stats = cbind(rho.gamma.mat,temp$stats,temp[[paste("pvals",aprx,sep=".")]]),
			comb = comb)
		colnames(out$stats) = c("rho","gamma","stat","pval")
		if (comb=="sum") out$sum.weights = sum.weights
		out$approx = aprx
		out$nsim = nsim
	}
	
	class(out) = "surv2.logrank"
	out
}

print.surv2.logrank = function(x,...)
{
	if (!inherits(x,"surv2.logrank"))
		stop("must be an object of class 'surv2.logrank'")
	
	if (is.null(x$comb)||(x$comb=="none")) {
		cat("\nTwo-sample weighted logrank test",sep="")
		cat("\nG(",x$rho.gamma[1],",",x$rho.gamma[2],") weight",sep="")
	} else {
		cat("\nCombination of two-sample weighted logrank tests",sep="")
		cat("\n",ifelse(x$com=="max","Maximum","Sum")," of ",
			paste("G(",x$stats[,1],",",x$stats[,2],")",sep="",collapse=", "),
			" statistics",sep="")
	}

	cat("\nTest statistic: ",x$stat,",   p-value: ",format.pval(x$pval),sep="")
	if (x$approx=="perm")
		cat("\np-value based on ",x$nsim," permutations",sep="")
	if (x$approx=="boot")
		cat("\np-value based on ",x$nsim," bootstrap samples",sep="")
	cat("\n\n")	

	invisible()
}
