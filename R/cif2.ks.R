### Kolmogorov--Smirnov two-sample test of equality of cumulative incidence functions
### (Lin, 1997)

cif2.ks = function(x, group, cause=1, nsim=2000, nsim.plot=50)
{
	if (!is.Survcomp(x))
		stop("x must be a Survcomp object")
	
	cause = as.integer(cause)
	if (is.na(match(as.integer(cause),1:attr(x,"ncauses"))))
		stop("incorrect cause")
	
	if (missing(group))
		stop("group missing")
	n = nrow(x)
	if (n!=length(group)) stop("lengths of 'x' and 'group' differ")
	
	if (missing(group))
		group = rep(1,n)
	if ( !identical(sort(as.integer(unique(group))),1:2) )
		stop("values of group must be 1,2")
	if ((!is.numeric(nsim))||(nsim<=0))
		stop("incorrect nsim (must be >0)")
	if ((!is.numeric(nsim.plot))||(nsim<0)||(nsim.plot>nsim))
		stop("incorrect nsim.plot (must be 0<=nsim.plot<=nsim)")
	
	tim = x[,1]
	r = order(tim)
	tim = tim[r]
	evt = x[r,2]
	evt = match(evt,c(0,cause),nomatch=3)-1 # replaces in evt 0 by 0, cause by 1, other by 2 (C routine tests cause 1)
	grp = group[r]
	
	temp = .C("twosample_incidence_ks",
		as.integer(evt),
		as.integer(grp),
		as.integer(n),
		as.integer(nsim),
		f11 = double(n),
		f12 = double(n),
		f21 = double(n),
		f22 = double(n),
		test.process = double(n),
		stat = double(1),
		pval.sim = double(1),
		test.process.sim = matrix(0,n,max(1,nsim.plot)),
		as.integer(nsim.plot),
		PACKAGE = "surv2sample")
	
	out = list()
	out$cause = cause
	out$stat = temp$stat
	out$time = tim
	out$test.process = temp$test.process
	out$nsim = nsim
	out$pval.sim = temp$pval.sim
	out$nsim.plot = nsim.plot
	out$test.process.sim = temp$test.process.sim
	
	class(out) = c("cif2.ks","lwy.test")
	out
}

print.cif2.ks = function(x,...)
{
	if (!inherits(x,"cif2.ks"))
		stop("must be an object of class 'cif2.ks'")
	
	cat("\nTwo-sample test of equality of cumulative incidence functions")
	cat("\nKolmogorov--Smirnov test")
	cat("\nTested cause of failure: ",x$cause,sep="")
	cat("\nStatistic: ",x$stat,",   p-value: ",format.pval(x$pval.sim),sep="")
	cat("\np-value based on",x$nsim,"simulations")
	cat("\n\n")
	
	invisible()
}

