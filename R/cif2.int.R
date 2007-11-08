### Integrated-difference two-sample test of equality of cumulative incidence functions
### (Pepe, 1991)

cif2.int = function(x, group, cause=1, tau, nsim=0)
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

	if ( !identical(sort(as.integer(unique(group))),1:2) )
		stop("values of group must be 1,2")
	
	tim = x[,1]
	r = order(tim)
	tim = tim[r]
	evt = x[r,2]
	evt = match(evt,c(0,cause),nomatch=3)-1 # replaces in evt 0 by 0, cause by 1, other by 2 (C routine tests cause 1)
	grp = group[r]
	
	if (missing(tau)) tau = tim[n] else if ((!is.numeric(tau))||(tau<=0)) stop("incorrect tau")
	if ((!is.numeric(nsim))||(nsim<0))
		stop("incorrect nsim (must be >=0)")
	
	temp = .C("twosample_incidence_pepe",
		as.double(tim),
		as.integer(evt),
		as.integer(grp),
		as.integer(n),
		as.double(tau),
		as.integer(nsim),
		f11 = double(n),
		f12 = double(n),
		f21 = double(n),
		f22 = double(n),
		stat. = double(1),
		pval.asympt = double(1),
		pval.sim = double(1),
		PACKAGE = "surv2sample")
	
	out = list()
	
	out$cause = cause
	out$tau = tau
	out$stat = temp$stat
	out$pval.asympt = temp$pval.asympt
	out$nsim = nsim
	out$pval.sim = temp$pval.sim
	
	class(out) = "cif2.int"
	out
}

print.cif2.int = function(x,...)
{
	if (!inherits(x,"cif2.int"))
		stop("must be an object of class 'cif2.int'")
	
	cat("\nTwo-sample test of equality of cumulative incidence functions")
	cat("\nPepe's integral test")
	cat("\nTested cause of failure: ",x$cause,",   tau = ",x$tau,sep="")
	if (x$nsim>0) {
		cat("\nStatistic: ",x$stat,",   p-value: ",format.pval(x$pval.sim),sep="")
		cat("\np-value based on",x$nsim,"simulations")
	} else {
		cat("\nStatistic: ",x$stat,",   p-value: ",format.pval(x$pval.asympt),sep="")
	}
	cat("\n\n")
	
	invisible()
}
