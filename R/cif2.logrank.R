### Gray's (1988) weighted _subdistribution_ logrank-type test of equality
### of cumulative incidence functions in two samples
### G^rho weights, i.e., (1-F0(t,cause))^rho

cif2.logrank = function(x, group, cause=1, rho=0)
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
	
	# the logrank-type test is computed by the same C-function as Neyman's test
	# we need some dummy parameters
	d = d0 = 1
	bas = 1
	penalty = log(n)
	choltol=1e-07
	
	# in the following output, values concerning Neyman's test meaningless, ignored further
	temp = .C("twosample_incidence_neyman",
		as.double(tim),
		as.integer(evt),
		as.integer(grp),
		as.integer(n),
		as.integer(d),
		as.integer(d0),
		as.integer(bas),
		as.double(penalty),
		as.integer(1), # 1 = Gray's test
		as.double(rho), # rho in the G^rho logrank-type weight
		as.double(choltol),
		score = double(d),
		sigma = matrix(double(d*d),d,d),
		f11 = double(n),
		f12 = double(n),
		f21 = double(n),
		f22 = double(n),
		stat.d = double(1),
		pval.d = double(1),
		stat.d0 = double(1),
		pval.d0 = double(1),
		stat.nested = double(1),
		stats.nested = double(d-d0+1),
		stats.nested.penal = double(d-d0+1),
		S = integer(1),
		pval.nested.asympt = double(1),
		pval.nested.twoterm = double(1),
		PACKAGE = "surv2sample")
	
	out = list(stat = temp$score/sqrt(temp$sigma[1,1]),
			pval = temp$pval.d,
			cause = cause,
			rho = rho)
			
	class(out) = "cif2.logrank"
	out
}

print.cif2.logrank = function(x,...)
{
	if (!inherits(x,"cif2.logrank"))
		stop("must be an object of class 'cif2.logrank'")
	
	cat("\nTwo-sample test of equality of cumulative incidence functions")
	cat("\nWeighted logrank test (G^rho weight with rho = ",x$rho,")",sep="")
	cat("\nTested cause of failure: ",x$cause,sep="")
	cat("\nTest statistic: ",x$stat,",   p-value: ",format.pval(x$pval),sep="")
	cat("\n\n")	

	invisible()
}
