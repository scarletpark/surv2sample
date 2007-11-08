### Smooth test of equality of cumulative incidence functions in two samples
### with fixed dimension or data-driven with nested subsets

cif2.neyman = function(x, group, cause=1, data.driven = FALSE,
	d=ifelse(data.driven,5,3), basis="legendre", choltol=1e-07)
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
	
	if (!(d>0)) stop("d must be >0")
	basis = tolower(basis)
	if (basis=="cosine") basis = "cos"
	bas = match(tolower(basis),c("legendre","cos"))
	if (is.na(bas)) stop("unknown basis")
	
	d0=ifelse(data.driven,1,d)
	penalty = log(n)
	
	d.star = max(d0,1) # minimal dimension of subset
	
	temp = .C("twosample_incidence_neyman",
		as.double(tim),
		as.integer(evt),
		as.integer(grp),
		as.integer(n),
		as.integer(d),
		as.integer(d0),
		as.integer(bas),
		as.double(penalty),
		as.integer(0), # 0 = not Gray's test
		as.double(0), # dummy argument (r in Gray's G^r)
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
	
	if (data.driven) {
		out = list(stat = temp$stat.nested,
			pval = temp$pval.nested.twoterm,
			pval.twoterm = temp$pval.nested.twoterm,
			pval.asympt = temp$pval.nested.asympt,
			cause = cause,
			basis = basis,
			data.driven = data.driven,
			subsets.class = "nested",
			approx = "twoterm",
			d = d,
			d0 = d0,
			d.star = d.star,
			df = d.star,
			S.ind = temp$S+1, # index of the selected set (+1 because C returns 0,1,...)
			S.dim = d.star+temp$S, # the smallest set (with 1 function) has S=0 in C
			S.set = 1:(d.star+temp$S),
			stats = temp$stats.nested,
			stats.penal = temp$stats.nested.penal)
	} else {
		out = list(stat = temp$stat.d,
			pval = temp$pval.d,
			pval.asympt = temp$pval.d,
			cause = cause,
			basis = basis,
			data.driven = data.driven,
			approx = "asympt",
			d = d,
			df = d)
	}
	
	class(out) = c("cif2.neyman","neyman.test")
	out
}

print.cif2.neyman = function(x,...)
{
	if (!inherits(x,"cif2.neyman"))
		stop("must be an object of class 'cif2.neyman'")
	
	cat("\nTwo-sample test of equality of cumulative incidence functions")
	cat("\nTested cause of failure: ",x$cause,sep="")
	NextMethod("print",x,...)
	
	invisible()
}

summary.cif2.neyman = function(object,...)
{
	if (!inherits(object,"cif2.neyman"))
		stop("must be an object of class 'cif2.neyman'")
	
	cat("\nTwo-sample test of equality of cumulative incidence functions")
	cat("\nTested cause of failure: ",object$cause,sep="")
	NextMethod("summary",object,...)
	
	invisible()
}
