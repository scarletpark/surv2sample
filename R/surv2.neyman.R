### Smooth test of equality of survival distributions in two samples
### with fixed dimension, or data-driven with nested or all subsets

surv2.neyman = function(x, group, data.driven=FALSE, subsets="nested",
	d=ifelse(data.driven,5,3), d0=0, basis="legendre", time.transf="F",
	approx="perm", nsim=2000, choltol=1e-07)
{
	# d0 means that functions phi_1,...,phi_d0 are always included
	
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
	
	tim = x[,1]
	r = order(tim)
	tim = tim[r]
	evt = x[r,2]
	grp = group[r]
	
	time.transf=toupper(time.transf)
	if (time.transf=="L") time.transf="A"
	tim.trnsf = match(time.transf,c("F","A","I"))
	if (is.na(tim.trnsf)) stop("unknown time transformation")
	basis = tolower(basis)
	if (basis=="cosine") basis = "cos"
	bas = match(tolower(basis),c("legendre","cos"))
	if (is.na(bas)) stop("unknown basis")
	if (d<1) stop("d must be >= 1")
	
	aprx = tolower(approx)
	if (is.na(match(aprx,c("asympt","perm","boot"))))
		stop("unknown approximation method (possible values 'asympt', 'perm', 'boot')")
	if ( ((aprx=="perm")||(aprx=="boot")) && (nsim<=0) ) # for perm, boot, or for comb. stat nsim must be >0
		stop("incorrect nsim (nsim>0 needed)")
	
	if (data.driven) {
		if (d0>d) stop("d0 must be <= d")
		subsets = tolower(subsets)
		dta.drvn = match(subsets,c("nested","all"))
		if (is.na(dta.drvn)) stop("unknown class of subsets (possible values 'nested', 'all')")
		if (aprx=="asympt") {
			# asymptotics is OK for "all" with d0=0 or "nested" with d0<=1 (in which case it' replaced by two-term approx)
			if ( ((subsets=="nested")&&(d0>1)) || ((subsets=="all")&&(d0>0)) )
				warning("asymptotics inaccurate, permutations or bootstrap recommended")
		}
	} else {
		dta.drvn = 0
	}
	nalt.nested = d-d0+1 - (d0==0)
	nalt.all = 2^(d-d0) - (d0==0) # if d0=0, empty set is excluded
	
	if (aprx=="asympt") {
		nsim.asympt = ifelse(data.driven&&(subsets=="all")&&(d0==0),nsim,0) # simulation from asympt. distribution (max of dep. chisq1) with all subs, d0=0
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
	
	penalty = log(n)
	
	temp = .C("twosample_neyman",
		as.double(tim),
		as.integer(evt),
		as.integer(grp),
		as.integer(n),
		as.integer(d),
		as.integer(d0),
		as.integer(tim.trnsf),
		as.integer(bas),
		as.integer(dta.drvn), # 0=fixed alt., 1=nested, 2=all
		as.double(penalty),
		as.integer(0), # no empty set
		as.double(choltol),
		as.integer(nsim.asympt),
		as.integer(nperm),
		as.integer(nboot),
		score = double(d),
		sigma = matrix(double(1),d,d),
		stat.d = double(1),
		pval.d.asympt = double(1),
		pval.d.perm = double(1),
		pval.d.boot = double(1),
		stat.nested = double(1),
		pval.nested.asympt = double(1),
		pval.nested.twoterm = double(1),
		pval.nested.perm = double(1),
		pval.nested.boot = double(1),
		stats.nested = double(nalt.nested),
		stats.penal.nested = double(nalt.nested),
		S.nested = integer(1),
		stat.all = double(1),
		pval.all.asympt = double(1),
		pval.all.perm = double(1),
		pval.all.boot = double(1),
		stats.all = double(nalt.all),
		stats.penal.all = double(nalt.all),
		S.all = integer(1),
		all.subsets = matrix(integer(1),nalt.all,d), # matrix(integer(1),2^d-1,d),
		PACKAGE = "surv2sample")
		
	if (data.driven) {
		d.star = max(d0,1) # dimension of the smallest set(s)
		if (subsets=="nested") {
			if ((approx=="asympt")&&(d.star==1)) approx="twoterm"
			out = list(stat = temp$stat.nested,
				pval = temp[[paste("pval.nested",approx,sep=".")]][1],
				pval.asympt = temp$pval.nested.asympt,
				data.driven = data.driven,
				subsets.class = "nested",
				approx = approx,
				d = d,
				d0 = d0,
				d.star = d.star,
				df = d.star,
				S.ind = temp$S.nested+1, # index of the selected set (+1 because C returns 0,1,...)
				S.dim = d.star+temp$S.nested, # the smallest set (with 1 function) has S=0 in C
				S.set = d.star:(d.star+temp$S.nested),
				stats = temp$stats.nested,
				stats.penal = temp$stats.penal.nested)
			if ((approx=="perm")||(approx=="boot")) out$nsim = nsim
			if (d.star==1) out$pval.twoterm = temp$pval.nested.twoterm
		} else { # subset=="all"
			S.ind = temp$S.all+1 # index of the selected set (+1 because C returns 0,1,...)
			S.set = (1:d)[as.logical(temp$all.subsets[S.ind,])]
			out = list(stat = temp$stat.all,
				pval = temp[[paste("pval.all",approx,sep=".")]][1],
				pval.asympt = temp$pval.all.asympt,
				data.driven = data.driven,
				subsets.class = "all",
				approx = approx,
				d = d,
				d0 = d0,
				d.star = d.star,
				S.ind = S.ind,
				S.dim = length(S.set),
				S.set = S.set,
				stats = temp$stats.all,
				stats.penal = temp$stats.penal.all,
				all.subsets = temp$all.subsets)
			if (d0>0)
				out$df = d0 # for d0>0 it's chisq with df=d0
			else
				out$df = NA # for d0=0 it's max of dependent chisq_1, hence no df
			if ((approx=="perm")||(approx=="boot")||(d0==0)) # otherwise no simulations performed
				out$nsim = nsim
		}
	} else {
		out = list(stat = temp$stat.d,
			pval = temp[[paste("pval.d",aprx,sep=".")]][1],
			pval.asympt = temp$pval.d.asympt,
			approx = aprx,
			data.driven = data.driven,
			d = d,
			df = d)
			if ((approx=="perm")||(approx=="boot")) out$nsim = nsim
	}
	
	out$basis = basis
	out$time.transf = time.transf
		
	class(out) = c("surv2.neyman","neyman.test")
	out
}

print.surv2.neyman = function(x,...)
{
	if (!inherits(x,"surv2.neyman"))
		stop("must be an object of class 'surv2.neyman'")
	
	cat("\nTwo-sample test of equality of survival distributions")
	NextMethod("print",x,...)
	
	invisible()
}

summary.surv2.neyman = function(object,...)
{
	if (!inherits(object,"surv2.neyman"))
		stop("must be an object of class 'surv2.neyman'")
	
	cat("\nTwo-sample test of equality of survival distributions")
	NextMethod("summary",object,...)
	
	invisible()
}
