### Kolmogorov--Smirnov two-sample test of equality of survival distributions

surv2.ks = function(x, group, process="w", approx="lwy", nsim=2000, nsim.plot=50)
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
	if ((!is.numeric(nsim))||(nsim<=0))
		stop("incorrect nsim (must be >0)")
	if ((!is.numeric(nsim.plot))||(nsim<0)||(nsim.plot>nsim))
		stop("incorrect nsim.plot (must be 0<=nsim.plot<=nsim)")
	aprx = tolower(approx)
	if (is.na(match(aprx,c("lwy","mart","perm","boot"))))
		stop("unknown approximation method (possible values 'lwy' or 'mart', 'perm', 'boot')")
	proc = tolower(process)
	if ((proc!="w")&&(proc!="b"))
		stop("unknown process type (possible values 'w', 'b')")
	
	tim = x[,1]
	r = order(tim)
	tim = tim[r]
	evt = x[r,2]
	grp = group[r]
	
	if ((aprx=="lwy")||(aprx=="mart")) {
		nlwy = nsim
		nlwy.plot = nsim.plot
		nperm = nboot = 0
		nperm.plot = nboot.plot = 0
	}
	if (aprx=="perm") {
		nperm = nsim
		nperm.plot = nsim.plot
		nlwy = nboot = 0
		nlwy.plot = nboot.plot = 0
	}
	if (aprx=="boot") {
		nboot = nsim
		nboot.plot = nsim.plot
		nperm = nlwy = 0
		nperm.plot = nlwy.plot = 0
	}
	
	temp = .C(".../src/twosample_ks_cm_ad",
		as.double(tim),
		as.integer(evt),
		as.integer(grp),
		as.integer(n),
		as.double(0),
		as.double(0),
		as.integer(nlwy),
		as.integer(nperm),
		as.integer(nboot),
		test.process = double(n), # increments of the test process ## pripadne se prida jeste nasim. trajekt. pro kresleni
		as.integer(nlwy.plot),
		as.integer(nperm.plot),
		as.integer(nboot.plot),
		test.process.lwy.plot = matrix(double(1),n,max(1,nlwy.plot)),
		test.process.perm.plot = matrix(double(1),n,max(1,nperm.plot)),
		test.process.boot.plot = matrix(double(1),n,max(1,nboot.plot)),
		stat.ks.w = double(1),
		pval.ks.w.lwy = double(1),
		pval.ks.w.perm = double(1),
		pval.ks.w.boot = double(1),
		pval.ks.w.asympt = double(1),
		stat.ks.b = double(1),
		pval.ks.b.lwy = double(1),
		pval.ks.b.perm = double(1),
		pval.ks.b.boot = double(1),
		pval.ks.b.asympt = double(1),
		stat.cm.w = double(1),
		pval.cm.w.lwy = double(1),
		pval.cm.w.perm = double(1),
		pval.cm.w.boot = double(1),
		stat.cm.b = double(1),
		pval.cm.b.lwy = double(1),
		pval.cm.b.perm = double(1),
		pval.cm.b.boot = double(1),
		stat.ad.w = double(1),
		pval.ad.w.lwy = double(1),
		pval.ad.w.perm = double(1),
		pval.ad.w.boot = double(1),
		stat.ad.b = double(1),
		pval.ad.b.lwy = double(1),
		pval.ad.b.perm = double(1),
		pval.ad.b.boot = double(1),
		PACKAGE = "surv2sample")
	
	if ((aprx=="lwy")||(aprx=="mart")) {
		test.process.sim = temp$test.process.lwy.plot
		if (proc=="w") {
			pval.ks = temp$pval.ks.w.lwy
			pval.ks.asympt = temp$pval.ks.w.asympt
			pval.cm = temp$pval.cm.w.lwy
			pval.ad = temp$pval.ad.w.lwy
		} else {
			pval.ks = temp$pval.ks.b.lwy
			pval.ks.asympt = temp$pval.ks.b.asympt
			pval.cm = temp$pval.cm.b.lwy
			pval.ad = temp$pval.ad.b.lwy
		}
	}
	if (aprx=="perm") {
		test.process.sim = temp$test.process.perm.plot
		if (proc=="w") {
			pval.ks = temp$pval.ks.w.perm
			pval.ks.asympt = temp$pval.ks.w.asympt
			pval.cm = temp$pval.cm.w.perm
			pval.ad = temp$pval.ad.w.perm
		} else {
			pval.ks = temp$pval.ks.b.perm
			pval.ks.asympt = temp$pval.ks.b.asympt
			pval.cm = temp$pval.cm.b.perm
			pval.ad = temp$pval.ad.b.perm
		}
	}
	if (aprx=="boot") {
		test.process.sim = temp$test.process.boot.plot
		if (proc=="w") {
			pval.ks = temp$pval.ks.w.boot
			pval.ks.asympt = temp$pval.ks.w.asympt
			pval.cm = temp$pval.cm.w.boot
			pval.ad = temp$pval.ad.w.boot
		} else {
			pval.ks = temp$pval.ks.b.boot
			pval.ks.asympt = temp$pval.ks.b.asympt
			pval.cm = temp$pval.cm.b.boot
			pval.ad = temp$pval.ad.b.boot
		}
	}
	test.process.sim = apply(test.process.sim,2,cumsum)
	
	if (proc=="w") {
		stat.ks = temp$stat.ks.w
		stat.cm = temp$stat.cm.w
		stat.ad = temp$stat.ad.w
	} else {
		stat.ks = temp$stat.ks.b
		stat.cm = temp$stat.cm.b
		stat.ad = temp$stat.ad.b
	}
	
	out = list(stat.ks = stat.ks,
		pval.ks = pval.ks,
		pval.ks.asympt = pval.ks.asympt,
		stat.cm = stat.cm,
		pval.cm = pval.cm,
		stat.ad = stat.ad,
		pval.ad = pval.ad,
		process = proc,
		approx = aprx,
		time = tim,
		test.process = cumsum(temp$test.process),
		test.process.sim = test.process.sim,
		nsim = nsim,
		nsim.plot = nsim.plot)
	
	class(out) = c("surv2.ks","lwy.test")
	out
}

print.surv2.ks = function(x,...)
{
	if (!inherits(x,"surv2.ks"))
		stop("must be an object of class 'surv2.ks'")
	
	cat("\nTwo-sample tests of equality of survival distributions",sep="")
	cat("\nKolmogorov--Smirnov statistic: ",x$stat.ks,",   p-value: ",format.pval(x$pval.ks),sep="")
	cat("\nCramer--von Mises statistic: ",x$stat.cm,",   p-value: ",format.pval(x$pval.cm),sep="")
	cat("\nAnderson--Darling statistic: ",x$stat.ad,",   p-value: ",format.pval(x$pval.ad),sep="")
	aprx = if ((x$approx=="lwy")||(x$approx=="mart")) "simulations"
		else if (x$approx=="perm") "permutations"
			else "bootstrap samples"
	cat("\np-values based on ",x$nsim," ",aprx,sep="")
	cat("\n\n")
	
	invisible()
}
