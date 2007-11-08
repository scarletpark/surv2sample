### Kolmogorov--Smirnov test of fit of the two-sample proportional rate model

proprate2.ks = function(x, group, model=0, nsim=2000, nsim.plot=50,
	beta.init=0, maxiter=20, eps=1e-09)
{
	if (!is.Surv(x))
		stop("x must be a Surv object")
	if (attr(x,"type")!="right")
		stop("right censored data allowed only")
	if (missing(group))
		stop("group missing")
	if ( !identical(sort(as.integer(unique(group))),1:2) )
		stop("values of group must be 1,2")
	if ( (model!=0)&&(model!=1) )
		stop("model must be 0 (prop. haz.) or 1 (prop.odds)")
	## TODO: user defined models
	
	n = nrow(x)
	if (n!=length(group)) stop("lengths of 'x' and 'group' differ")

	tim = x[,1]
	r = order(tim)
	tim = tim[r]
	evt = x[r,2]
	grp = group[r]
	
	temp = .C("twosample_proprate_ks",
		as.integer(evt),
		as.integer(grp),
		as.integer(n),
		as.integer(model),
		beta = as.double(beta.init),
		iter = as.integer(maxiter),
		as.double(eps),
		flag = integer(1),
		as.integer(nsim),
		test.process = double(n),
		stat = double(1),
		pval.sim = double(1),
		test.process.sim = matrix(0,n,max(1,nsim.plot)),
		as.integer(nsim.plot),
		PACKAGE="surv2sample")
	
	if (temp$flag == 1000)
		warning("ran out of iterations and did not converge")
	
	out = list(stat = temp$stat,
		time = tim,
		test.process = cumsum(temp$test.process),
		nsim = nsim,
		pval.sim = temp$pval.sim,
		nsim.plot = nsim.plot,
		test.process.sim = apply(temp$test.process.sim,2,cumsum),
		model = model,
		model.name = ifelse(model==0,"proportional hazards","proportional odds"),
		converged = (temp$flag!=1000))
	
	class(out) = c("proprate2.ks","lwy.test")
	out
}


print.proprate2.ks = function(x,...)
{
	if (!inherits(x,"proprate2.ks"))
		stop("must be an object of class 'proprate2.ks'")
	
	cat("\nTwo-sample proportional rate model (",x$model.name,")",sep="")
	cat("\nKolmogorov--Smirnov test of fit")
	cat("\nStatistic: ",x$stat,",   p-value: ",format.pval(x$pval.sim),sep="")
	cat("\np-value based on",x$nsim,"simulations")
	cat("\n\n")
	
	invisible()
}

