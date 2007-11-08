### Estimation of the two-sample proportional rate transformation model
### Based on the simplified partial likelihood (Bagdonavicius & Nikulin, 2000) (for PH it's Cox's PL)

proprate2 = function(x, group, model=0, beta.init=0, maxiter=20, eps=1e-09)
{
	if (!is.Surv(x))
		stop("x must be a Survcomp object")
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
	
	temp = .C("twosample_proprate_estim_r_call",
		as.integer(evt),
		as.integer(grp),
		as.integer(n),
		beta = as.double(beta.init),
		iter = as.integer(maxiter),
		as.double(eps),
		flag = integer(1),
		as.integer(model),
		q1 = double(n),
		q2 = double(n),
		q1dot = double(n),
		q2dot = double(n),
		score.process.incr = double(n),
		score = double(1),
		double(n),
		d11 = double(1),
		sigma11 = double(1),
		logliks = double(2),
		PACKAGE="surv2sample")
	
	if (temp$flag == 1000)
		warning("ran out of iterations and did not converge")
	
	n1 = sum(grp==1)
	n2 = n-n1
	Y1 = (n1-cumsum(grp==1)+1)
	Y2 = (n2-cumsum(grp==2)+1)
	q1 = temp$q1
	q2 = temp$q2
	dA1 = (evt*(grp==1))/Y1
	dA2 = (evt*(grp==2))/Y2
	G1 = cumsum(dA1/q1)
	G2 = cumsum(dA2/q2)
	G0 = cumsum(evt/(Y1*q1+exp(temp$beta)*Y2*q2))
	
	out = list(beta = temp$beta,
		var = temp$sigma11/(temp$d11)^2,
		model = model,
		model.name = ifelse(model==0,"proportional hazards","proportional odds"),
		iter = temp$iter,
		converged = (temp$flag!=1000),
		beta.init = beta.init,
		n = n,
		time = tim,
		d11 = temp$d11, # derivative of score
		sigma11 = temp$sigma11, # variance of score
		loglik.init = temp$logliks[1],
		loglik = temp$logliks[2],
		G1 = G1,
		G2 = G2,
		G0 = G0)
	
	class(out) = "proprate2.fit"
	out
}

print.proprate2.fit = function(x, digits=max(3,getOption("digits")-3),...)
{
	if (!inherits(x,"proprate2.fit"))
		stop("must be an object of class 'proprate2.fit'")
	
	beta = x$beta
	se = sqrt(x$var)
	temp = matrix(c(beta,exp(beta),se,beta/se,1-pchisq((beta/se)^2,1)), nrow=1)
	dimnames(temp) = list("",c("beta","exp(beta)","se(beta)","z","p"))
	
	cat("\nTwo-sample proportional rate model (",x$model.name,")",sep="")
	if (x$converged) {
		cat("\nn =",x$n)
		cat("\n")
		printCoefmat(temp)
		lr = 2*(x$loglik-x$loglik.init)
		cat("Likelihood ratio test of beta = ",x$beta.init,": statistic = ",lr,
			", p = ",format.pval(1-pchisq(lr,1)),sep="")
	} else {
		cat("\nEstimation did not converge")
	}
	cat("\n\n")
	
	invisible()
}

coef.proprate2.fit = function(object, ...)
{
	if (!inherits(fit,"proprate2.fit"))
		stop("must be an object of class 'proprate2.fit'")
	object$beta
}

plot.proprate2.fit = function(x, log.transform=FALSE, diff=FALSE,
	lwds=1, cols=1, ltys, ...)
{
	if (!inherits(x,"proprate2.fit"))
		stop("must be an object of class 'proprate2.fit'")
	
	if (!x$converged)
		stop("cannot plot, estimation did not converge")
	
	if (diff==TRUE) log.transform=TRUE
	
	if (missing(ltys)) {
		if ((log.transform==TRUE) && (diff==TRUE)) ltys=1 else ltys=c(1,1,2,2)
	}
	xx = c(0,x$time)
	yy = rbind(rep(0,4),cbind(x$G1,x$G2,x$G0,exp(x$beta)*x$G0))
	if (log.transform) yy = log(yy)
	if ((log.transform==TRUE) && (diff==TRUE)) yy = cbind(yy[,2]-yy[,1],yy[,4]-yy[,3])
	multiplot(xx,yy,lwds=lwds,cols=cols,ltys=ltys,...)
	invisible()
}
