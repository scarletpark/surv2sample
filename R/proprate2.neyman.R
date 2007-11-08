### Data-driven Neyman's smooth test of proportional rates in two samples
### with fixed dimension or data-driven with nested subsets

proprate2.neyman = function(x, group, model=0, data.driven=TRUE,
	d=ifelse(data.driven,5,3), basis="legendre", time.transf="F",
	beta.init=0, maxiter=20, eps=1e-09, choltol=1e-07)
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
	
	if (!(d>0)) stop("d must be >0")
	time.transf=toupper(time.transf)
	if (time.transf=="L") time.transf="A"
	tim.trnsf = match(time.transf,c("F","A","I"))
	if (is.na(tim.trnsf)) stop("unknown time transformation")
	basis = tolower(basis)
	if (basis=="cosine") basis = "cos"
	bas = match(tolower(basis),c("legendre","cos"))
	if (is.na(bas)) stop("unknown basis")
	
	n = nrow(x)
	if (n!=length(group)) stop("lengths of 'x' and 'group' differ")

	tim = x[,1]
	r = order(tim)
	tim = tim[r]
	evt = x[r,2]
	grp = group[r]
	
	d0=ifelse(data.driven,1,d)
	penalty = log(n)
	
	d.star = max(d0,1) # minimal dimension of subset
	
	temp = .C("twosample_proprate_neyman",
		as.double(tim),
		as.integer(evt),
		as.integer(grp),
		as.integer(n),
		as.integer(model),
		beta = as.double(beta.init),
		as.integer(maxiter),
		as.double(eps), # eps
		flag = integer(1),
		as.integer(d),
		as.integer(d0),
		as.integer(tim.trnsf),
		as.integer(bas),
		as.double(penalty),
		as.double(choltol),
		score = double(d),
		sigma = matrix(double(1),d,d),
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

	if (temp$flag == 1000)
		warning("ran out of iterations and did not converge")
	
	if (data.driven) {
		out = list(stat = temp$stat.nested,
			pval = temp$pval.nested.twoterm,
			pval.twoterm = temp$pval.nested.twoterm,
			pval.asympt = temp$pval.nested.asympt,
			basis = basis,
			time.transf = time.transf,
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
			stats.penal = temp$stats.nested.penal,
			model = model,
			model.name = ifelse(model==0,"proportional hazards","proportional odds"),
			converged = (temp$flag!=1000))
	} else {
		out = list(stat = temp$stat.nested,
			pval = temp$pval.d,
			pval.asympt = temp$pval.d,
			basis = basis,
			time.transf = time.transf,
			data.driven = data.driven,
			approx = "asympt",
			d = d,
			df = d,
			model = model,
			model.name = ifelse(model==0,"proportional hazards","proportional odds"),
			converged = (temp$flag!=1000))
	}
	
	class(out) = c("proprate2.neyman","neyman.test")
	out
}

print.proprate2.neyman = function(x,...)
{
	if (!inherits(x,"proprate2.neyman"))
		stop("must be an object of class 'proprate2.neyman'")
	
	cat("\nTwo-sample proportional rate model (",x$model.name,")",sep="")
	if (x$converged)
  		NextMethod("print",x,...)
	else {
		cat("\nEstimation did not converge")
		cat("\n\n")
	}
	
	invisible()
}

summary.proprate2.neyman = function(object,...)
{
	if (!inherits(object,"proprate2.neyman"))
		stop("must be an object of class 'proprate2.neyman'")
	
	cat("\nTwo-sample proportional rate model (",object$model.name,")",sep="")
	if (object$converged)
  		NextMethod("summary",object,...)
	else {
		cat("\nEstimation did not converge")
		cat("\n\n")
	}
	
	invisible()
}
