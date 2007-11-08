### Gill--Schumacher test of fit of the two-sample proportional rate model

proprate2.gs = function(x, group, model=0, weight1="logrank", weight2="prentice")
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
	
	wtcodes = c("logrank"=0,"prentice"=1,"gehan"=-1)
	weight1 = tolower(weight1)
	wt1 = wtcodes[weight1]
	if (is.na(wt1)) stop("unknown weight1")
	weight2 = tolower(weight2)
	wt2 = wtcodes[weight2]
	if (is.na(wt2)) stop("unknown weight2")
	if (wt1==wt2) stop("weight1 and weight2 cannot be the same")
	
	n = nrow(x)
	if (n!=length(group)) stop("lengths of 'x' and 'group' differ")
	
	tim = x[,1]
	r = order(tim)
	tim = tim[r]
	evt = x[r,2]
	grp = group[r]

	temp = .C("twosample_proprate_gs",
		as.integer(evt),
		as.integer(grp),
		as.integer(n),
		as.integer(model),
		as.integer(wt1),
		as.integer(wt2),
		eta1 = double(1),
		eta2 = double(1),
		stat = double(1),
		pval = double(1),
		PACKAGE = "surv2sample")
	
	
	out = list(stat = temp$stat,
		pval = temp$pval,
		model = model,
		model.name = ifelse(model==0,"proportional hazards","proportional odds"),
		eta1 = temp$eta1,
		eta2 = temp$eta2,
		weight1 = weight1,
		weight2 = weight2)
	
	class(out) = "proprate2.gs"
	out
}

print.proprate2.gs = function(x,...)
{
	if (!inherits(x,"proprate2.gs"))
		stop("must be an object of class 'proprate2.gs'")
	
	cat("\nTwo-sample proportional rate model (",x$model.name,")",sep="")
	cat("\nGill--Schumacher test")
	cat("\neta1 = ",x$eta1," (",x$weight1," weight), eta2 = ",
		x$eta2," (",x$weight2," weight)",sep="")
	cat("\nStatistic: ",x$stat,",   p-value: ",format.pval(x$pval),sep="")
	cat("\n\n")
	
	invisible()
}

