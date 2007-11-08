### Surv object for competing risks data

Survcomp = function (time, event)
{
	if (missing(event)) event = rep(1,length(time))
	if (length(time)!=length(event))
		stop("time and event have different lengths")
    if ( (!identical(sort(as.integer(unique(event))),min(event):max(event))) || (match(min(event),c(0,1),nomatch=0)==0) )
    	stop("incorrect values in event (must be integers 1,2,..., or 0 for censoring)")

	out = cbind(time,event)
	attr(out, "ncauses") = max(event)
	colnames(out) = c("time","event")
	class(out) = "Survcomp"
	out
}

is.Survcomp = function(x)
{
	inherits(x, "Survcomp")
}

print.Survcomp = function(x, quote=FALSE, ...)
{
	out = paste(format(x[,1])," (",ifelse(x[,2]==0,"+",x[,2]),")",sep="")
	print(out, quote=quote, ...)
	invisible(out)
}
