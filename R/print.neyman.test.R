### Summary and print method for Neyman's tests (fixed or data-driven)

print.neyman.test = function(x, detail=FALSE,...)
{
	if (!inherits(x,"neyman.test"))
		stop("must be an object of class 'neyman.test'")
	
	if (x$data.driven) {
		cat("\nData-driven Neyman's smooth test")
		cat("\nBasis of smooth functions: ",x$basis,sep="")
		if (x$subsets.class=="nested") {
			cat("\nSubsets: ",x$subsets.class,", d = ",x$d,", d0 = ",x$d0,
				", selected dimension: ",x$S.dim,sep="")
			cat("\nTest statistic: ",x$stat,",   p-value: ",format.pval(x$pval)," ",sep="")
			if (x$approx=="twoterm")
				cat("(two-term approximation)",sep="")
			else 
				if (x$approx=="boot")
					cat("(",x$nsim," bootstrap samples)",sep="")
				else
					if (x$approx=="perm")
						cat("(",x$nsim," permutations)",sep="")
					else # "asympt"
						cat("(",x$df," df)",sep="")
			if (detail) {
				dims = x$d.star:x$d
				selected = rep("",length(dims))
				selected[x$S.ind] = "<="
				temp = cbind(dim=dims,round(cbind(stat=x$stats,penal.stat=x$stats.penal),4),selected=selected)
				rownames(temp)=rep("",length(dims))
				cat("\n")
				print(temp,quote=FALSE)
			} else
				cat("\n")
		} else { # all subsets
			cat("\nSubsets: ",x$subsets.class,", d = ",x$d,", d0 = ",x$d0,
				", selected set: ",paste(x$S.set,collapse=" "),sep="")
			cat("\nTest statistic: ",x$stat,",   p-value: ",format.pval(x$pval)," ",sep="")
			if (x$approx=="boot")
				cat("(",x$nsim," bootstrap samples)",sep="")
			else
				if (x$approx=="perm")
					cat("(",x$nsim," permutations)",sep="")
				else # "asympt"
					if (x$d0>0) # for d0==0 it's max of chisq1, no df
						cat("(",x$df," df)",sep="")
			if (detail) {
				nalt = length(x$stats)
				selected = rep("",nalt)
				selected[x$S.ind] = "<="
				star.mat = ifelse(x$all.subset==1,"*","")
				colnames(star.mat) = paste("phi_",1:x$d,sep="")
				temp = cbind(star.mat,round(cbind(stat=x$stats,penal.stat=x$stats.penal),4),selected=selected)
				rownames(temp)=rep("",nalt)
				cat("\n")
				print(temp,quote=FALSE)
			} else
				cat("\n")
		}
		cat("\n")
	} else {
		cat("\nNeyman's smooth test")
		cat("\nBasis of smooth functions: ",x$basis,", d = ",x$d,sep="")
		cat("\nTest statistic: ",x$stat,",   p-value: ",format.pval(x$pval)," ",sep="")
		if (x$approx=="twoterm")
			cat("(two-term approximation)",sep="")
		else 
			if (x$approx=="boot")
				cat("(",x$nsim," bootstrap samples)",sep="")
			else
				if (x$approx=="perm")
					cat("(",x$nsim," permutations)",sep="")
				else # "asympt"
					cat("(",x$df," df)",sep="")
		cat("\n\n")
	}
	
	invisible()
}

summary.neyman.test = function(object, ...)
{
	if (!inherits(object,"neyman.test"))
		stop("must be an object of class 'neyman.test'")
	
	print.neyman.test(object,detail=TRUE,...)

	invisible()
}

