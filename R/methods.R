setValidity("corObject", function(object) {

	retval<-NULL

	#if ( (sum(dim(object@dat.miRNA))<3) || (is.null(object@dat.mRNA))) {
	#	retval<-c(retval,"Must define a miRNA and mRNA dataset")
	#comprovar els fenotips
	#	}

	#if (!all(rownames(object@pheno.miRNA)==colnames(object@miRNAdat))) {
	#	retval<-c(retval, "MiRNA samples should match")
	#	}

	#if (!all(rownames(object@pheno.mRNA)==colnames(object@mRNAdat))) {
	#	retval<-c(retval, "MRNA samples should match")
	#	}

	
	if(is.null(retval)) return(TRUE) 
	else return(retval)
}
)



setMethod("initialize", signature ="corObject", definition=function(.Object,...) {

	.Object<-callNextMethod()
	validObject(.Object)
	.Object
	}
)









