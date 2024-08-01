#' @name loglinearbart
#' @title loglinearbart
#' @useDynLib loglinearBART, .registration = TRUE
NULL
#' @param x.train Description of arg1
#' @param y.train Description of arg2
#' @param de.train Description of arg1
#' @param wi.train Description of arg2
#' @param v.train Description of arg1
#' @param x.test Description of arg2
#' @param power Description of arg2
#' @param base Description of arg1
#' @param ntree Description of arg1
#' @param ndpost Description of arg2
#' @param nskip Description of arg1
#' @param printevery Description of arg2
#' @param keepevery Description of arg1
#' @param keeptrainfits Description of arg2
#' @param usequants Description of arg1
#' @param numcut Description of arg2
#' @param printcutoffs Description of arg1
#' @param verbose Description of arg2
#' @param aa Description of arg2
#' @param bb Description of arg2

#' this is a function
#' @export
#'
loglinearbart = function(
   x.train, y.train,de.train,wi.train,v.train, x.test=matrix(0.0,0,0),power=2.0, base=.95,ntree=200,ndpost=1000, nskip=100,
   printevery=100, keepevery=1,keeptrainfits=TRUE,
   usequants=FALSE, numcut=100, printcutoffs=1,
   verbose=TRUE,aa=4.5,bb=4
)
{
   binary=FALSE

   if(is.factor(y.train)) {
      if(length(levels(y.train)) != 2) stop("y.train is a factor with number of levels != 2")
      binary = TRUE
      y.train = as.numeric(y.train)-1
   } else {
      if((length(unique(y.train)) == 2) & (max(y.train) == 1) & (min(y.train) == 0)) {
         cat('NOTE: assumming numeric response is binary\n')
         binary = TRUE
      }
   }

   if(is.vector(x.train) | is.factor(x.train)) x.train = data.frame(x=x.train)
   if(is.vector(x.test) | is.factor(x.test)) x.test = data.frame(x=x.test)

   if(is.data.frame(x.train)) {
      if(nrow(x.test)) {
         if(!is.data.frame(x.test)) stop('x.train is a data frame so x.test must be also')
	 xtemp = rbind(x.train,x.test)
	 xmtemp = makeind(xtemp)
	 x.train = xmtemp[1:nrow(x.train),,drop=FALSE]
	 x.test = xmtemp[nrow(x.train) + 1:nrow(x.test),,drop=FALSE]
      } else {
         x.train = makeind(x.train)
      }
   }

   #check input arguments:
   if((!is.matrix(x.train)) || (typeof(x.train)!="double")) stop("argument x.train must be a double matrix")
   if((!is.matrix(x.test)) || (typeof(x.test)!="double")) stop("argument x.test must be a double matrix")
   if(!binary) {
     if((!is.vector(y.train)) || (typeof(y.train)!="double")) stop("argument y.train must be a double vector")
   }
   if(nrow(x.train) != length(y.train)) stop("number of rows in x.train must equal length of y.train")
   if((nrow(x.test) >0) && (ncol(x.test)!=ncol(x.train))) stop("input x.test must have the same number of columns as x.train")
   if((mode(printevery)!="numeric") || (printevery<0)) stop("input printevery must be a positive number")
   if((mode(keepevery)!="numeric") || (keepevery<0)) stop("input keepevery must be a positive number")
   if((mode(ntree)!="numeric") || (ntree<0)) stop("input ntree must be a positive number")
   if((mode(ndpost)!="numeric") || (ndpost<0)) stop("input ndpost must be a positive number")
   if((mode(nskip)!="numeric") || (nskip<0)) stop("input nskip must be a positive number")
   if(mode(numcut)!="numeric") stop("input numcut must be a numeric vector")
   if(length(numcut)==1) numcut = rep(numcut,ncol(x.train))
   if(length(numcut) != ncol(x.train)) stop("length of numcut must equal number of columns of x.train")
   numcut = as.integer(numcut)
   if(min(numcut)<1) stop("numcut must be >= 1")
   if(typeof(usequants)  != "logical") stop("input usequants must a logical variable")
   if(typeof(keeptrainfits)  != "logical") stop("input keeptrainfits must a logical variable")
   if(typeof(verbose)  != "logical") stop("input verbose must a logical variable")
   if(mode(printcutoffs)  != "numeric") stop("input printcutoffs must be numeric")
   printcutoffs = as.integer(printcutoffs)
   if(printcutoffs <0) stop("input printcutoffs must be >=0")
   if(power <= 0) stop("power must be positive")
   if(base <= 0) stop("base must be positive")

   ncskip = floor(nskip/keepevery)
   ncpost = floor(ndpost/keepevery)
   nctot = ncskip + ncpost
   totnd = keepevery*nctot

   y=y.train
   de=de.train
   wi=wi.train
   v=v.train

   cres = .C('loglinearbart',as.integer(nrow(x.train)), as.integer(ncol(x.train)), as.integer(nrow(x.test)),
                   as.double(x.train), as.double(y),as.double(de),as.double(wi),as.double(v),
                   as.double(x.test), as.double(power), as.double(base),as.integer(ntree),as.integer(totnd),
             as.integer(printevery), as.integer(keepevery), as.integer(keeptrainfits),
             as.integer(numcut), as.integer(usequants), as.integer(printcutoffs),
             as.integer(verbose),as.double(aa),as.double(bb),
                   trdraw=double(nrow(x.train)*nctot),
                   tedraw=double(nrow(x.test)*nctot))


   yhat.train = yhat.test = yhat.train.mean = yhat.test.mean = NULL

  if (keeptrainfits) {
     yhat.train = matrix(cres$trdraw,nrow=nctot,byrow=T)[(ncskip+1):nctot,]
     #if(!is.null(dim(yhat.train))){yhat.train.mean = apply(yhat.train,2,mean)}
     #else{yhat.train.mean =mean(yhat.train)}
   }

   if (nrow(x.test)) {
     yhat.test = matrix(cres$tedraw,nrow=nctot,byrow=T)[(ncskip+1):nctot,]
     yhat.test.mean = apply(yhat.test,2,mean)
   }

      #yprint = matrix(cres$trdraw,nrow=nctot,byrow=T)[1,]
     retval = list(
       call=match.call(),
       yhat.train=yhat.train,
       #yhat.train.mean=yhat.train.mean,
       yhat.test=yhat.test,
       yhat.test.mean=yhat.test.mean,
       y = y.train
     )
   class(retval) = 'loglinearbart'
   return(invisible(retval))
}
