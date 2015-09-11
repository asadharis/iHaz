###################################################################################

# NAME: ASAD HARIS
# DATE: 2015/07/10
# SUMMARY: A SUPPORT REDUCTION ALGORITHM FOR ESTIMATING HAZARD FOR
#           LEFT TRUNCATED INTERVAL CENCSORED DATA.
#
#           THIS FILE CONTAINS THE S3 METHODS FOR FITTED OBJECTS

###################################################################################
#This file contains the main algorithm which aims to find the hazard function
#using a support reduction algorithm
#
#data: A list of three vectors labelled a, b and t. 't' is left truncation time
#       and (a,b) is the censoring interval
#ini.index: The inital index vector for support reduction algorithm. This
#           needs to be a subvector of 1:k where k is the length of the parameter vector.
#           This vector must contain 1.
#inter.maxiter: The maximum iterations for each intermediate optimization problem we do with
#             our reduced support.
#inter.tol: The tolerance for the intermediate optimization algorithms.
#main.maxiter: The maximum number of iterations in the main support reduction
#             algorithm
#check.condition: A string which specifies which stopping criteria to use.
#                 Two options are KKT(Karush Kan Tucker) conditions and
#                 derv (based on lemma 3.1 Wellner and Zhan 1997).

iHaz<- function(data, ini.index = 1:3, inter.maxiter = 1000, inter.tol = 1e-4,
                             main.maxiter = 1000, main.tol = 1e-3,
                            check.condition = c("derv","KKT"), verbose = FALSE){

  #Initalize initial index. This could just be the scaler 1.
  index<- ini.index

  #Begin building the D matrix.
  #This operation can be slow
  if(verbose)cat("Building D matrix...")
  a<- data$a
  t<- data$t
  b<- data$b
  z<- sort( unique(c(a,b,t)) )
  Dmatrix<- cpp_buildMatrixD(t,a, b, z)
  if(verbose)cat("done.\n\n")

  for(i in 1:main.maxiter){
    if(verbose)cat(paste0("Iteration Number: ", i, "\n"))

    #Find the MLE with a reduced support
    myans<- InterAlg(ini= 1, data, index = index, epsilon = inter.tol,
                     maxiter = inter.maxiter, Dmatrix = Dmatrix)

    #Check the conditions to find which support points to add
    if(check.condition[1]=="derv"){
      temp<- check.derv(myans)
      #This function finds the local maximum which are positive
      newindex<- cpp_findLocalMax(temp, main.tol)
    }else{
      u_i<- check.KKT(myans)
      #We need negative local minima hence we use -u_i
      newindex<- cpp_findLocalMax(-u_i, main.tol)
    }
    #Update the index set
    index2<- sort(unique(c(index, newindex)))

    #Check for convergence. Either there will be no new index
    #or the new index will already be in our index set. In which
    #case we'll continue on an infinite loop.
    if(length(index2)-length(index)==0){
      #Return object
      k<- length(myans$z)
      myf<- stepfun(z[myans$index],c(0,myans$lam) )
      LAMBDA<- myf(z[1:(k-1)])

      object<- list()
      object$hazard<- myf
      survivalF<- function(x){ sapply(x,EstimateSurv,
                                            lambda = LAMBDA, z = myans$z) }

      class(survivalF)<- "SurvivalFunction"
      object$survival<- survivalF
      object$index<- myans$index
      object$a<- a
      object$b<- b
      object$t<- t
      object$conv<- myans$conv
      object$call<- match.call()
      class(object)<- "iHaz"
      return(object)
    }else{
      #Else update index
      index<- index2
    }

  }
  #Return object without convergence
  k<- length(myans$z)
  myf<- stepfun(z[myans$index],c(0,myans$lam) )
  LAMBDA<- myf(z[1:(k-1)])
  object<- list()
  object$hazard<- myf
  survivalF<- function(x){ sapply(x,EstimateSurv,
                                  lambda = LAMBDA, z = myans$z) }

  class(survivalF)<- "SurvivalFunction"
  object$survival<- survivalF
  object$index<- myans$index
  object$a<- a
  object$b<- b
  object$t<- t
  object$call<- match.call()
  object$conv<- FALSE
  class(object)<- "iHaz"
  return(object)
}

#A print function for the objects of class. S3 generic method.
print.iHaz<- function(x,...){
  obj<- x
  if(obj$conv){
    cat("Call:\n")
    print(obj$call)
    cat("\nThe algorithm converged giving us MLEs for the hazard and survival functions.\n")
    cat("\nUsage for the hazard: $hazard(x) for a vector x." )
    cat("\nUsage for the survival function: $surv(x) for a vector x." )
  }else{
    cat("Call:\n")
    print(obj$call)
    cat("\nThe algorithm did not converge for either the main loop or the secondary loop.")
    cat("\nIncrease either 'inter.maxiter' or 'main.maxiter' (or both).")
    cat("\nCurrent estimates can still be accessed but may not be useful.\n")
    cat("\nUsage for the hazard: $hazard(x) for a vector x." )
    cat("\nUsage for the survival function: $surv(x) for a vector x." )
  }
}

plot.iHaz<- function(x,...){
  obj<- x
  z<- sort(unique( c(obj$a,obj$b,obj$t)))
  upper.lim<- max(z[is.finite(z)])
  lower.lim<- min(z[is.finite(z)])
  mys<- seq(lower.lim+1e-3, upper.lim-1e-3,length = 100)
  plot(mys, obj$hazard(mys), main = "Estimated hazard function for data",
       xlab = "x", ylab = expression(lambda*"(x)") ,...)
  par(ask = TRUE)
  plot(mys, obj$surv(mys), main = "Estimated survival function for data",
       xlab = "x", ylab = "S(x)" ,...)
  par(ask = FALSE)
}

print.SurvivalFunction<- function(x,...){
  cat("A survival function object.\n")
  cat("\nThe Maximum Likelihood Estimator of the survival function for our data.")
  cat("\nThis can be used as a regular R function with vector input.")
}


#This function estimats the survival curve at a scaler t
#Given an estimate for the hazard, lambda. This function evaluates the
#survival function at a scaler time point t.
EstimateSurv<- function(t,lambda,z){
  ind<- which(t<z)[1]
  k<- length(z)
  zdiff<- z[2:k]-z[1:(k-1)]
  if(ind==1){
    ans<- (t*lambda[ind])
  }else{
    ans<- sum(lambda[1:(ind-1)]*zdiff[1:(ind-1)])
    if(ind< length(z)){
      ans<- ans+  (t-z[ind-1])*lambda[ind]
    }
  }

  exp(-ans)
}


