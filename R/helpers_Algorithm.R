###################################################################################

# NAME: ASAD HARIS
# DATE: 2015/07/10
# SUMMARY: SOME HELPER FUNCTIONS USED BY OUR SUPPORT REDUCTION ALGORITHM

###################################################################################

#A function to evaluate the log-likelihood when we do not have
#a full lambda vector
#dsz dz: are matrices of dim n*length(lambda).
#These are used to evaluate the log likelihhod
logLikeRestricted<- function(lambda,dsz, dz){
  #k<- length(z)
  sum(log(1-exp(-as.vector(lambda%*%dsz))) - as.vector(lambda%*%dz))

}

#A function to evaluate the gradient of the logLikelihood
#where lambda has reduced support
nablaLogLikeRestrict<- function(lambda, dsz, dz){
  #k<- length(z)
  const<- as.vector(1/(exp(lambda%*%dsz)-1))
  as.vector(dsz%*%const)-apply(dz,1,sum)

}

#A cpp implementation of tapply. Can be faster in
#some cases
tapply_fast<- function(x,index){
  unlist(lapply(split(x, fast_factor(index)), sum))
}

#The intermediate algorithm used in the support reduction algorihthm
#This is equivalent to MainAlg but with lower dimensions
InterAlg<- function(ini= 1, data,index, epsilon = 1e-3,maxiter = 500,
                    Dmatrix = NULL){
  a<- data$a
  b<- data$b
  t<- data$t
  z<- sort( unique(c(a,b,t)) )
  k<- length(z)
  n<- length(a)

  lambda<- rep(ini,length(index))

  if(is.null(Dmatrix)){
    #cat("Building D matrix...")
    Dmatrix<- cpp_buildMatrixD(t,a,b,z)
    #cat("done.\n")
  }

  #For this we will build the Dstar*(zj-zj-1) and D*(zj-zj-1)
  myf2<- stepfun(index, c(0,1:length(index)))
  fullIndex<- myf2(1:(k-1))
  diffz<- z[2:k]-z[1:(k-1)]
  diffz[which(diffz==Inf)]<- 1e+10

  #DMatrix
  Dstar<- Dmatrix$Dstar
  D<- Dmatrix$D

  #DsZ<- scale(Dstar, center = FALSE, scale = 1/diffz)
  DsZ<- cpp_scale(Dstar, diffz)
  dsz<- apply(DsZ, 1, function(x){tapply_fast(x, fullIndex)} )

  #DZ<- scale(D, center = FALSE, scale = 1/diffz)
  DZ<- cpp_scale(D, diffz)
  dz<- apply(DZ, 1, function(x){tapply_fast(x, fullIndex)} )

  #Index of lambda values which will be set at Infinity
  #indK<- which(apply(Dmatrix[[1]],2,sum)==0)
  for(i in 1:maxiter){
    #print(i)
    pi<- projLambda( lambda+ cpp_nablaLogLikeRestrict(lambda, dsz,dz ) ) - lambda
    pi[is.na(pi)]<- 1e+5

    newLambda <- cpp_newLambda(pi, lambda, dsz, dz)

    diff<- abs(cpp_logLikeRestricted(newLambda, dsz,dz) -
                 cpp_logLikeRestricted(lambda, dsz, dz))

    if( diff <= epsilon ){
      return(list(lambda = newLambda, conv = TRUE,dsz = dsz, dz = dz,
                  index = index, Dmatrix = Dmatrix, z = z))
    }else{
      lambda<- newLambda
    }
  }

  return(list(lambda = newLambda, conv = FALSE, dsz = dsz, dz = dz,
              index = index, Dmatrix = Dmatrix, z = z))
}


#This function evaluates the logLikelihood of an estimated hazard
#obj: an object given by the outpit of the function InterAlg
getlogLike<- function(obj){
  cpp_logLikeRestricted(obj$lam, obj$dsz, obj$dz)
}

#This function checks the KKT conditions for our algorithm
#by building a vector of dual variables which satisfy the
#KKT conditions. Some of the dual variables will be negative and hence
#we will later add those negative points to our extended support.
check.KKT<- function(obj){
  lam<- obj$lambda
  index<- obj$index
  Dmatrix<- obj$Dmatrix
  z<- obj$z

  k<- length(z)
  myf<- stepfun(index,c(0,lam) )
  LAMBDA<- myf(1:(k-1))
  h_1<- -1*LAMBDA[1]
  h_i<-  LAMBDA[1:(k-2)]-LAMBDA[2:(k-1)]
  h_i<- c(h_1,h_i)

  #Find gradient
  gradLogL<- cpp_nablaLogLike(LAMBDA, Dmatrix,z)
  gradHi<- matrix(0, ncol = k-1, nrow = k-1)
  diag(gradHi)<- -1
  diag(gradHi[-1,-ncol(gradHi)])<- 1

  #Find non-zero h_i functions
  nZeroHi<- which(h_i!=0)
  u_i<- rep(1,k-1)
  u_i[nZeroHi]<- 0

  a<- t(gradHi[-nZeroHi, -nZeroHi])
  b<- gradLogL[-nZeroHi]

  myu<- cpp_solve(b,a)
  u_i[-nZeroHi]<- myu
  return(u_i)
}

#This function returns a similar vector as the check.KKT function
#In this case we will need to look for support points which are positive.
check.derv<- function(obj){
  lam<- obj$lambda
  index<- obj$index
  Dmatrix<- obj$Dmatrix
  z<- obj$z

  k<- length(z)
  myf<- stepfun(index,c(0,lam) )
  LAMBDA<- myf(1:(k-1))

  #Find gradient
  derv<- nablaLogLike(LAMBDA, Dmatrix,z)
  vec<- rev(cumsum(rev(derv)))
}

#This function finds local maxima and returns an index vector
#for the local max which are ALSO POSITIVE. This is a helper function
#for support reduction used with the check.derv and check.KKT functions
#to add support points.
find.local.max<- function(vec, tol = 1e-3){
  index<- which(diff(sign(diff(vec)))== - 2) +1
  values<- vec[index]
  index<- index[(values>tol)]
  index
}
