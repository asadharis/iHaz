###################################################################################

# NAME: ASAD HARIS
# DATE: 2015/07/10
# SUMMARY: HELPER FUNCTIONS FOR ALGORITHMS USED

###################################################################################


#Function to build matrices for terms d_ij and d_ij*
#As done in Pan and Chappell 1998 the likelihood is expressed
#in terms of values d_{i,j} and d_{i,j}*. We evaluate these in matrix
#form by this function here
#t: The vector of left truncated times
#(a,b): vectors for the censoring intervals
#z: sort(unique(c(t,a,b)))
buildMatrixD<- function(t,a,b,z){
  k<- length(z)
  n<- length(t)
  Dstar<- D<- matrix(0, nrow = n, ncol = k-1)
  
  z<- c(z)
  
  for(i in 1:n){
    for(j in 2:(k)){
      D[i,j-1]<- 1*(z[j-1]>=t[i] & z[j] <= a[i])
      Dstar[i,j-1]<- 1*(z[j-1]>=a[i] & z[j] <= b[i])
    }
  }
  return(list(D = D, Dstar = Dstar))
}

#The projection function used in the algorithms for maximum likelihood
#In our case it is simply a non-negative isotonic regression
#We use the default function in R
projLambda<- function(lambda){
  whichInf<- which(lambda==Inf)
  
  c(pmax(isoreg(x = lambda[lambda!=Inf])$yf,0), rep(Inf,length(whichInf)))
  
}

#This function evaluates the logLikehood for a given lambda value
#Dmatrix: the output of the function buildMatrixD
#z: as defined in buildMatrixD function
logLike<- function(lambda, Dmatrix, z){
  
  D<- Dmatrix[[1]]
  Dstar<- Dmatrix[[2]]
  k<- length(z)
  
  z0<- z[1:(k-1)]
  #Obtain the vector lambda_z in report
  lambda_z<- lambda*(z[2:k]-z0)
  lambda_z[lambda_z==Inf]<- 1e+10
  
  #Return the sum A+B
  sum( log(1-exp(-Dstar%*%lambda_z) ) - D%*%lambda_z )
}

#This function evaluates the derivative of the log-Likelihood
#The parameters of the function are the same as that for the logLike function
nablaLogLike<- function(lambda,Dmatrix,z){
  D<- Dmatrix[[1]]
  Dstar<- Dmatrix[[2]]
  k<- length(z)
  
  z0<- z[1:(k-1)]
  zdiff<- z[2:k]-z0
  zdiff[zdiff==Inf]<- 1e+5
  #Obtain the vector lambda_z in report
  lambda_z<- lambda*zdiff
  lambda_z[lambda_z==Inf]<- 1e+5
  
  #We evaluate nablaB first
  nablaB<- apply(D,2,sum)
  nablaB<- nablaB*zdiff
  
  #We now evaluate nablaA
  nablaA<- as.vector(exp(Dstar%*%lambda_z)-1)
  
  myfunc<- function(x,nablaA){
    ans<- x/nablaA
    ans[which(is.na(ans))]<- 0
    ans
  }
  
  nablaA<- zdiff*apply(apply(Dstar,2,myfunc, nablaA),2,sum)
  
  #Return the vector of vector
  as.vector(nablaA-nablaB)
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
