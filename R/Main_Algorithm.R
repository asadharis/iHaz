###################################################################################

# NAME: ASAD HARIS
# DATE: 2015/07/10
# SUMMARY: A SUPPORT REDUCTION ALGORITHM FOR ESTIMATING HAZARD FOR 
#           LEFT TRUNCATED INTERVAL CENCSORED DATA.

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

SupportReduction<- function(data, ini.index = 1:3, inter.maxiter = 1000, inter.tol = 1e-4, 
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
      return(myans)
    }else{
      #Else update index
      index<- index2
    }
    
  }
  #Return object without convergence
  myans$conv<- FALSE
  return(myans)
}
