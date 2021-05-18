#' solveols Function
#'
#' This function solves a linear system using Gauss-Seidel or Jacobi method, allows user to specify how many cores to use for parallel implementation.
#' @param method The method to solve a linear system, choose from "GaussSeidel","Jacobi"
#' @param cores Defaults to NULL, how many cores to use for parallel implementation.
#' @param A is n×n tridiagonal with alpha on the diagonals and -1 on theimmediate off diagonals.
#' @param v a vector repeated with 1 and 0.
#' @param iter the maximum iteration
#'
#'
#'
#' @return  return a list of error and runtime for each iteration
#'
#' @export
#' @examples
#' alpha3=3
#' n <- 100
#' A3 <- diag(alpha3, n)
#' A3[abs(row(A3) - col(A3)) == 1] <- -1
#' v=rep(c(1,0),n/2)
#' solveols("GaussSeidel",A=A3,v=v,iter=10)
#'
#'
solveols=function(method,cores=NULL,A,v,iter=10000){

if (method=="GaussSeidel"){
  res=GS.solve(A,v,iter)

}

else if (method=="Jacobi"){
  if(is.null(cores)){
    res=JacobiSeq.solve(A,v,iter)
  }
  else{
    res=JacobiPara.solve(A,v,iter,cores)
  }

}
else
  print("unsupported method type")

return(res)
}



GS.solve=function(A,v,iter=10000){
  n=dim(A)[1]

  L<-lower.tri(A)*A  # lower triang. A
  U<-upper.tri(A)*A  # upper triang. A
  D<-diag(diag(A))   # diag of A
  b=A%*%v

  LDinv=(MASS::ginv((L+D)))
  x0=rep(0,n)
  result.gs=matrix(NA,iter,2)

  i=1
  ptm0 <- proc.time()[3]
  while(i<=iter){
    #print(x0)
    # Gauss-Seidel formula
    x1<-LDinv%*%(b-U%*%x0)
    x0<-x1

    ptm=proc.time()[3] - ptm0
    result.gs[i,2]=as.numeric(ptm)
    error=norm(x1-v,"2")/norm(v,"2")
    #print(error)
    result.gs[i,1]=error
    i=i+1
  }

  return(result.gs)

}


JacobiSeq.solve=function(A,v,iter=10000){
  n=dim(A)[1]
  L<-lower.tri(A)*A  # lower triang. A
  U<-upper.tri(A)*A  # upper triang. A
  D<-diag(diag(A))   # diag of A
  Dinv=MASS::ginv(D)

  b=A%*%v
  x0=rep(0,n)
  x1=rep(0,n)
  result.gs=matrix(NA,iter,2)


  i=1
  ptm0 <- proc.time()[3]
  while(i<=iter){
    #print(x0)
    #Jacobi Sequential formula

    x1<-Dinv%*%(b-(L+U)%*%x0)
    x0<-x1
    ptm=proc.time()[3] - ptm0
    result.gs[i,2]=as.numeric(ptm)
    error=norm(x1-v,"2")/norm(v,"2")
    #print(error)
    result.gs[i,1]=error
    i=i+1
  }
  return(result.gs)
}

`%dopar%` <- foreach::`%dopar%`
JacobiPara.solve=function(A,v,iter=10000,cores=4){
  n=dim(A)[1]
  L<-lower.tri(A)*A  # lower triang. A
  U<-upper.tri(A)*A  # upper triang. A
  D<-diag(diag(A))   # diag of A
  LU=L+U
  Dinv=MASS::ginv(D)
  b=A%*%v
  x0=rep(0,n)
  result.gs=matrix(0,iter,2)

  cl=parallel::makeCluster(cores, setup_timeout = 0.5)
  doParallel::registerDoParallel(cl)

  i=1
  const=dim(A)[1]/cores
  ptm0 <- proc.time()[3]
  while(i<=iter){
    outlist = foreach::foreach(j=1:cores, .combine = cbind,.multicombine = TRUE) %dopar%{
      begin=const*(j-1)+1
      end=const*j
      x1<-Dinv[begin:end,begin:end]%*%(b[begin:end]-L[begin:end,]%*%x0-U[begin:end,]%*%x0)
      x1
    }

    x0=as.vector(as.numeric(outlist))
    ptm=proc.time()[3] - ptm0
    result.gs[i,2]=as.numeric(ptm)
    error=norm(x0-v,"2")/norm(v,"2")
    #print(error)
    result.gs[i,1]=error
    i=i+1
  }
  parallel::stopCluster(cl)
  return(result.gs)
}

