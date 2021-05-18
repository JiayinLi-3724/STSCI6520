#' elasticnet Function
#'
#' This function fits elastic net to data using coordinate descent algorithm.
#' @param X dataset X
#' @param Y dataset Y
#' @param alpha the parameter determining the ratio of L1 penalty to add
#' @param n sample size of data set
#' @param p dimension
#' @param nlambda the length of the sequence of the lambda
#' @param tol tolerance; when the difference less than tolerance, stop the algorithm
#'
#'
#'
#' @return  return a list of estimated coefficients for each lambda
#' @seealso \link[glmnet]{glmnet}
#' @export
#' @examples
#' p=20
#' n=20
#' mu <- rep(0,20)
#' Sigma <- diag(rep(1,20))
#' Sigma[1,2]=0.8
#' Sigma[2,1]=0.8
#' Sigma[5,6]=0.8
#' Sigma[6,5]=0.8
#' X=matrix(NA,n,p)
#' set.seed(65206520)
#' for (i in 1:n){X[i,]=MASS::mvrnorm(n=1, mu=mu, Sigma=Sigma)}
#' beta=c(2,0,-2,0,1,0,-1,0,rep(0,12))
#' set.seed(65206520)
#' eps=rnorm(n,0,1)
#' Y=X%*%beta+eps
#' elasticnet(X,Y,alpha=0,n,p,nlambda=100,tol=1e-5)


elasticnet=function(X,Y,alpha,n,p,nlambda=100,tol=1e-5){
  lambda.min.ratio=0.0005
  if (alpha==0){alpha_0=0.001}else{alpha_0=alpha}
  lambda.max=max( abs(t(Y - mean(Y)*(1-mean(Y))) %*% X ) )/ (alpha_0 * n)
  lambda.min=lambda.max*lambda.min.ratio
  log.lambda=seq(log(lambda.min),log(lambda.max),length.out = nlambda)
  lambda.seq=sort(exp(log.lambda),decreasing = TRUE)
  result=matrix(NA,p+1,nlambda)
  result[1,]= lambda.seq

  for (i in 1:nlambda){
    lambda= lambda.seq[i]
    beta_li=descent(X,Y,n,lambda,alpha,p,tol)
    result[2:(p+1),i]=beta_li
    #print(beta_li)
  }
  return(result)
}


descent=function(X,Y,n,lambda,alpha,p,tol=1e-5){
  beta_old=rep(0,p)
  beta_new=rep(0,p)
  run=TRUE

  while(run){
    beta_old=beta_new
    for (j in 1:p){
      X.j=X[,-j]
      beta.j=beta_new[-j]
      Ytilde.j=X.j%*%beta.j

      r.j=(Y- Ytilde.j)
      x.j=X[,j]
      z=(t(x.j)%*%r.j)/n
      if (alpha==0){
        S_operator=n*z
        beta_j= S_operator/(sum(x.j^2)+lambda)
      }
      else{
        S_operator=sign(z)*max(0, abs(z)-lambda*alpha)
        beta_j= S_operator/(sum(x.j^2)/n+lambda*(1-alpha))
      }

      #print(beta_j)
      beta_new[j]=beta_j
      #print(beta_new)
    }
    normbeta=sum(abs(beta_new-beta_old))
    #print(normbeta)
    #print(beta_new)

    if (normbeta < tol){run=FALSE}
  }
  return(beta_new)

}


